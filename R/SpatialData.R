#' @include zzz.R
#' @include AnnDataEncoding.R
#' @include LoadZarr.R
#' @include SaveZarr.R
#' @include SpatialConversion.R
#' @importFrom Matrix sparseMatrix t
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject Images
#'   GetTissueCoordinates CreateFOV scalefactors DefaultAssay
#' @importFrom SeuratObject Cells CreateFOV AddMetaData
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SpatialData zarr format support
#
# SpatialData is the scverse standard for spatial omics data.
# It uses zarr stores with a specific layout combining OME-NGFF images
# and anndata tables.
#
# Layout:
#   spatialdata.zarr/
#   +-- .zattrs       {"spatialdata_attrs": {"version": "0.2.0"}}
#   +-- .zgroup       {"zarr_format": 2}
#   +-- images/       OME-NGFF multiscale images
#   +-- labels/       segmentation masks
#   +-- points/       point annotations (transcript coords)
#   +-- shapes/       geometric shapes (spots, cells)
#   +-- tables/       anndata expression tables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Read SpatialData
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read a SpatialData zarr store into a Seurat object
#'
#' Reads a SpatialData zarr store following the scverse SpatialData
#' specification. The expression data is read from the \code{tables/}
#' subdirectory (standard anndata zarr), spatial coordinates from
#' \code{shapes/} or \code{points/}, and tissue images from \code{images/}
#' in OME-NGFF format.
#'
#' @param path Path to the SpatialData .zarr directory
#' @param table Name of the table within \code{tables/} to read
#'   (default: "table")
#' @param images Logical; if TRUE (default), attempt to read OME-NGFF images
#'   from the \code{images/} directory
#' @param verbose Show progress messages
#'
#' @return A \code{\link[SeuratObject]{Seurat}} object with spatial
#'   coordinates and optionally tissue images
#'
#' @export
#'
readSpatialData <- function(path, table = "table", images = TRUE,
                            verbose = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required for readSpatialData. ",
         "Install with: install.packages('jsonlite')", call. = FALSE)
  }

  if (!dir.exists(path)) {
    stop("SpatialData store not found: ", path, call. = FALSE)
  }

  # Validate root zarr metadata (v2: .zgroup, v3: zarr.json)
  is_zarr_v2 <- file.exists(file.path(path, ".zgroup"))
  is_zarr_v3 <- file.exists(file.path(path, "zarr.json"))
  if (!is_zarr_v2 && !is_zarr_v3) {
    stop("Not a valid zarr store (missing .zgroup or zarr.json): ", path, call. = FALSE)
  }

  # Check for spatialdata_attrs
  zattrs_path <- file.path(path, ".zattrs")
  if (!file.exists(zattrs_path) && is_zarr_v3) {
    # zarr v3: attrs may be in zarr.json
    zj <- .zarr_read_json(file.path(path, "zarr.json"))
    root_attrs <- zj$attributes %||% list()
  } else {
    root_attrs <- .zarr_read_json(zattrs_path)
  }
  sd_version <- root_attrs$spatialdata_attrs$version
  if (verbose) {
    if (!is.null(sd_version)) {
      message("Reading SpatialData store (v", sd_version, "): ", path)
    } else {
      message("Reading SpatialData store: ", path)
    }
  }

  # ---------------------------------------------------------------------------
  # 1. Read the anndata table
  # ---------------------------------------------------------------------------
  table_path <- file.path(path, "tables", table)
  if (!dir.exists(table_path)) {
    # Try listing available tables
    tables_dir <- file.path(path, "tables")
    if (dir.exists(tables_dir)) {
      avail <- list.dirs(tables_dir, recursive = FALSE, full.names = FALSE)
      stop("Table '", table, "' not found in ", tables_dir,
           ". Available: ", paste(avail, collapse = ", "), call. = FALSE)
    }
    stop("No tables directory found in SpatialData store: ", path,
         call. = FALSE)
  }

  if (verbose) message("Reading expression table: tables/", table)
  obj <- readZarr(table_path, verbose = verbose)

  # ---------------------------------------------------------------------------
  # 2. Parse spatialdata_attrs from the table
  # ---------------------------------------------------------------------------
  table_attrs <- .zarr_read_json(file.path(table_path, ".zattrs"))
  sd_attrs <- table_attrs$spatialdata_attrs
  region <- NULL
  region_key <- "region"
  instance_key <- "instance_id"


  if (!is.null(sd_attrs)) {
    region <- sd_attrs$region
    if (is.list(region)) region <- unlist(region)
    region_key <- sd_attrs$region_key %||% "region"
    instance_key <- sd_attrs$instance_key %||% "instance_id"
    if (verbose && !is.null(region)) {
      message("SpatialData region(s): ", paste(region, collapse = ", "))
      message("  region_key: ", region_key, ", instance_key: ", instance_key)
    }
  }

  # Store spatialdata attributes in misc for roundtrip
  obj@misc[["__spatialdata_attrs__"]] <- sd_attrs
  obj@misc[["__spatialdata_version__"]] <- sd_version

  # ---------------------------------------------------------------------------
  # 3. Read spatial coordinates from shapes/ or points/
  # ---------------------------------------------------------------------------
  cell_names <- Cells(obj)
  coords_added <- FALSE

  if (!is.null(region) && length(region) > 0) {
    for (region_name in region) {
      # Try shapes first
      shapes_path <- file.path(path, "shapes", region_name)
      if (dir.exists(shapes_path)) {
        if (verbose) message("Reading shapes: shapes/", region_name)
        shape_data <- .read_spatialdata_shapes(shapes_path, verbose = verbose)
        if (!is.null(shape_data)) {
          obj <- .add_spatialdata_coords(
            obj, shape_data$coords, shape_data$radius,
            region_name = region_name, cell_names = cell_names,
            instance_key = instance_key, verbose = verbose
          )
          coords_added <- TRUE
          next
        }
      }

      # Try points (only if shapes not found for this region)
      points_path <- file.path(path, "points", region_name)
      if (dir.exists(points_path)) {
        if (verbose) message("Reading points: points/", region_name)
        point_data <- .read_spatialdata_points(points_path, verbose = verbose)
        if (!is.null(point_data)) {
          # Points are typically transcript-level; store in misc
          obj@misc[[paste0("points_", region_name)]] <- point_data
          if (verbose) message("  Stored ", nrow(point_data),
                               " points in misc$points_", region_name)
        }
      }
    }
  }

  # If no region specified, scan shapes/ and points/ directories
  if (!coords_added) {
    shapes_dir <- file.path(path, "shapes")
    if (dir.exists(shapes_dir)) {
      shape_names <- list.dirs(shapes_dir, recursive = FALSE,
                               full.names = FALSE)
      for (sn in shape_names) {
        if (verbose) message("Reading shapes: shapes/", sn)
        shape_data <- .read_spatialdata_shapes(file.path(shapes_dir, sn),
                                               verbose = verbose)
        if (!is.null(shape_data)) {
          obj <- .add_spatialdata_coords(
            obj, shape_data$coords, shape_data$radius,
            region_name = sn, cell_names = cell_names,
            instance_key = instance_key, verbose = verbose
          )
          coords_added <- TRUE
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 4. Read OME-NGFF images
  # ---------------------------------------------------------------------------
  if (images) {
    images_dir <- file.path(path, "images")
    if (dir.exists(images_dir)) {
      image_names <- list.dirs(images_dir, recursive = FALSE,
                               full.names = FALSE)
      for (img_name in image_names) {
        if (verbose) message("Reading image: images/", img_name)
        img_data <- .read_ome_ngff_image(file.path(images_dir, img_name),
                                         verbose = verbose)
        if (!is.null(img_data)) {
          # Store in misc for now; if spatial coords present, try to
          # build a proper Seurat image object
          if (coords_added) {
            obj <- .attach_spatialdata_image(
              obj, img_data, img_name = img_name, verbose = verbose
            )
          } else {
            obj@misc[[paste0("image_", img_name)]] <- img_data
            if (verbose) message("  Stored image in misc$image_", img_name)
          }
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 5. Read labels (segmentation masks) metadata
  # ---------------------------------------------------------------------------
  labels_dir <- file.path(path, "labels")
  if (dir.exists(labels_dir)) {
    label_names <- list.dirs(labels_dir, recursive = FALSE,
                             full.names = FALSE)
    if (length(label_names) > 0) {
      obj@misc[["__spatialdata_labels__"]] <- label_names
      if (verbose) message("Found label layers: ",
                           paste(label_names, collapse = ", "))
    }
  }

  if (verbose) {
    message("\nSuccessfully loaded SpatialData store")
    message("  Cells: ", ncol(obj))
    message("  Features: ", nrow(obj))
    if (length(Images(obj)) > 0) {
      message("  Images: ", paste(Images(obj), collapse = ", "))
    }
  }

  obj
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Write SpatialData
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write a Seurat object to SpatialData zarr format
#'
#' Writes a Seurat object to the scverse SpatialData zarr specification.
#' The expression data is written as an anndata table to \code{tables/table/},
#' spatial coordinates are written to \code{shapes/}, and tissue images
#' are written to \code{images/} in OME-NGFF format.
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param path Path for the output .zarr directory
#' @param table Name for the table (default: "table")
#' @param region Name for the spatial region (default: auto-detect from
#'   image names, or "spots")
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{path}
#'
#' @export
#'
writeSpatialData <- function(object, path, table = "table",
                             region = NULL, overwrite = FALSE,
                             verbose = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required for writeSpatialData. ",
         "Install with: install.packages('jsonlite')", call. = FALSE)
  }

  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }

  if (dir.exists(path) || file.exists(path)) {
    if (!overwrite) {
      stop("Output exists: ", path,
           ". Use overwrite = TRUE to replace.", call. = FALSE)
    }
    unlink(path, recursive = TRUE)
  }

  # Determine region name
  image_names <- Images(object)
  if (is.null(region)) {
    if (length(image_names) > 0) {
      region <- image_names[1]
    } else {
      region <- "spots"
    }
  }

  if (verbose) message("Writing SpatialData store: ", path)

  # ---------------------------------------------------------------------------
  # 1. Root group with spatialdata_attrs
  # ---------------------------------------------------------------------------
  sd_version <- object@misc[["__spatialdata_version__"]] %||% "0.2.0"

  .zarr_create_group(path, attrs = list(
    spatialdata_attrs = list(version = sd_version)
  ))

  # ---------------------------------------------------------------------------
  # 2. Write the anndata table
  # ---------------------------------------------------------------------------
  table_path <- file.path(path, "tables", table)
  dir.create(file.path(path, "tables"), recursive = TRUE,
             showWarnings = FALSE)

  if (verbose) message("Writing expression table: tables/", table)
  writeZarr(object, table_path, verbose = verbose)

  # Add spatialdata_attrs to the table
  sd_attrs <- object@misc[["__spatialdata_attrs__"]]
  if (is.null(sd_attrs)) {
    sd_attrs <- list(
      region = region,
      region_key = "region",
      instance_key = "instance_id"
    )
  }

  # Merge spatialdata_attrs into existing table .zattrs
  table_zattrs_path <- file.path(table_path, ".zattrs")
  table_attrs <- .zarr_read_json(table_zattrs_path)
  table_attrs$spatialdata_attrs <- sd_attrs
  json <- jsonlite::toJSON(table_attrs, auto_unbox = TRUE, null = "null",
                           pretty = TRUE)
  writeLines(json, table_zattrs_path)

  # ---------------------------------------------------------------------------
  # 3. Write shapes (spatial coordinates)
  # ---------------------------------------------------------------------------
  has_spatial <- FALSE

  if (length(image_names) > 0) {
    for (img_id in image_names) {
      img_obj <- object[[img_id]]
      coords <- tryCatch(GetTissueCoordinates(img_obj), error = function(e) NULL)
      if (is.null(coords)) next

      # Determine x, y columns
      if (all(c("x", "y") %in% colnames(coords))) {
        xy <- as.matrix(coords[, c("x", "y")])
      } else if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
        xy <- as.matrix(coords[, c("imagecol", "imagerow")])
      } else {
        num_cols <- vapply(coords, is.numeric, logical(1))
        if (sum(num_cols) >= 2) {
          xy <- as.matrix(coords[, which(num_cols)[1:2]])
        } else {
          next
        }
      }
      colnames(xy) <- c("x", "y")

      # Determine radius
      radius <- NULL
      sf <- tryCatch(GetScaleFactors(img_obj), error = function(e) NULL)
      if (!is.null(sf)) {
        spot_diam <- sf[["spot"]] %||% sf[["spot_diameter_fullres"]]
        if (!is.null(spot_diam) && is.finite(as.numeric(spot_diam))) {
          radius <- as.numeric(spot_diam) / 2.0
        }
      }
      if (is.null(radius)) radius <- 1.0

      shapes_name <- if (img_id == image_names[1]) region else img_id
      shapes_path <- file.path(path, "shapes", shapes_name)
      if (verbose) message("Writing shapes: shapes/", shapes_name)
      .write_spatialdata_shapes(xy, radius = radius, path = shapes_path,
                                verbose = verbose)
      has_spatial <- TRUE
    }
  }

  # Fallback: check meta.data for spatial_x / spatial_y
  if (!has_spatial) {
    if (all(c("spatial_x", "spatial_y") %in% colnames(object@meta.data))) {
      xy <- as.matrix(object@meta.data[, c("spatial_x", "spatial_y")])
      colnames(xy) <- c("x", "y")
      shapes_path <- file.path(path, "shapes", region)
      if (verbose) message("Writing shapes from metadata: shapes/", region)
      .write_spatialdata_shapes(xy, radius = 1.0, path = shapes_path,
                                verbose = verbose)
      has_spatial <- TRUE
    }
  }

  if (!has_spatial) {
    # Create empty shapes group
    .zarr_create_group(file.path(path, "shapes"))
    if (verbose) message("No spatial coordinates found; shapes/ is empty")
  }

  # ---------------------------------------------------------------------------
  # 4. Write images in OME-NGFF format
  # ---------------------------------------------------------------------------
  images_written <- FALSE

  if (length(image_names) > 0) {
    for (img_id in image_names) {
      img_obj <- object[[img_id]]
      img_array <- tryCatch(img_obj@image, error = function(e) NULL)

      if (!is.null(img_array) && length(dim(img_array)) == 3) {
        img_zarr_name <- gsub("[^[:alnum:]_]", "_", img_id)
        img_path <- file.path(path, "images", img_zarr_name)
        if (verbose) message("Writing image: images/", img_zarr_name)
        .write_ome_ngff_image(img_array, img_path, name = img_zarr_name,
                              verbose = verbose)
        images_written <- TRUE

        # Also write hires if available
        hires <- attr(img_array, "hires.image")
        if (!is.null(hires) && length(dim(hires)) == 3) {
          hires_name <- paste0(img_zarr_name, "_hires")
          hires_path <- file.path(path, "images", hires_name)
          if (verbose) message("Writing image: images/", hires_name)
          .write_ome_ngff_image(hires, hires_path, name = hires_name,
                                verbose = verbose)
        }
      }
    }
  }

  if (!images_written) {
    .zarr_create_group(file.path(path, "images"))
  }

  # ---------------------------------------------------------------------------
  # 5. Create empty labels and points groups
  # ---------------------------------------------------------------------------
  .zarr_create_group(file.path(path, "labels"))
  .zarr_create_group(file.path(path, "points"))

  # Write points data from misc if present
  point_keys <- grep("^points_", names(object@misc), value = TRUE)
  for (pk in point_keys) {
    pt_name <- sub("^points_", "", pk)
    pt_data <- object@misc[[pk]]
    if (is.data.frame(pt_data) || is.matrix(pt_data)) {
      pt_path <- file.path(path, "points", pt_name)
      if (verbose) message("Writing points: points/", pt_name)
      .write_spatialdata_points(pt_data, pt_path, verbose = verbose)
    }
  }

  if (verbose) {
    message("\nSuccessfully saved SpatialData store")
    message("  Cells: ", ncol(object))
    message("  Features: ", nrow(object))
    message("  Path: ", path)
  }

  invisible(path)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers: reading
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read shapes from a SpatialData shapes directory
#'
#' SpatialData shapes are stored as zarr DataFrames with x, y columns and
#' optionally a radius column. This follows the anndata on-disk DataFrame
#' encoding (vlen-utf8 string index + numeric columns).
#'
#' @param path Path to the shapes zarr group (e.g., \code{shapes/spots})
#' @param verbose Show progress
#'
#' @return A list with \code{coords} (numeric matrix with x, y columns) and
#'   \code{radius} (numeric scalar or NULL). Returns NULL if the shapes
#'   cannot be read.
#'
#' @keywords internal
#'
.read_spatialdata_shapes <- function(path, verbose = TRUE) {
  if (!dir.exists(path)) return(NULL)

  attrs <- .zarr_read_json(file.path(path, ".zattrs"))
  encoding <- attrs[["encoding-type"]] %||% ""

  if (encoding != "dataframe" && encoding != "") {
    if (verbose) message("  Unknown shapes encoding: ", encoding)
    return(NULL)
  }

  # Read the index
  index <- tryCatch(
    .zarr_read_anndata_index_from_attrs(path, attrs),
    error = function(e) NULL
  )

  # Read x, y coordinates
  x_vals <- .sd_read_column(path, "x")
  y_vals <- .sd_read_column(path, "y")

  if (is.null(x_vals) || is.null(y_vals)) {
    if (verbose) message("  Could not find x/y columns in shapes")
    return(NULL)
  }

  n <- length(x_vals)
  coords <- matrix(c(x_vals, y_vals), ncol = 2, dimnames = list(index, c("x", "y")))

  # Read radius if present
  radius <- .sd_read_column(path, "radius")
  if (is.null(radius)) {
    radius <- 1.0
  } else {
    # Use median radius as the representative value
    radius <- stats::median(radius, na.rm = TRUE)
  }

  if (verbose) message("  Read ", n, " shapes (radius: ",
                        round(radius, 2), ")")

  list(coords = coords, radius = radius)
}


#' Read points from a SpatialData points directory
#'
#' Points are stored as zarr DataFrames. They typically represent transcript
#' locations with x, y, and gene columns. Because there can be millions of
#' points, these are stored in the Seurat object's \code{misc} slot rather
#' than as coordinates.
#'
#' @param path Path to the points zarr group
#' @param verbose Show progress
#'
#' @return A data.frame with columns read from the zarr DataFrame, or NULL.
#'
#' @keywords internal
#'
.read_spatialdata_points <- function(path, verbose = TRUE) {
  if (!dir.exists(path)) return(NULL)

  attrs <- .zarr_read_json(file.path(path, ".zattrs"))
  column_order <- attrs[["column-order"]]

  # Read index
  index <- tryCatch(
    .zarr_read_anndata_index_from_attrs(path, attrs),
    error = function(e) NULL
  )

  # Read columns
  if (is.null(column_order)) {
    column_order <- .zarr_list_children(path, "")
    column_order <- setdiff(column_order, c("_index", "index", "__categories"))
  }

  result <- data.frame(row.names = if (!is.null(index)) index else NULL)

  for (col in column_order) {
    col_full <- file.path(path, col)
    if (!dir.exists(col_full) && !file.exists(file.path(col_full, ".zarray"))) {
      next
    }
    tryCatch({
      vals <- .zarr_read_anndata_column(path, col)
      if (!is.null(vals)) {
        if (nrow(result) == 0 && is.null(index)) {
          result <- data.frame(x_placeholder__ = seq_along(vals))
        }
        if (length(vals) == nrow(result) || nrow(result) == 0) {
          result[[col]] <- vals
        }
      }
    }, error = function(e) {
      if (verbose) message("  Could not read points column '", col, "': ",
                           e$message)
    })
  }

  # Clean up placeholder
  result$x_placeholder__ <- NULL

  if (ncol(result) == 0) return(NULL)

  if (verbose) message("  Read ", nrow(result), " points with ",
                        ncol(result), " columns")
  result
}


#' Read an OME-NGFF image from a zarr directory
#'
#' Reads the highest resolution level from an OME-NGFF zarr image. The image
#' is stored as a chunked zarr array with multiscales metadata describing
#' axes and resolution levels. Only the first dataset (highest resolution)
#' is read.
#'
#' @param path Path to the OME-NGFF zarr image group
#' @param verbose Show progress
#'
#' @return A numeric array with dimensions (height, width, channels) normalized
#'   to \[0, 1\], suitable for a Seurat image slot. Returns NULL on failure.
#'
#' @keywords internal
#'
.read_ome_ngff_image <- function(path, verbose = TRUE) {
  if (!dir.exists(path)) return(NULL)

  attrs <- .zarr_read_json(file.path(path, ".zattrs"))

  # Parse multiscales metadata
  multiscales <- attrs$multiscales
  if (is.null(multiscales)) {
    if (verbose) message("  No multiscales metadata; not an OME-NGFF image")
    return(NULL)
  }

  # multiscales is a list; take the first element
  if (is.list(multiscales) && length(multiscales) > 0) {
    ms <- multiscales[[1]]
    if (is.null(ms)) ms <- multiscales
  } else {
    ms <- multiscales
  }

  # Get highest resolution dataset path
  datasets <- ms$datasets
  if (is.null(datasets) || length(datasets) == 0) {
    if (verbose) message("  No datasets in multiscales metadata")
    return(NULL)
  }

  # datasets may be a list of lists or a data.frame
  ds_path <- NULL
  if (is.data.frame(datasets)) {
    ds_path <- as.character(datasets$path[1])
  } else if (is.list(datasets)) {
    if (!is.null(datasets[[1]]$path)) {
      ds_path <- datasets[[1]]$path
    } else if (is.character(datasets[[1]])) {
      ds_path <- datasets[[1]]
    }
  }
  if (is.null(ds_path)) ds_path <- "0"

  # Parse axes metadata
  axes <- ms$axes
  axes_names <- NULL
  if (!is.null(axes)) {
    if (is.data.frame(axes)) {
      axes_names <- as.character(axes$name)
    } else if (is.list(axes)) {
      axes_names <- vapply(axes, function(a) a$name %||% "", character(1))
    }
  }

  # Read the zarr array
  array_path <- file.path(path, ds_path)
  if (!dir.exists(array_path)) {
    if (verbose) message("  Image array not found at: ", ds_path)
    return(NULL)
  }

  tryCatch({
    meta <- .zarr_read_json(file.path(array_path, ".zarray"))
    shape <- if (is.list(meta$shape)) unlist(meta$shape) else meta$shape
    ndim <- length(shape)

    if (ndim < 2 || ndim > 5) {
      if (verbose) message("  Unsupported image dimensions: ", ndim)
      return(NULL)
    }

    # For large images, only read if manageable (< 100M pixels)
    total_pixels <- prod(shape)
    if (total_pixels > 1e8) {
      if (verbose) message("  Image too large (", total_pixels,
                           " elements); skipping")
      return(NULL)
    }

    # Read image data
    raw_array <- .zarr_read_numeric(path, ds_path)

    # Reshape and reorder axes to (H, W, C)
    img <- .ome_ngff_to_hwc(raw_array, shape, axes_names)

    if (is.null(img)) {
      if (verbose) message("  Could not reshape image to (H, W, C)")
      return(NULL)
    }

    # Normalize to [0, 1]
    img_max <- max(img, na.rm = TRUE)
    if (img_max > 1.0) {
      if (img_max <= 255) {
        img <- img / 255.0
      } else if (img_max <= 65535) {
        img <- img / 65535.0
      } else {
        img <- img / img_max
      }
    }
    img[img < 0] <- 0
    img[img > 1] <- 1
    storage.mode(img) <- "double"

    if (verbose) message("  Read image: ", paste(dim(img), collapse = " x "))
    img
  }, error = function(e) {
    if (verbose) message("  Failed to read image: ", e$message)
    NULL
  })
}


#' Convert raw OME-NGFF array to height x width x channels format
#'
#' @param arr Numeric vector or matrix from zarr read
#' @param shape Original zarr shape
#' @param axes_names Character vector of axis names (e.g., c("c", "y", "x"))
#'
#' @return 3D array (H x W x C) or NULL
#'
#' @keywords internal
#'
.ome_ngff_to_hwc <- function(arr, shape, axes_names = NULL) {
  ndim <- length(shape)

  # Convert C-order flat vector to R array with correct element layout.
  # Zarr C-order: last dim varies fastest; R column-major: first dim varies fastest.
  # Fix: reverse dims for fill, then reverse back with aperm.
  img <- array(arr, dim = rev(shape))
  img <- aperm(img, rev(seq_len(ndim)))
  # Now img[i1, i2, ...] correctly maps to the logical zarr element

  if (ndim == 2) {
    # Grayscale: (Y, X) -> (Y, X, 1)
    dim(img) <- c(shape[1], shape[2], 1L)
    return(img)
  }

  if (ndim == 3) {
    # Could be (C, Y, X) or (Y, X, C) or (Y, C, X)
    # Use axes metadata if available
    if (!is.null(axes_names) && length(axes_names) == 3) {
      c_pos <- which(axes_names == "c")
      y_pos <- which(axes_names == "y")
      x_pos <- which(axes_names == "x")

      if (length(c_pos) == 1 && length(y_pos) == 1 && length(x_pos) == 1) {
        img <- aperm(img, c(y_pos, x_pos, c_pos))
        return(img)
      }
    }

    # Heuristic: if first dim is small (1-4), assume (C, Y, X)
    if (shape[1] <= 4) {
      img <- aperm(img, c(2, 3, 1))
      return(img)
    }

    # If last dim is small (1-4), assume (Y, X, C)
    if (shape[3] <= 4) {
      return(img)
    }

    # Default: assume (C, Y, X)
    img <- aperm(img, c(2, 3, 1))
    return(img)
  }

  if (ndim == 4) {
    # (T, C, Y, X) or (C, Z, Y, X) -- take first T/Z slice
    slice <- img[1, , , , drop = FALSE]
    slice_shape <- shape[-1]
    slice <- array(slice, dim = slice_shape)
    # Apply 3D -> HWC permutation
    if (slice_shape[1] <= 4) return(aperm(slice, c(2, 3, 1)))
    if (slice_shape[3] <= 4) return(slice)
    return(aperm(slice, c(2, 3, 1)))
  }

  if (ndim == 5) {
    # (T, C, Z, Y, X) -- take first T and Z slice
    slice <- img[1, , 1, , , drop = FALSE]
    slice_shape <- c(shape[2], shape[4], shape[5])
    slice <- array(slice, dim = slice_shape)
    if (slice_shape[1] <= 4) return(aperm(slice, c(2, 3, 1)))
    if (slice_shape[3] <= 4) return(slice)
    return(aperm(slice, c(2, 3, 1)))
  }

  NULL
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers: writing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write spatial coordinates as a SpatialData shapes directory
#'
#' Writes a zarr DataFrame with x, y, and radius columns following the
#' SpatialData shapes specification.
#'
#' @param coords Numeric matrix with x, y columns (rows = cells/spots)
#' @param radius Numeric scalar for spot radius
#' @param path Output path for the shapes zarr group
#' @param verbose Show progress
#'
#' @keywords internal
#'
.write_spatialdata_shapes <- function(coords, radius = 1.0, path,
                                      verbose = TRUE) {
  n <- nrow(coords)
  compressor <- list(id = "zlib", level = GetCompressionLevel())

  # Index: use rownames if available, else generate
  index <- rownames(coords)
  if (is.null(index)) {
    index <- as.character(seq_len(n) - 1L)
  }

  # Column names for the DataFrame
  col_names <- c("x", "y", "radius")

  .zarr_create_group(path, attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = col_names
  ))

  # Write _index
  .zarr_write_strings(
    dir = file.path(path, "_index"),
    strings = as.character(index),
    compressor = compressor,
    attrs = list(
      `encoding-type` = "string-array",
      `encoding-version` = "0.2.0"
    )
  )

  # Write x column
  .zarr_write_numeric(
    dir = file.path(path, "x"),
    data = as.numeric(coords[, 1]),
    dtype = "<f8",
    compressor = compressor
  )

  # Write y column
  .zarr_write_numeric(
    dir = file.path(path, "y"),
    data = as.numeric(coords[, 2]),
    dtype = "<f8",
    compressor = compressor
  )

  # Write radius column (constant for all spots)
  .zarr_write_numeric(
    dir = file.path(path, "radius"),
    data = rep(as.numeric(radius), n),
    dtype = "<f8",
    compressor = compressor
  )

  if (verbose) message("  Wrote ", n, " shapes with radius ", round(radius, 2))
}


#' Write point coordinates as a SpatialData points directory
#'
#' Writes a zarr DataFrame containing point annotations (e.g., transcript
#' positions) following the SpatialData points specification.
#'
#' @param points Data frame or matrix with point data
#' @param path Output path for the points zarr group
#' @param verbose Show progress
#'
#' @keywords internal
#'
.write_spatialdata_points <- function(points, path, verbose = TRUE) {
  if (!is.data.frame(points)) {
    points <- as.data.frame(points)
  }

  n <- nrow(points)
  compressor <- list(id = "zlib", level = GetCompressionLevel())

  index <- rownames(points) %||% as.character(seq_len(n) - 1L)
  col_names <- colnames(points)

  .zarr_create_group(path, attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = col_names
  ))

  # Write _index
  .zarr_write_strings(
    dir = file.path(path, "_index"),
    strings = as.character(index),
    compressor = compressor,
    attrs = list(
      `encoding-type` = "string-array",
      `encoding-version` = "0.2.0"
    )
  )

  # Write each column
  for (col in col_names) {
    vals <- points[[col]]
    col_dir <- file.path(path, col)

    if (is.numeric(vals)) {
      dtype <- if (is.integer(vals)) "<i4" else "<f8"
      .zarr_write_numeric(dir = col_dir, data = vals, dtype = dtype,
                          compressor = compressor)
    } else if (is.character(vals) || is.factor(vals)) {
      .zarr_write_strings(dir = col_dir, strings = as.character(vals),
                          compressor = compressor,
                          attrs = list(
                            `encoding-type` = "string-array",
                            `encoding-version` = "0.2.0"
                          ))
    }
  }

  if (verbose) message("  Wrote ", n, " points")
}


#' Write an image as an OME-NGFF zarr store
#'
#' Writes a 3D array (H x W x C) as an OME-NGFF zarr image with
#' multiscales metadata. This creates a single resolution level.
#'
#' @param image Numeric 3D array (height, width, channels)
#' @param path Output path for the image zarr group
#' @param name Name for the multiscales metadata
#' @param verbose Show progress
#'
#' @keywords internal
#'
.write_ome_ngff_image <- function(image, path, name = "image",
                                   verbose = TRUE) {
  if (!is.array(image) || length(dim(image)) != 3) {
    if (verbose) message("  Skipping: image must be a 3D array")
    return(invisible(NULL))
  }

  h <- dim(image)[1]
  w <- dim(image)[2]
  c_dim <- dim(image)[3]

  # Create group with multiscales metadata
  multiscales_attrs <- list(
    multiscales = list(
      list(
        axes = list(
          list(name = "c", type = "channel"),
          list(name = "y", type = "space"),
          list(name = "x", type = "space")
        ),
        datasets = list(
          list(
            path = "0",
            coordinateTransformations = list(
              list(type = "identity")
            )
          )
        ),
        name = name,
        version = "0.4"
      )
    )
  )

  .zarr_create_group(path, attrs = multiscales_attrs)

  # Write the image array as (C, Y, X)
  # Convert from R's (H, W, C) to OME-NGFF's (C, Y, X)
  img_cyx <- aperm(image, c(3, 1, 2))

  compressor <- list(id = "zlib", level = GetCompressionLevel())

  array_dir <- file.path(path, "0")
  dir.create(array_dir, recursive = TRUE, showWarnings = FALSE)

  shape <- as.integer(c(c_dim, h, w))
  # Use reasonable chunk sizes
  chunk_y <- min(h, 256L)
  chunk_x <- min(w, 256L)
  chunks <- as.integer(c(c_dim, chunk_y, chunk_x))

  # Write .zarray metadata
  meta <- list(
    zarr_format = 2L,
    shape = as.list(shape),
    chunks = as.list(chunks),
    dtype = "<f8",
    compressor = compressor,
    fill_value = 0.0,
    order = "C",
    filters = NULL,
    dimension_separator = "."
  )
  json <- jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null",
                           pretty = TRUE)
  writeLines(json, file.path(array_dir, ".zarray"))

  # Write chunks
  n_chunks_y <- ceiling(h / chunk_y)
  n_chunks_x <- ceiling(w / chunk_x)

  for (cy in seq_len(n_chunks_y) - 1L) {
    for (cx in seq_len(n_chunks_x) - 1L) {
      y_start <- cy * chunk_y + 1L
      y_end <- min((cy + 1L) * chunk_y, h)
      x_start <- cx * chunk_x + 1L
      x_end <- min((cx + 1L) * chunk_x, w)

      chunk_data <- img_cyx[, y_start:y_end, x_start:x_end, drop = FALSE]

      # Pad to full chunk size if needed
      actual_cy <- y_end - y_start + 1L
      actual_cx <- x_end - x_start + 1L
      if (actual_cy < chunk_y || actual_cx < chunk_x) {
        padded <- array(0.0, dim = c(c_dim, chunk_y, chunk_x))
        padded[, seq_len(actual_cy), seq_len(actual_cx)] <- chunk_data
        chunk_data <- padded
      }

      # Flatten to C order (row-major): for (C, Y, X), X varies fastest
      # R column-major as.vector makes dim[1] fastest; aperm(c(3,2,1)) puts X first
      flat <- as.double(as.vector(aperm(chunk_data, c(3, 2, 1))))

      raw_data <- writeBin(flat, raw(), size = 8L, endian = "little")
      raw_data <- .zarr_compress(raw_data, compressor)

      chunk_file <- file.path(array_dir, paste0("0.", cy, ".", cx))
      writeBin(raw_data, chunk_file)
    }
  }

  if (verbose) message("  Wrote OME-NGFF image: ",
                        c_dim, " x ", h, " x ", w, " (C x Y x X)")
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers: coordinate / image attachment
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Add spatial coordinates from SpatialData to Seurat object
#'
#' @param obj Seurat object
#' @param coords Coordinate matrix (n x 2) with x, y columns
#' @param radius Spot radius
#' @param region_name Name of the spatial region
#' @param cell_names Cell names in the Seurat object
#' @param instance_key Column in obs linking to shapes index
#' @param verbose Show progress
#'
#' @return Modified Seurat object
#'
#' @keywords internal
#'
.add_spatialdata_coords <- function(obj, coords, radius = 1.0,
                                    region_name = "spots",
                                    cell_names = NULL, instance_key = NULL,
                                    verbose = TRUE) {
  if (is.null(cell_names)) cell_names <- Cells(obj)
  n_cells <- length(cell_names)

  # Try to match coords to cells via instance_key or row order
  if (!is.null(instance_key) && instance_key %in% colnames(obj@meta.data)) {
    # Map via instance key
    instance_ids <- as.character(obj@meta.data[[instance_key]])
    coord_idx <- match(instance_ids, rownames(coords))
    if (any(!is.na(coord_idx))) {
      matched_coords <- coords[coord_idx, , drop = FALSE]
      rownames(matched_coords) <- cell_names
    } else {
      matched_coords <- NULL
    }
  } else if (nrow(coords) == n_cells) {
    # Direct row-order match
    matched_coords <- coords
    rownames(matched_coords) <- cell_names
  } else if (!is.null(rownames(coords))) {
    # Try matching by rownames
    common <- intersect(cell_names, rownames(coords))
    if (length(common) > 0 && length(common) == n_cells) {
      matched_coords <- coords[cell_names, , drop = FALSE]
    } else if (length(common) > n_cells * 0.5) {
      # Partial match - use what we have
      matched_coords <- matrix(NA_real_, nrow = n_cells, ncol = 2,
                               dimnames = list(cell_names, c("x", "y")))
      matched_coords[common, ] <- coords[common, ]
      if (verbose) message("  Matched ", length(common), "/", n_cells,
                           " cells to shape coordinates")
    } else {
      matched_coords <- NULL
    }
  } else {
    matched_coords <- NULL
  }

  if (is.null(matched_coords)) {
    if (verbose) message("  Could not match shapes to cells; ",
                         "storing in meta.data")
    # Store coordinates in meta.data as fallback
    if (nrow(coords) == n_cells) {
      obj@meta.data$spatial_x <- coords[, 1]
      obj@meta.data$spatial_y <- coords[, 2]
    }
    return(obj)
  }

  # Create FOV object
  coord_df <- data.frame(
    x = matched_coords[, 1],
    y = matched_coords[, 2],
    row.names = cell_names,
    stringsAsFactors = FALSE
  )

  assay_name <- DefaultAssay(obj)
  image_key <- SanitizeImageKey(region_name)

  tryCatch({
    fov <- CreateFOV(
      coords = coord_df[, c("x", "y"), drop = FALSE],
      type = "centroids",
      radius = as.numeric(radius),
      assay = assay_name,
      key = image_key
    )
    obj[[region_name]] <- fov
    if (verbose) message("  Added FOV '", region_name, "' with ",
                         nrow(coord_df), " centroids")
  }, error = function(e) {
    if (verbose) message("  Could not create FOV (", e$message,
                         "); storing coordinates in meta.data")
    obj@meta.data$spatial_x <<- matched_coords[, 1]
    obj@meta.data$spatial_y <<- matched_coords[, 2]
  })

  obj
}


#' Attach an image array to a Seurat object
#'
#' @param obj Seurat object (must already have spatial coordinates)
#' @param img_data 3D array (H, W, C) normalized to [0, 1]
#' @param img_name Name for the image
#' @param verbose Show progress
#'
#' @return Modified Seurat object
#'
#' @keywords internal
#'
.attach_spatialdata_image <- function(obj, img_data, img_name = "image",
                                     verbose = TRUE) {
  # Store in misc since constructing full Visium image objects requires

  # scale factors that may not be available from SpatialData
  obj@misc[[paste0("image_", img_name)]] <- img_data
  if (verbose) message("  Stored image '", img_name, "' in misc (",
                        paste(dim(img_data), collapse = " x "), ")")
  obj
}


#' Read an anndata index given pre-read attrs
#'
#' @param store_path Path to zarr group
#' @param attrs Pre-parsed .zattrs list
#'
#' @return Character vector of index values
#'
#' @keywords internal
#'
.zarr_read_anndata_index_from_attrs <- function(store_path, attrs) {
  index_col <- attrs[["_index"]] %||% "_index"

  index_path <- file.path(store_path, index_col)
  if (dir.exists(index_path) || file.exists(file.path(index_path, ".zarray"))) {
    return(as.character(.zarr_read_strings(store_path, index_col)))
  }

  # Fallback
  index_path <- file.path(store_path, "_index")
  if (dir.exists(index_path)) {
    return(as.character(.zarr_read_strings(store_path, "_index")))
  }

  index_path <- file.path(store_path, "index")
  if (dir.exists(index_path)) {
    return(as.character(.zarr_read_strings(store_path, "index")))
  }

  NULL
}


#' Read a single column from a zarr DataFrame
#'
#' Dispatches to numeric or string reading based on dtype.
#'
#' @param store_path Path to the zarr DataFrame group
#' @param col_name Column name
#'
#' @return Vector of values, or NULL if column not found
#'
#' @keywords internal
#'
.sd_read_column <- function(store_path, col_name) {
  col_path <- file.path(store_path, col_name)

  if (.zarr_node_type(store_path, col_name) == "missing") {
    return(NULL)
  }

  # Use the general anndata column reader
  .zarr_read_anndata_column(store_path, col_name)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Direct conversion functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert a SpatialData zarr store to h5ad format
#'
#' Extracts the anndata table from a SpatialData store and writes it as an
#' h5ad file. Spatial coordinates and images are embedded in the h5ad using
#' standard scanpy/squidpy conventions (\code{obsm/spatial}, \code{uns/spatial}).
#'
#' @param source Path to input SpatialData .zarr directory
#' @param dest Path for output .h5ad file
#' @param table Name of the table within the SpatialData store (default: "table")
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
SpatialDataToH5AD <- function(source, dest, table = "table",
                               overwrite = FALSE, verbose = TRUE) {
  if (!dir.exists(source)) {
    stop("SpatialData store not found: ", source, call. = FALSE)
  }
  if (file.exists(dest)) {
    if (!overwrite) {
      stop("Output file exists: ", dest,
           ". Use overwrite = TRUE.", call. = FALSE)
    }
    file.remove(dest)
  }

  if (verbose) message("Converting SpatialData -> h5ad")
  obj <- readSpatialData(source, table = table, verbose = verbose)
  writeH5AD(obj, dest, overwrite = overwrite, verbose = verbose)

  invisible(dest)
}


#' Convert an h5ad file to SpatialData zarr format
#'
#' Reads an h5ad file and writes it as a SpatialData zarr store. If the h5ad
#' contains spatial information (\code{obsm/spatial}, Visium-style
#' \code{uns/spatial}), it is placed in the appropriate SpatialData elements.
#'
#' @param source Path to input .h5ad file
#' @param dest Path for output SpatialData .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5ADToSpatialData <- function(source, dest, overwrite = FALSE,
                               verbose = TRUE) {
  if (!file.exists(source)) {
    stop("Input file not found: ", source, call. = FALSE)
  }
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite) {
      stop("Output exists: ", dest,
           ". Use overwrite = TRUE.", call. = FALSE)
    }
    unlink(dest, recursive = TRUE)
  }

  if (verbose) message("Converting h5ad -> SpatialData")
  obj <- readH5AD(source, verbose = verbose)
  writeSpatialData(obj, dest, overwrite = overwrite, verbose = verbose)

  invisible(dest)
}


#' Convert a SpatialData zarr store to h5Seurat format
#'
#' Reads a SpatialData store and writes it as an h5Seurat file.
#'
#' @param source Path to input SpatialData .zarr directory
#' @param dest Path for output .h5Seurat file
#' @param table Name of the table within the SpatialData store (default: "table")
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
SpatialDataToH5Seurat <- function(source, dest, table = "table",
                                   overwrite = FALSE, verbose = TRUE) {
  if (!dir.exists(source)) {
    stop("SpatialData store not found: ", source, call. = FALSE)
  }
  if (file.exists(dest)) {
    if (!overwrite) {
      stop("Output file exists: ", dest,
           ". Use overwrite = TRUE.", call. = FALSE)
    }
    file.remove(dest)
  }

  if (verbose) message("Converting SpatialData -> h5Seurat")
  obj <- readSpatialData(source, table = table, verbose = verbose)
  writeH5Seurat(obj, dest, overwrite = overwrite, verbose = verbose)

  invisible(dest)
}


#' Convert an h5Seurat file to SpatialData zarr format
#'
#' Reads an h5Seurat file and writes it as a SpatialData zarr store.
#'
#' @param source Path to input .h5Seurat file
#' @param dest Path for output SpatialData .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5SeuratToSpatialData <- function(source, dest, overwrite = FALSE,
                                   verbose = TRUE) {
  if (!file.exists(source)) {
    stop("Input file not found: ", source, call. = FALSE)
  }
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite) {
      stop("Output exists: ", dest,
           ". Use overwrite = TRUE.", call. = FALSE)
    }
    unlink(dest, recursive = TRUE)
  }

  if (verbose) message("Converting h5Seurat -> SpatialData")
  obj <- readH5Seurat(source, verbose = verbose)
  writeSpatialData(obj, dest, overwrite = overwrite, verbose = verbose)

  invisible(dest)
}


#' Convert a SpatialData zarr store to a standard anndata zarr store
#'
#' Extracts the anndata table from a SpatialData store and writes it as
#' a standalone anndata zarr store (without the SpatialData wrapper).
#'
#' @param source Path to input SpatialData .zarr directory
#' @param dest Path for output .zarr directory
#' @param table Name of the table within the SpatialData store (default: "table")
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
SpatialDataToZarr <- function(source, dest, table = "table",
                               overwrite = FALSE, verbose = TRUE) {
  if (!dir.exists(source)) {
    stop("SpatialData store not found: ", source, call. = FALSE)
  }
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite) {
      stop("Output exists: ", dest,
           ". Use overwrite = TRUE.", call. = FALSE)
    }
    unlink(dest, recursive = TRUE)
  }

  if (verbose) message("Converting SpatialData -> Zarr (extracting table)")
  obj <- readSpatialData(source, table = table, verbose = verbose)
  writeZarr(obj, dest, overwrite = overwrite, verbose = verbose)

  invisible(dest)
}


#' Wrap a standard anndata zarr store as a SpatialData store
#'
#' Takes an existing anndata zarr store and wraps it as a SpatialData store
#' by creating the SpatialData directory layout and placing the anndata
#' store as a table.
#'
#' @param source Path to input anndata .zarr directory
#' @param dest Path for output SpatialData .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToSpatialData <- function(source, dest, overwrite = FALSE,
                               verbose = TRUE) {
  if (!dir.exists(source)) {
    stop("Zarr store not found: ", source, call. = FALSE)
  }
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite) {
      stop("Output exists: ", dest,
           ". Use overwrite = TRUE.", call. = FALSE)
    }
    unlink(dest, recursive = TRUE)
  }

  if (verbose) message("Converting Zarr -> SpatialData (wrapping as table)")
  obj <- readZarr(source, verbose = verbose)
  writeSpatialData(obj, dest, overwrite = overwrite, verbose = verbose)

  invisible(dest)
}
