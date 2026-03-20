#' @include zzz.R
#' @include Convert.R
#' @importFrom hdf5r H5File h5attr h5types
#' @importFrom Seurat CreateSeuratObject CreateAssayObject Images GetTissueCoordinates scalefactors
#' @importFrom SeuratObject AddMetaData Cells CreateFOV
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spatial Data Conversion Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert spatial coordinates from h5ad to Seurat format
#'
#' @param h5ad_file H5AD file handle or path
#' @param seurat_obj Seurat object to add spatial data to
#' @param assay_name Name of the assay to associate spatial data with
#' @param verbose Print progress messages
#'
#' @return Modified Seurat object with spatial data
#' @export
H5ADSpatialToSeurat <- function(h5ad_file, seurat_obj = NULL,
                                       assay_name = "Spatial", verbose = TRUE) {

  # Open h5ad file if path provided
  if (is.character(h5ad_file)) {
    h5ad <- H5File$new(h5ad_file, mode = "r")
    on.exit(h5ad$close_all())
  } else {
    h5ad <- h5ad_file
  }

  # Check for spatial coordinates in obsm
  has_spatial_coords <- FALSE
  spatial_coords <- NULL

  cell_names <- NULL

  if (!is.null(seurat_obj)) {
    cell_names <- Cells(seurat_obj)
  } else if (h5ad$exists("obs")) {
    obs_group <- h5ad[["obs"]]
    if (obs_group$exists("_index")) {
      cell_names <- as.character(obs_group[["_index"]][])
    } else if (obs_group$exists("index")) {
      cell_names <- as.character(obs_group[["index"]][])
    }
  }

  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    spatial_obj <- h5ad[["obsm"]][["spatial"]]

    # Handle different storage formats for obsm/spatial
    spatial_coords <- tryCatch({
      if (inherits(spatial_obj, "H5Group")) {
        enc <- tryCatch(h5attr(spatial_obj, "encoding-type"), error = function(e) "")
        if (enc == "dataframe") {
          # DataFrame format: read X/Y (or first two numeric) columns
          sp_names <- names(spatial_obj)
          xy_cols <- intersect(c("X", "Y", "x", "y"), sp_names)
          if (length(xy_cols) >= 2L) {
            cbind(spatial_obj[[xy_cols[1]]][], spatial_obj[[xy_cols[2]]][])
          } else {
            # Use first two numeric columns
            numeric_cols <- character(0)
            for (cn in sp_names) {
              if (inherits(spatial_obj[[cn]], "H5D")) {
                numeric_cols <- c(numeric_cols, cn)
              }
              if (length(numeric_cols) >= 2L) break
            }
            if (length(numeric_cols) >= 2L) {
              cbind(spatial_obj[[numeric_cols[1]]][], spatial_obj[[numeric_cols[2]]][])
            } else {
              NULL
            }
          }
        } else {
          # Sparse array — try ReadH5ADMatrix-style read
          NULL
        }
      } else if (inherits(spatial_obj, "H5D")) {
        # Dense array
        spatial_obj[,]
      } else {
        NULL
      }
    }, error = function(e) {
      if (verbose) warning("Could not read obsm/spatial: ", e$message, immediate. = TRUE)
      NULL
    })

    if (!is.null(spatial_coords)) {
      spatial_coords <- as.matrix(spatial_coords)
      if (!is.null(dim(spatial_coords))) {
        if (!is.null(cell_names)) {
          if (nrow(spatial_coords) == length(cell_names) && ncol(spatial_coords) == 2L) {
            rownames(spatial_coords) <- cell_names
          } else if (ncol(spatial_coords) == length(cell_names) && nrow(spatial_coords) == 2L) {
            spatial_coords <- t(spatial_coords)
            rownames(spatial_coords) <- cell_names
          }
        } else if (nrow(spatial_coords) == 2L) {
          spatial_coords <- t(spatial_coords)
        }
      }

      has_spatial_coords <- TRUE
      if (verbose) message("Found spatial coordinates in obsm['spatial']")
    }
  }

  # Check for Visium-style spatial data in uns
  has_visium_data <- FALSE
  visium_data <- list()

  if (h5ad$exists("uns") && "spatial" %in% names(h5ad[["uns"]])) {
    spatial_uns <- h5ad[["uns/spatial"]]
    library_ids <- names(spatial_uns)
    has_visium_data <- length(library_ids) > 0

    if (has_visium_data && verbose) {
      message("Found Visium spatial data for libraries: ", paste(library_ids, collapse = ", "))
    }

    # Process each library
    for (lib_id in library_ids) {
      lib_data <- list()
      lib_group <- spatial_uns[[lib_id]]

      # Get scale factors
      if ("scalefactors" %in% names(lib_group)) {
        sf_group <- lib_group[["scalefactors"]]
        lib_data$scalefactors <- list()

        for (sf_name in names(sf_group)) {
          lib_data$scalefactors[[sf_name]] <- sf_group[[sf_name]][]
        }
      }

      # Get image metadata (not loading actual images here)
      if ("images" %in% names(lib_group)) {
        lib_data$has_images <- TRUE
        image_names <- names(lib_group[["images"]])
        lib_data$image_names <- image_names
      }

      # Get metadata
      if ("metadata" %in% names(lib_group)) {
        meta_group <- lib_group[["metadata"]]
        lib_data$metadata <- list()

        for (meta_name in names(meta_group)) {
          lib_data$metadata[[meta_name]] <- meta_group[[meta_name]][]
        }
      }

      visium_data[[lib_id]] <- lib_data
    }
  }

  # Create or modify Seurat object
  if (is.null(seurat_obj)) {
    if (!has_spatial_coords) {
      stop("No spatial coordinates found and no Seurat object provided")
    }

    # Create minimal Seurat object with spatial coordinates
    n_cells <- nrow(spatial_coords)
    n_genes <- 100  # Placeholder

    counts <- matrix(0, nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))

    seurat_obj <- CreateSeuratObject(counts = counts, assay = assay_name)
    cell_names <- Cells(seurat_obj)
  }

  # Add spatial coordinates to Seurat object
  if (has_spatial_coords) {
    if (!is.null(cell_names) && nrow(spatial_coords) == length(cell_names)) {
      spatial_coords <- spatial_coords[cell_names, , drop = FALSE]
    }
    if (ncol(spatial_coords) == 2L) {
      colnames(spatial_coords) <- c('imagerow', 'imagecol')
    }

    # Determine technology type based on data structure
    technology <- DetectSpatialTechnology(spatial_coords, visium_data)

    if (technology == "Visium") {
      # Create Visium-specific spatial object
      seurat_obj <- AddVisiumSpatialData(
        seurat_obj,
        spatial_coords,
        visium_data,
        h5ad,
        assay_name = assay_name,
        verbose = verbose
      )
    } else if (technology == "SlideSeq") {
      # Create SlideSeq-specific spatial object
      seurat_obj <- AddSlideSeqSpatialData(
        seurat_obj,
        spatial_coords,
        assay_name = assay_name,
        verbose = verbose
      )
    } else {
      # Generic spatial data
      seurat_obj <- AddGenericSpatialData(
        seurat_obj,
        spatial_coords,
        assay_name = assay_name,
        verbose = verbose
      )
    }
  }

  return(seurat_obj)
}

#' Convert Seurat spatial data to h5ad format
#'
#' @param seurat_obj Seurat object with spatial data
#' @param h5ad_file H5AD file handle or path to write to
#' @param library_id Library ID for spatial data (default: "library_1")
#' @param verbose Print progress messages
#'
#' @export
SeuratSpatialToH5AD <- function(seurat_obj, h5ad_file,
                                       library_id = "library_1",
                                       verbose = TRUE) {

  # Open h5ad file if path provided
  if (is.character(h5ad_file)) {
    h5ad <- H5File$new(h5ad_file, mode = "r+")
    on.exit(h5ad$close_all())
  } else {
    h5ad <- h5ad_file
  }

  # Check for spatial data in Seurat object
  images <- Images(seurat_obj)

  if (length(images) == 0) {
    if (verbose) message("No spatial data found in Seurat object")
    return(invisible(NULL))
  }

  if (verbose) message("Found ", length(images), " spatial image(s)")

  # Use actual library ID from the image name
  if (library_id == "library_1" && length(images) > 0) {
    library_id <- images[1]
  }

  # Process first image (extend for multiple images later)
  img_obj <- seurat_obj[[images[1]]]

  # Extract coordinates
  coords <- GetTissueCoordinates(img_obj)

  # Convert to h5ad format (cells x 2 matrix)
  # AnnData/squidpy convention: [X, Y]
  # VisiumV2 returns x/y/cell; VisiumV1 returns imagecol/imagerow
  if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    spatial_matrix <- as.matrix(coords[, c("imagecol", "imagerow")])
  } else if (all(c("x", "y") %in% colnames(coords))) {
    spatial_matrix <- as.matrix(coords[, c("x", "y")])
  } else {
    # Fallback: use first two numeric columns
    num_cols <- sapply(coords, is.numeric)
    spatial_matrix <- as.matrix(coords[, which(num_cols)[1:2]])
  }

  # Create obsm group if not exists
  if (!h5ad$exists("obsm")) {
    h5ad$create_group("obsm")
  }

  # Write spatial coordinates
  if (h5ad[["obsm"]]$exists("spatial")) {
    h5ad[["obsm"]]$link_delete("spatial")
  }

  # scTranspose: hdf5r writes R (n_cells, 2) as HDF5 (2, n_cells);
  # anndata/squidpy expect (n_cells, 2), so write t(spatial_matrix)
  h5ad[["obsm"]]$create_dataset(
    name = "spatial",
    robj = t(spatial_matrix),
    dtype = h5types$H5T_NATIVE_DOUBLE
  )
  # anndata requires encoding attributes on obsm datasets
  h5ad[["obsm"]][["spatial"]]$create_attr(
    attr_name = 'encoding-type', robj = 'array',
    dtype = CachedGuessDType(x = 'array'), space = ScalarSpace()
  )
  h5ad[["obsm"]][["spatial"]]$create_attr(
    attr_name = 'encoding-version', robj = '0.2.0',
    dtype = CachedGuessDType(x = '0.2.0'), space = ScalarSpace()
  )

  if (verbose) message("Wrote spatial coordinates to obsm['spatial']")

  # Create uns/spatial structure for Visium-like data
  if (!h5ad$exists("uns")) {
    h5ad$create_group("uns")
  }

  if (!h5ad[["uns"]]$exists("spatial")) {
    h5ad[["uns"]]$create_group("spatial")
  }

  # Create library-specific group
  spatial_group <- h5ad[["uns/spatial"]]

  if (spatial_group$exists(library_id)) {
    spatial_group$link_delete(library_id)
  }

  lib_group <- spatial_group$create_group(library_id)

  # Add scale factors if available
  if (inherits(img_obj, "VisiumV1") || inherits(img_obj, "VisiumV2") || inherits(img_obj, "FOV")) {
    scalefactors <- GetScaleFactors(img_obj)

    if (!is.null(scalefactors) && length(scalefactors) > 0) {
      sf_group <- lib_group$create_group("scalefactors")

      # Map Seurat scale factor names to h5ad/scanpy convention
      sf_name_map <- c(
        spot     = "spot_diameter_fullres",
        fiducial = "fiducial_diameter_fullres",
        hires    = "tissue_hires_scalef",
        lowres   = "tissue_lowres_scalef"
      )

      for (sf_name in names(scalefactors)) {
        val <- as.numeric(scalefactors[[sf_name]])
        if (is.finite(val)) {
          # Use mapped name if available, otherwise keep original
          dst_name <- if (!is.null(sf_name_map[[sf_name]])) sf_name_map[[sf_name]] else sf_name
          # Write as HDF5 scalar (shape=()) not 1-element array (shape=(1,)).
          # squidpy/scanpy expect scalars; arrays cause plotting issues.
          sf_group$create_dataset(
            name = dst_name,
            robj = val,
            dtype = h5types$H5T_NATIVE_DOUBLE,
            space = ScalarSpace(),
            chunk_dims = NULL
          )
        }
      }

      if (verbose) message("Wrote scale factors")
    }

    # Write tissue images
    # R image arrays are (height, width, channels) in column-major order.
    # hdf5r reverses dimensions for HDF5 row-major, so Python would read
    # (channels, width, height). aperm reverses the R dims so that after
    # HDF5 transposition Python gets (height, width, channels).
    tryCatch({
      img_data <- img_obj@image
      if (!is.null(img_data) && length(dim(img_data)) == 3) {
        images_group <- lib_group$create_group("images")
        images_group$create_dataset(
          name = "lowres",
          robj = aperm(img_data, c(3, 2, 1)),
          dtype = h5types$H5T_NATIVE_DOUBLE
        )
        if (verbose) message("Wrote lowres tissue image")

        # Write hires if available
        hires <- attr(img_data, "hires.image")
        if (!is.null(hires)) {
          images_group$create_dataset(
            name = "hires",
            robj = aperm(hires, c(3, 2, 1)),
            dtype = h5types$H5T_NATIVE_DOUBLE
          )
          if (verbose) message("Wrote hires tissue image")
        }
      }
    }, error = function(e) {
      if (verbose) message("Could not write tissue images: ", e$message)
    })
  }

  # Add metadata
  meta_group <- lib_group$create_group("metadata")
  meta_group$create_dataset(
    name = "technology",
    robj = class(img_obj)[1]
  )

  if (verbose) message("Spatial data conversion complete")

  return(invisible(NULL))
}

#' Detect spatial technology from data structure
#'
#' @param coords Spatial coordinates matrix
#' @param metadata Additional metadata
#' @return Technology type string
#' @keywords internal
DetectSpatialTechnology <- function(coords, metadata = NULL) {
  # Check metadata for Visium indicators
  if (!is.null(metadata) && length(metadata) > 0) {
    for (lib in metadata) {
      sf_names <- names(lib$scalefactors)
      # Check for both full names and short names
      visium_sf <- c("tissue_hires_scalef", "hires", "tissue_lowres_scalef", "lowres")
      if (any(visium_sf %in% sf_names)) return("Visium")
      # Check for images (Visium has tissue images, Slide-seq does not)
      if (isTRUE(lib$has_images)) return("Visium")
    }
  }

  # Check coordinate patterns (gridded vs continuous)
  if (ncol(coords) != 2) return("Generic")

  n_unique <- sapply(1:2, function(i) length(unique(coords[, i])))
  n_cells <- nrow(coords)
  # Visium: gridded spots with limited unique coordinate values
  if (all(n_unique < n_cells / 2)) {
    return("Visium")
  }
  # Classify by cell count and coordinate density
  if (n_cells > 50000) {
    return("HighDensitySpatial")  # Stereo-seq, MERFISH, Slide-seq V2
  }
  "Generic"  # IMC, CODEX, small Slide-seq, etc.
}

#' Read HDF5 image dataset with proper dimension handling
#'
#' @param dataset HDF5 dataset containing image data
#' @return Array with proper dimensions (height x width x channels)
#' @keywords internal
ReadH5ImageDataset <- function(dataset) {
  arr <- dataset$read()
  dims <- dataset$dims

  # Reorder dimensions if needed
  # HDF5 image stored as (h, w, 3); hdf5r reverses to R dims (3, w, h).
  # Need aperm to get Seurat's expected (h, w, 3).
  d <- dim(arr)
  if (length(d) == 3L && d[1] <= 4L) {
    arr <- aperm(arr, c(3L, 2L, 1L))
  }

  # Normalize to [0, 1]
  arr[arr < 0] <- 0
  arr[arr > 1] <- 1
  storage.mode(arr) <- 'double'
  arr
}

#' Sanitize string to valid Seurat image key
#'
#' @param x String to sanitize
#' @return Valid Seurat key ending with underscore
#' @keywords internal
SanitizeImageKey <- function(x) {
  key <- gsub(pattern = '\\.+', replacement = '_', x = make.names(x))
  key <- gsub(pattern = '_+', replacement = '_', x = key)
  key <- gsub(pattern = '^_|_$', replacement = '', x = key)
  if (nchar(key) == 0L) key <- 'spatial'
  paste0(key, '_')
}

#' Add Visium spatial data to Seurat object
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param visium_data Visium-specific metadata
#' @param assay_name Assay name
#' @param verbose Print messages
#'
#' @return Modified Seurat object
#' @keywords internal
AddVisiumSpatialData <- function(seurat_obj, coords, visium_data, h5ad,
                                 assay_name = "Spatial", verbose = TRUE) {

  coords <- as.matrix(coords)
  cell_names <- Cells(seurat_obj)

  if (nrow(coords) != length(cell_names) || ncol(coords) != 2L) {
    warning("Coordinate count doesn't match cell count")
    return(seurat_obj)
  }

  coords <- coords[cell_names, , drop = FALSE]
  colnames(coords) <- c('imagerow', 'imagecol')

  seurat_obj@meta.data$spatial_x <- coords[, 'imagecol']
  seurat_obj@meta.data$spatial_y <- coords[, 'imagerow']

  seurat_obj@misc$spatial_technology <- "Visium"
  seurat_obj@misc$spatial_metadata <- visium_data

  if (is.null(visium_data) || length(visium_data) == 0L || !h5ad$exists("uns") ||
      !h5ad[["uns"]]$exists("spatial")) {
    if (verbose) {
      message("Stored Visium coordinates; no image data available in h5ad")
    }
    return(seurat_obj)
  }

  spatial_uns <- tryCatch(h5ad[["uns/spatial"]], error = function(e) NULL)
  if (is.null(spatial_uns)) {
    return(seurat_obj)
  }

  for (lib_id in names(visium_data)) {
    lib_info <- visium_data[[lib_id]]
    lib_group <- tryCatch(spatial_uns[[lib_id]], error = function(e) NULL)
    if (is.null(lib_group) || !lib_group$exists("images")) {
      next
    }

    images_group <- lib_group[["images"]]
    lowres_img <- if (images_group$exists("lowres")) {
      ReadH5ImageDataset(images_group[["lowres"]])
    } else {
      NULL
    }
    if (is.null(lowres_img)) {
      next
    }

    hires_img <- if (images_group$exists("hires")) {
      ReadH5ImageDataset(images_group[["hires"]])
    } else {
      NULL
    }

    sf_values <- lib_info$scalefactors %||% list()
    spot <- sf_values[['spot_diameter_fullres']] %||% sf_values[['spot']]
    fiducial <- sf_values[['fiducial_diameter_fullres']] %||% sf_values[['fiducial']]
    hires_sf <- sf_values[['tissue_hires_scalef']] %||% sf_values[['hires']]
    lowres_sf <- sf_values[['tissue_lowres_scalef']] %||% sf_values[['lowres']]

    scales <- scalefactors(
      spot = as.numeric(spot %||% NA_real_),
      fiducial = as.numeric(fiducial %||% NA_real_),
      hires = as.numeric(hires_sf %||% NA_real_),
      lowres = as.numeric(lowres_sf %||% NA_real_)
    )

    radius <- as.numeric(scales[['spot']]) / 2.0
    if (!is.finite(radius)) {
      radius <- 1
    }

    coord_df <- data.frame(
      imagerow = coords[, 'imagerow'],
      imagecol = coords[, 'imagecol'],
      row.names = cell_names,
      stringsAsFactors = FALSE
    )

    image_key <- SanitizeImageKey(lib_id)
    fov <- CreateFOV(
      coords = coord_df[, c('imagerow', 'imagecol'), drop = FALSE],
      type = 'centroids',
      radius = radius,
      assay = assay_name,
      key = image_key
    )

    # Try VisiumV2 first (primary path), fall back to VisiumV1 (legacy)
    visium_image <- tryCatch({
      # VisiumV2 extends FOV: uses boundaries/centroids, no coordinates slot
      img <- new(
        Class = 'VisiumV2',
        image = lowres_img,
        scale.factors = scales,
        molecules = fov@molecules,
        boundaries = fov@boundaries,
        coords_x_orientation = 'horizontal',
        assay = assay_name,
        key = image_key
      )
      img <- img[cell_names]
      if (!is.null(hires_img)) {
        attr(img@image, 'hires.image') <- hires_img
      }
      img
    }, error = function(e) {
      # Legacy fallback: VisiumV1 has coordinates + spot.radius slots
      tryCatch({
        img <- new(
          Class = 'VisiumV1',
          image = lowres_img,
          scale.factors = scales,
          coordinates = coord_df,
          spot.radius = radius,
          assay = assay_name,
          key = image_key
        )
        if (!is.null(hires_img)) {
          attr(img@image, 'hires.image') <- hires_img
        }
        if (verbose) message("  Using VisiumV1 (legacy) for library '", lib_id, "'")
        img
      }, error = function(e2) {
        # Final fallback: plain FOV
        if (verbose) message("  Using FOV fallback for library '", lib_id, "'")
        fov
      })
    })
    seurat_obj[[lib_id]] <- visium_image

    if (verbose) {
      message("  Added Visium image for library '", lib_id, "'")
    }
  }

  if (verbose) {
    message("Added Visium spatial data:")
    message("  - Spots: ", nrow(coords))
    message("  - Libraries: ", paste(names(visium_data), collapse = ", "))
  }

  return(seurat_obj)
}

#' Add SlideSeq spatial data to Seurat object
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param assay_name Assay name
#' @param verbose Print messages
#'
#' @return Modified Seurat object
#' @keywords internal
AddSlideSeqSpatialData <- function(seurat_obj, coords,
                                   assay_name = "Spatial", verbose = TRUE) {

  # Add continuous coordinates to metadata
  if (nrow(coords) == ncol(seurat_obj)) {
    seurat_obj@meta.data$spatial_x <- coords[, 1]
    seurat_obj@meta.data$spatial_y <- coords[, 2]

    seurat_obj@misc$spatial_technology <- "SlideSeq"

    if (verbose) {
      message("Added SlideSeq spatial data:")
      message("  - Beads: ", nrow(coords))
      message("  - X range: ", round(min(coords[,1]), 2), " - ", round(max(coords[,1]), 2))
      message("  - Y range: ", round(min(coords[,2]), 2), " - ", round(max(coords[,2]), 2))
    }
  } else {
    warning("Coordinate count doesn't match cell count")
  }

  return(seurat_obj)
}

#' Add generic spatial data to Seurat object
#'
#' @param seurat_obj Seurat object
#' @param coords Spatial coordinates
#' @param assay_name Assay name
#' @param verbose Print messages
#'
#' @return Modified Seurat object
#' @keywords internal
AddGenericSpatialData <- function(seurat_obj, coords,
                                  assay_name = "Spatial", verbose = TRUE) {

  # Add coordinates to metadata
  if (nrow(coords) == ncol(seurat_obj)) {
    seurat_obj@meta.data$spatial_x <- coords[, 1]
    seurat_obj@meta.data$spatial_y <- coords[, 2]

    seurat_obj@misc$spatial_technology <- "Generic"

    if (verbose) {
      message("Added generic spatial data:")
      message("  - Points: ", nrow(coords))
    }
  } else {
    warning("Coordinate count doesn't match cell count")
  }

  return(seurat_obj)
}

#' Get scale factors from Visium object
#'
#' @param visium_obj Visium image object
#'
#' @return List of scale factors, or NULL if unavailable
#' @keywords internal
#'
#' @importFrom Seurat scalefactors
#'
GetScaleFactors <- function(visium_obj) {
  tryCatch({
    # VisiumV2 stores scale.factors as a named list with class "scalefactors"
    sf <- visium_obj@scale.factors
    if (is.list(sf) && length(sf) > 0) {
      return(as.list(sf))
    }
    # Fallback: try scalefactors() generic
    sf <- scalefactors(visium_obj)
    if (!is.null(sf)) {
      return(as.list(sf))
    }
  }, error = function(e) NULL)

  NULL
}
