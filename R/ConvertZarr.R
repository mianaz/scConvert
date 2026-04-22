#' @include Convert.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Direct format converters: zarr
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an h5ad file to zarr format
#'
#' Converts an AnnData h5ad file to zarr format. By default uses streaming
#' conversion that copies fields directly without loading into R memory.
#' Both formats follow the same AnnData on-disk specification.
#'
#' @param source Path to input .h5ad file
#' @param dest Path for output .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly without creating
#'   a Seurat intermediate. If FALSE, load as Seurat then save as zarr.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5ADToZarr <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_h5ad_to_zarr(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5AD(source, verbose = verbose)
    writeZarr(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_h5ad_to_zarr <- function(h5ad_path, zarr_path, gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for H5ADToZarr", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for H5ADToZarr", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
  on.exit(h5$close_all())

  compressor <- list(id = "zlib", level = as.integer(gzip))

  # Root group
  .zarr_create_group(zarr_path, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # Helper: stream a DataFrame group (obs or var)
  .stream_h5_df_to_zarr <- function(h5_group, zarr_group_path) {
    if (!inherits(h5_group, "H5Group")) return()

    # Read index
    index <- NULL
    index_name <- "_index"
    if (h5_group$attr_exists("_index")) {
      index_name <- hdf5r::h5attr(h5_group, "_index")
    }
    if (h5_group$exists(index_name)) {
      index <- as.character(h5_group[[index_name]][])
    } else if (h5_group$exists("_index")) {
      index <- as.character(h5_group[["_index"]][])
    } else if (h5_group$exists("index")) {
      index <- as.character(h5_group[["index"]][])
    }
    if (is.null(index)) return()

    # Get column-order
    col_order <- if (h5_group$attr_exists("column-order")) {
      hdf5r::h5attr(h5_group, "column-order")
    } else {
      setdiff(names(h5_group), c("_index", "index", "__categories"))
    }

    # Create group with DataFrame attrs
    .zarr_create_group(zarr_group_path, attrs = list(
      `encoding-type` = "dataframe",
      `encoding-version` = "0.2.0",
      `_index` = "_index",
      `column-order` = col_order
    ))

    # Write index
    .zarr_write_strings(
      dir = file.path(zarr_group_path, "_index"),
      strings = index,
      compressor = compressor,
      attrs = list(`encoding-type` = "string-array", `encoding-version` = "0.2.0")
    )

    # Cache legacy categories
    has_cats <- h5_group$exists("__categories")
    cats <- if (has_cats) h5_group[["__categories"]] else NULL
    cat_names <- if (has_cats) names(cats) else character(0)

    # Stream each column
    for (col in col_order) {
      if (!h5_group$exists(col)) next
      col_obj <- h5_group[[col]]
      col_dir <- file.path(zarr_group_path, col)

      if (inherits(col_obj, "H5Group")) {
        # Modern categorical
        enc_type <- tryCatch(hdf5r::h5attr(col_obj, "encoding-type"),
                             error = function(e) "")
        if (enc_type == "categorical" &&
            col_obj$exists("categories") && col_obj$exists("codes")) {
          codes <- col_obj[["codes"]]$read()
          categories <- as.character(col_obj[["categories"]]$read())
          .zarr_create_group(col_dir, attrs = list(
            `encoding-type` = "categorical",
            `encoding-version` = "0.2.0",
            ordered = FALSE
          ))
          .zarr_write_numeric(
            dir = file.path(col_dir, "codes"),
            data = as.integer(codes),
            dtype = "|i1",
            compressor = compressor
          )
          .zarr_write_strings(
            dir = file.path(col_dir, "categories"),
            strings = categories,
            compressor = compressor,
            attrs = list(`encoding-type` = "string-array",
                         `encoding-version` = "0.2.0")
          )
        }
      } else if (inherits(col_obj, "H5D")) {
        if (has_cats && col %in% cat_names) {
          # Legacy categorical
          codes <- col_obj$read()
          categories <- as.character(cats[[col]]$read())
          .zarr_create_group(col_dir, attrs = list(
            `encoding-type` = "categorical",
            `encoding-version` = "0.2.0",
            ordered = FALSE
          ))
          .zarr_write_numeric(
            dir = file.path(col_dir, "codes"),
            data = as.integer(codes),
            dtype = "|i1",
            compressor = compressor
          )
          .zarr_write_strings(
            dir = file.path(col_dir, "categories"),
            strings = categories,
            compressor = compressor,
            attrs = list(`encoding-type` = "string-array",
                         `encoding-version` = "0.2.0")
          )
        } else {
          vals <- col_obj$read()
          if (is.character(vals)) {
            .zarr_write_strings(
              dir = col_dir, strings = vals, compressor = compressor,
              attrs = list(`encoding-type` = "string-array",
                           `encoding-version` = "0.2.0")
            )
          } else if (is.logical(vals)) {
            .zarr_write_numeric(dir = col_dir, data = as.integer(vals),
                                dtype = "|b1", compressor = compressor)
          } else if (is.integer(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<i4", compressor = compressor)
          } else if (is.numeric(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<f8", compressor = compressor)
          }
        }
      }
    }
  }

  # Helper: stream a matrix (sparse or dense) from h5ad to zarr
  .stream_h5_matrix_to_zarr <- function(h5_obj, zarr_mat_path, transpose) {
    if (inherits(h5_obj, "H5Group") &&
        h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
      # Sparse
      encoding <- tryCatch(hdf5r::h5attr(h5_obj, "encoding-type"),
                           error = function(e) "csr_matrix")
      shape <- tryCatch(hdf5r::h5attr(h5_obj, "shape"), error = function(e) NULL)

      data_vals <- h5_obj[["data"]][]
      indices <- h5_obj[["indices"]][]
      indptr <- h5_obj[["indptr"]][]

      is_csr <- identical(encoding, "csr_matrix")
      is_csc <- identical(encoding, "csc_matrix")

      if (is.null(shape)) {
        if (is_csr) {
          shape <- c(length(indptr) - 1L,
                     if (length(indices) > 0) max(indices) + 1L else 0L)
        } else {
          shape <- c(if (length(indices) > 0) max(indices) + 1L else 0L,
                     length(indptr) - 1L)
        }
      }

      if (transpose && is_csr) {
        # CSR cellsxgenes -> reinterpret as CSC genesxcells -> write back as CSR cellsxgenes
        # Just copy directly -- the data is already in cellsxgenes CSR
        out_encoding <- "csr_matrix"
        out_shape <- as.integer(shape)
      } else if (!transpose) {
        # Keep original encoding (obsp graphs)
        out_encoding <- encoding
        out_shape <- as.integer(shape)
      } else {
        # CSC + transpose: need to re-encode as CSR
        # Build dgCMatrix from CSC, then DeconstructSparseCSR for transposed CSR
        mat <- Matrix::sparseMatrix(
          i = as.integer(indices) + 1L,
          p = as.integer(indptr),
          x = as.numeric(data_vals),
          dims = as.integer(shape),
          repr = "C"
        )
        csr <- DeconstructSparseCSR(mat)
        data_vals <- csr$data
        indices <- csr$indices
        indptr <- csr$indptr
        out_shape <- as.integer(csr$shape)
        out_encoding <- "csr_matrix"
      }

      .zarr_create_group(zarr_mat_path, attrs = list(
        `encoding-type` = out_encoding,
        `encoding-version` = "0.1.0",
        shape = as.list(out_shape)
      ))
      .zarr_parallel_write(
        write_specs = list(
          list(dir = file.path(zarr_mat_path, "data"),
               data = as.numeric(data_vals), dtype = "<f8"),
          list(dir = file.path(zarr_mat_path, "indices"),
               data = as.integer(indices), dtype = "<i4"),
          list(dir = file.path(zarr_mat_path, "indptr"),
               data = as.integer(indptr), dtype = "<i4")
        ),
        compressor = compressor
      )
    } else if (inherits(h5_obj, "H5D")) {
      # Dense -- check dimensionality: 1D datasets need [] not [,]
      ndims <- length(h5_obj$dims)
      mat <- if (ndims == 1L) {
        matrix(h5_obj[], ncol = 1L)
      } else {
        h5_obj[,]
      }
      if (transpose) mat <- t(mat)
      .zarr_write_numeric(
        dir = zarr_mat_path, data = mat, dtype = "<f8",
        compressor = compressor,
        attrs = list(`encoding-type` = "array", `encoding-version` = "0.2.0")
      )
    }
  }

  # 1. obs
  if (verbose) message("Streaming obs...")
  if (h5$exists("obs")) {
    .stream_h5_df_to_zarr(h5[["obs"]], file.path(zarr_path, "obs"))
  } else {
    .zarr_create_group(file.path(zarr_path, "obs"))
  }

  # 2. var
  if (verbose) message("Streaming var...")
  if (h5$exists("var")) {
    .stream_h5_df_to_zarr(h5[["var"]], file.path(zarr_path, "var"))
  } else {
    .zarr_create_group(file.path(zarr_path, "var"))
  }

  # 3. X
  if (h5$exists("X")) {
    if (verbose) message("Streaming X...")
    .stream_h5_matrix_to_zarr(h5[["X"]], file.path(zarr_path, "X"),
                               transpose = FALSE)
  }

  # 4. layers
  .zarr_create_group(file.path(zarr_path, "layers"))
  if (h5$exists("layers")) {
    for (ln in names(h5[["layers"]])) {
      if (verbose) message("  Streaming layer: ", ln)
      tryCatch({
        .stream_h5_matrix_to_zarr(
          h5[["layers"]][[ln]],
          file.path(zarr_path, "layers", ln),
          transpose = FALSE
        )
      }, error = function(e) {
        if (verbose) warning("Could not stream layer ", ln, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 5. obsm
  .zarr_create_group(file.path(zarr_path, "obsm"))
  if (h5$exists("obsm") && inherits(h5[["obsm"]], "H5Group")) {
    for (rn in names(h5[["obsm"]])) {
      if (verbose) message("  Streaming obsm: ", rn)
      tryCatch({
        item <- h5[["obsm"]][[rn]]
        if (inherits(item, "H5Group")) {
          # Sparse embedding (rare)
          .stream_h5_matrix_to_zarr(item, file.path(zarr_path, "obsm", rn),
                                     transpose = FALSE)
        } else {
          # Dense: read and write as-is (already cells x dims in h5ad)
          ndims_item <- length(item$dims)
          emb <- if (ndims_item == 1L) {
            matrix(item[], ncol = 1L)
          } else {
            item[,]
          }
          # h5ad stores (n_obs, n_dims) but hdf5r may return transposed
          if (ndims_item > 1L && ncol(emb) > nrow(emb) && item$dims[1] > item$dims[2]) {
            emb <- t(emb)
          }
          .zarr_write_numeric(
            dir = file.path(zarr_path, "obsm", rn),
            data = emb, dtype = "<f8", compressor = compressor
          )
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream obsm ", rn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 6. obsp
  .zarr_create_group(file.path(zarr_path, "obsp"))
  if (h5$exists("obsp") && inherits(h5[["obsp"]], "H5Group")) {
    for (gn in names(h5[["obsp"]])) {
      if (verbose) message("  Streaming obsp: ", gn)
      tryCatch({
        .stream_h5_matrix_to_zarr(h5[["obsp"]][[gn]],
                                   file.path(zarr_path, "obsp", gn),
                                   transpose = FALSE)
      }, error = function(e) {
        if (verbose) warning("Could not stream obsp ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 7. uns
  .zarr_create_group(file.path(zarr_path, "uns"))
  if (h5$exists("uns") && inherits(h5[["uns"]], "H5Group")) {
    for (un in names(h5[["uns"]])) {
      tryCatch({
        item <- h5[["uns"]][[un]]
        if (inherits(item, "H5D")) {
          vals <- item[]
          if (is.character(vals)) {
            .zarr_write_strings(
              dir = file.path(zarr_path, "uns", un),
              strings = vals, compressor = compressor
            )
          } else if (is.numeric(vals) || is.logical(vals)) {
            .zarr_write_numeric(
              dir = file.path(zarr_path, "uns", un),
              data = vals, compressor = compressor
            )
          }
        }
        # Skip complex groups silently
      }, error = function(e) {
        if (verbose) warning("Could not stream uns ", un, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 8. varp
  if (h5$exists("varp") && inherits(h5[["varp"]], "H5Group")) {
    if (verbose) message("Streaming varp...")
    .zarr_create_group(file.path(zarr_path, "varp"))
    for (vn in names(h5[["varp"]])) {
      tryCatch({
        .stream_h5_matrix_to_zarr(h5[["varp"]][[vn]],
                                   file.path(zarr_path, "varp", vn),
                                   transpose = FALSE)
      }, error = function(e) {
        if (verbose) warning("Could not stream varp ", vn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 9. raw (if present, write as raw/ group with X and var)
  if (h5$exists("raw")) {
    if (verbose) message("Streaming raw...")
    raw_grp <- h5[["raw"]]
    .zarr_create_group(file.path(zarr_path, "raw"))
    if (raw_grp$exists("X")) {
      .stream_h5_matrix_to_zarr(raw_grp[["X"]],
                                 file.path(zarr_path, "raw", "X"),
                                 transpose = FALSE)
    }
    if (raw_grp$exists("var")) {
      .stream_h5_df_to_zarr(raw_grp[["var"]], file.path(zarr_path, "raw", "var"))
    }
  }

  if (verbose) message("Streaming h5ad -> zarr complete: ", zarr_path)
  invisible(zarr_path)
}

#' Convert a zarr store to h5ad format
#'
#' Converts an AnnData zarr store to h5ad (HDF5) format. By default uses
#' streaming conversion that copies fields directly without loading into R
#' memory.
#'
#' @param source Path to input .zarr directory
#' @param dest Path for output .h5ad file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly without creating
#'   a Seurat intermediate. If FALSE, load as Seurat then save as h5ad.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToH5AD <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!dir.exists(source))
    stop("Zarr store not found: ", source, call. = FALSE)
  if (file.exists(dest) || dir.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_zarr_to_h5ad(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readZarr(source, verbose = verbose)
    writeH5AD(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_zarr_to_h5ad <- function(zarr_path, h5ad_path, gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for ZarrToH5AD", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for ZarrToH5AD", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5ad_path, mode = "w")
  on.exit(h5$close_all())

  # Root attrs
  h5$create_attr(attr_name = "encoding-type", robj = "anndata",
                 dtype = CachedGuessDType("anndata"), space = ScalarSpace())
  h5$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                 dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())

  # Helper: stream a zarr DataFrame to h5ad
  .stream_zarr_df_to_h5 <- function(store_path, zarr_group, h5_parent, h5_name) {
    grp_path <- zarr_group
    if (.zarr_node_type(store_path, grp_path) != "group") {
      # Write empty group
      g <- h5_parent$create_group(h5_name)
      g$create_attr(attr_name = "encoding-type", robj = "dataframe",
                    dtype = CachedGuessDType("dataframe"), space = ScalarSpace())
      g$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                    dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
      g$create_attr(attr_name = "_index", robj = "_index",
                    dtype = CachedGuessDType("_index"), space = ScalarSpace())
      g$create_attr(attr_name = "column-order", robj = character(0),
                    dtype = CachedUtf8Type())
      return()
    }

    attrs <- .zarr_read_attrs(store_path, grp_path)
    index <- .zarr_read_anndata_index(store_path, zarr_group)
    col_order <- attrs[["column-order"]]
    if (is.null(col_order)) {
      children <- .zarr_list_children(store_path, grp_path)
      col_order <- setdiff(children, c("_index", "index"))
    }

    # Build a data.frame for WriteDFGroup
    df <- data.frame(row.names = if (!is.null(index)) index else seq_along(index))
    for (col in col_order) {
      col_path <- file.path(grp_path, col)
      if (.zarr_node_type(store_path, col_path) == "missing") next
      tryCatch({
        vals <- .zarr_read_anndata_column(store_path, col_path)
        if (!is.null(vals)) df[[col]] <- vals
      }, error = function(e) {
        if (verbose) warning("Could not read column ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }

    idx <- if (!is.null(index)) index else rownames(df)
    WriteDFGroup(h5_parent, h5_name, df, idx, gzip = gzip)
  }

  # Helper: stream a zarr matrix to h5ad
  .stream_zarr_matrix_to_h5 <- function(store_path, zarr_mat_path,
                                         h5_parent, h5_name) {
    node_type <- .zarr_node_type(store_path, zarr_mat_path)
    if (node_type == "missing") return()

    if (node_type == "group") {
      attrs <- .zarr_read_attrs(store_path, zarr_mat_path)
      enc <- attrs[["encoding-type"]] %||% ""

      if (enc %in% c("csr_matrix", "csc_matrix")) {
        data_vals <- .zarr_read_numeric(store_path,
                                         file.path(zarr_mat_path, "data"))
        indices <- .zarr_read_numeric(store_path,
                                       file.path(zarr_mat_path, "indices"))
        indptr <- .zarr_read_numeric(store_path,
                                      file.path(zarr_mat_path, "indptr"))
        shape <- attrs[["shape"]]
        if (is.list(shape)) shape <- unlist(shape)

        chunk_size <- 65536L
        grp <- h5_parent$create_group(h5_name)
        grp$create_dataset("data", robj = as.numeric(data_vals),
                           chunk_dims = min(length(data_vals), chunk_size),
                           gzip_level = gzip)
        grp$create_dataset("indices", robj = as.integer(indices),
                           chunk_dims = min(length(indices), chunk_size),
                           gzip_level = gzip)
        grp$create_dataset("indptr", robj = as.integer(indptr),
                           chunk_dims = min(length(indptr), chunk_size),
                           gzip_level = gzip)
        grp$create_attr(attr_name = "encoding-type", robj = enc,
                        dtype = CachedGuessDType(enc), space = ScalarSpace())
        grp$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                        dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())
        grp$create_attr(attr_name = "shape", robj = as.integer(shape),
                        dtype = GuessDType(as.integer(shape)))
      }
    } else if (node_type == "array") {
      mat <- .zarr_read_numeric(store_path, zarr_mat_path)
      if (is.matrix(mat)) {
        h5_parent$create_dataset(h5_name, robj = mat,
                                 chunk_dims = dim(mat), gzip_level = gzip)
      } else {
        h5_parent$create_dataset(h5_name, robj = mat,
                                 chunk_dims = length(mat), gzip_level = gzip)
      }
      ds <- h5_parent[[h5_name]]
      ds$create_attr(attr_name = "encoding-type", robj = "array",
                     dtype = CachedGuessDType("array"), space = ScalarSpace())
      ds$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                     dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
    }
  }

  # 1. obs
  if (verbose) message("Streaming obs...")
  .stream_zarr_df_to_h5(zarr_path, "obs", h5, "obs")

  # 2. var
  if (verbose) message("Streaming var...")
  .stream_zarr_df_to_h5(zarr_path, "var", h5, "var")

  # 3. X
  if (.zarr_node_type(zarr_path, "X") != "missing") {
    if (verbose) message("Streaming X...")
    .stream_zarr_matrix_to_h5(zarr_path, "X", h5, "X")
  }

  # 4. layers
  if (.zarr_node_type(zarr_path, "layers") == "group") {
    layer_names <- .zarr_list_children(zarr_path, "layers")
    if (length(layer_names) > 0) {
      layers_grp <- h5$create_group("layers")
      for (ln in layer_names) {
        if (verbose) message("  Streaming layer: ", ln)
        tryCatch({
          .stream_zarr_matrix_to_h5(zarr_path, file.path("layers", ln),
                                     layers_grp, ln)
        }, error = function(e) {
          if (verbose) warning("Could not stream layer ", ln, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 5. obsm
  if (.zarr_node_type(zarr_path, "obsm") == "group") {
    obsm_names <- .zarr_list_children(zarr_path, "obsm")
    if (length(obsm_names) > 0) {
      obsm_grp <- h5$create_group("obsm")
      for (rn in obsm_names) {
        if (verbose) message("  Streaming obsm: ", rn)
        tryCatch({
          node <- .zarr_node_type(zarr_path, file.path("obsm", rn))
          if (node == "group") {
            .stream_zarr_matrix_to_h5(zarr_path, file.path("obsm", rn),
                                       obsm_grp, rn)
          } else if (node == "array") {
            emb <- .zarr_read_numeric(zarr_path, file.path("obsm", rn))
            if (is.matrix(emb)) {
              obsm_grp$create_dataset(rn, robj = emb,
                                      chunk_dims = dim(emb), gzip_level = gzip)
            }
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream obsm ", rn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 6. obsp
  if (.zarr_node_type(zarr_path, "obsp") == "group") {
    obsp_names <- .zarr_list_children(zarr_path, "obsp")
    if (length(obsp_names) > 0) {
      obsp_grp <- h5$create_group("obsp")
      for (gn in obsp_names) {
        if (verbose) message("  Streaming obsp: ", gn)
        tryCatch({
          .stream_zarr_matrix_to_h5(zarr_path, file.path("obsp", gn),
                                     obsp_grp, gn)
        }, error = function(e) {
          if (verbose) warning("Could not stream obsp ", gn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 7. uns
  if (.zarr_node_type(zarr_path, "uns") == "group") {
    uns_names <- .zarr_list_children(zarr_path, "uns")
    if (length(uns_names) > 0) {
      uns_grp <- h5$create_group("uns")
      for (un in uns_names) {
        tryCatch({
          node <- .zarr_node_type(zarr_path, file.path("uns", un))
          if (node == "array") {
            meta <- .zarr_read_json(file.path(zarr_path, "uns", un, ".zarray"))
            if (!is.null(meta$dtype) && meta$dtype == "|O") {
              vals <- .zarr_read_strings(zarr_path, file.path("uns", un))
              uns_grp$create_dataset(un, robj = vals,
                                     dtype = CachedUtf8Type(),
                                     chunk_dims = length(vals), gzip_level = gzip)
            } else {
              vals <- .zarr_read_numeric(zarr_path, file.path("uns", un))
              uns_grp$create_dataset(un, robj = vals,
                                     chunk_dims = length(vals), gzip_level = gzip)
            }
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream uns ", un, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 8. varp
  if (.zarr_node_type(zarr_path, "varp") == "group") {
    varp_names <- .zarr_list_children(zarr_path, "varp")
    if (length(varp_names) > 0) {
      varp_grp <- h5$create_group("varp")
      for (vn in varp_names) {
        if (verbose) message("  Streaming varp: ", vn)
        tryCatch({
          .stream_zarr_matrix_to_h5(zarr_path, file.path("varp", vn),
                                     varp_grp, vn)
        }, error = function(e) {
          if (verbose) warning("Could not stream varp ", vn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 9. raw
  if (.zarr_node_type(zarr_path, "raw") == "group") {
    if (verbose) message("Streaming raw...")
    raw_grp <- h5$create_group("raw")
    if (.zarr_node_type(zarr_path, "raw/X") != "missing") {
      .stream_zarr_matrix_to_h5(zarr_path, "raw/X", raw_grp, "X")
    }
    if (.zarr_node_type(zarr_path, "raw/var") == "group") {
      .stream_zarr_df_to_h5(zarr_path, "raw/var", raw_grp, "var")
    }
  }

  if (verbose) message("Streaming zarr -> h5ad complete: ", h5ad_path)
  invisible(h5ad_path)
}

#' Convert an h5seurat file to zarr (AnnData) format
#'
#' Converts an h5seurat file to AnnData zarr format, translating the h5seurat
#' layout (levels/values factors, genes x cells CSC, reductions as
#' cell.embeddings) to AnnData conventions (categories/codes, cells x genes
#' CSR, obsm arrays). By default uses streaming that avoids a Seurat
#' intermediate.
#'
#' @param source Path to input .h5seurat file
#' @param dest Path for output .zarr directory
#' @param assay Assay name (default "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as zarr.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5SeuratToZarr <- function(source, dest, assay = "RNA", overwrite = FALSE,
                            gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_h5seurat_to_zarr(source, dest, assay = assay, gzip = gzip,
                              verbose = verbose)
  } else {
    obj <- readH5Seurat(source, verbose = verbose)
    writeZarr(obj, dest, assay = assay, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_h5seurat_to_zarr <- function(h5seurat_path, zarr_path, assay = "RNA",
                                      gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for H5SeuratToZarr", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for H5SeuratToZarr", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5seurat_path, mode = "r")
  on.exit(h5$close_all())

  compressor <- list(id = "zlib", level = as.integer(gzip))

  # Detect active assay
  if (h5$attr_exists("active.assay")) {
    assay <- hdf5r::h5attr(h5, "active.assay")
  }

  # Root zarr group
  .zarr_create_group(zarr_path, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # Read cell names (SeuratDisk may store as Nx1 2D dataset)
  cell_names <- if (h5$exists("cell.names")) {
    cn_ds <- h5[["cell.names"]]
    if (length(cn_ds$dims) == 2L) as.character(cn_ds[, 1]) else as.character(cn_ds[])
  } else {
    NULL
  }

  # Read feature names (Seurat v5 then v4 paths)
  assay_path <- paste0("assays/", assay)
  feature_names <- NULL
  if (h5$exists(assay_path)) {
    ag <- h5[[assay_path]]
    if (ag$exists("features")) {
      feat_ds <- ag[["features"]]
      # SeuratDisk may store features as Nx1 (2D), handle both 1D and 2D
      feature_names <- if (length(feat_ds$dims) == 2L) {
        as.character(feat_ds[, 1])
      } else {
        as.character(feat_ds[])
      }
    }
  }

  n_cells <- length(cell_names)
  n_features <- length(feature_names)

  # --- obs (meta.data -> AnnData DataFrame) ---
  if (verbose) message("Streaming obs (meta.data)...")
  if (h5$exists("meta.data")) {
    md <- h5[["meta.data"]]
    md_cols <- character(0)
    # Get column order from colnames attr
    if (md$attr_exists("colnames")) {
      md_cols <- hdf5r::h5attr(md, "colnames")
    } else {
      md_cols <- names(md)
    }

    .zarr_create_group(file.path(zarr_path, "obs"), attrs = list(
      `encoding-type` = "dataframe",
      `encoding-version` = "0.2.0",
      `_index` = "_index",
      `column-order` = md_cols
    ))

    # Write _index
    if (!is.null(cell_names)) {
      .zarr_write_strings(
        dir = file.path(zarr_path, "obs", "_index"),
        strings = cell_names,
        compressor = compressor,
        attrs = list(`encoding-type` = "string-array",
                     `encoding-version` = "0.2.0")
      )
    }

    # Write each column
    for (col in md_cols) {
      if (!md$exists(col)) next
      col_obj <- md[[col]]
      col_dir <- file.path(zarr_path, "obs", col)

      tryCatch({
        if (inherits(col_obj, "H5Group")) {
          # Factor: levels/values (1-based) -> categories/codes (0-based)
          if (col_obj$exists("levels") && col_obj$exists("values")) {
            levs <- as.character(col_obj[["levels"]][])
            vals <- col_obj[["values"]][]
            codes <- as.integer(vals) - 1L  # 1-based -> 0-based

            .zarr_create_group(col_dir, attrs = list(
              `encoding-type` = "categorical",
              `encoding-version` = "0.2.0",
              ordered = FALSE
            ))
            .zarr_write_numeric(
              dir = file.path(col_dir, "codes"),
              data = codes, dtype = "|i1", compressor = compressor
            )
            .zarr_write_strings(
              dir = file.path(col_dir, "categories"),
              strings = levs, compressor = compressor,
              attrs = list(`encoding-type` = "string-array",
                           `encoding-version` = "0.2.0")
            )
          }
        } else if (inherits(col_obj, "H5D")) {
          vals <- col_obj$read()
          if (is.character(vals)) {
            .zarr_write_strings(
              dir = col_dir, strings = vals, compressor = compressor,
              attrs = list(`encoding-type` = "string-array",
                           `encoding-version` = "0.2.0")
            )
          } else if (is.logical(vals)) {
            .zarr_write_numeric(dir = col_dir, data = as.integer(vals),
                                dtype = "|b1", compressor = compressor)
          } else if (is.integer(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<i4", compressor = compressor)
          } else if (is.numeric(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<f8", compressor = compressor)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream obs column ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  } else {
    .zarr_create_group(file.path(zarr_path, "obs"), attrs = list(
      `encoding-type` = "dataframe", `encoding-version` = "0.2.0",
      `_index` = "_index", `column-order` = character(0)
    ))
    if (!is.null(cell_names)) {
      .zarr_write_strings(
        dir = file.path(zarr_path, "obs", "_index"),
        strings = cell_names, compressor = compressor,
        attrs = list(`encoding-type` = "string-array",
                     `encoding-version` = "0.2.0")
      )
    }
  }

  # --- var (features -> AnnData DataFrame) ---
  if (verbose) message("Streaming var...")
  var_cols <- character(0)
  if (h5$exists("var")) {
    var_grp <- h5[["var"]]
    if (var_grp$attr_exists("colnames")) {
      var_cols <- hdf5r::h5attr(var_grp, "colnames")
    } else {
      var_cols <- names(var_grp)
    }
  }
  .zarr_create_group(file.path(zarr_path, "var"), attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = var_cols
  ))
  if (!is.null(feature_names)) {
    .zarr_write_strings(
      dir = file.path(zarr_path, "var", "_index"),
      strings = feature_names, compressor = compressor,
      attrs = list(`encoding-type` = "string-array",
                   `encoding-version` = "0.2.0")
    )
  }
  # Write var metadata columns
  if (h5$exists("var") && length(var_cols) > 0) {
    var_grp <- h5[["var"]]
    for (col in var_cols) {
      if (!var_grp$exists(col)) next
      col_obj <- var_grp[[col]]
      col_dir <- file.path(zarr_path, "var", col)
      tryCatch({
        if (inherits(col_obj, "H5Group") &&
            col_obj$exists("levels") && col_obj$exists("values")) {
          levs <- as.character(col_obj[["levels"]][])
          vals <- col_obj[["values"]][]
          codes <- as.integer(vals) - 1L
          .zarr_create_group(col_dir, attrs = list(
            `encoding-type` = "categorical", `encoding-version` = "0.2.0",
            ordered = FALSE))
          .zarr_write_numeric(dir = file.path(col_dir, "codes"),
                              data = codes, dtype = "|i1", compressor = compressor)
          .zarr_write_strings(dir = file.path(col_dir, "categories"),
                              strings = levs, compressor = compressor,
                              attrs = list(`encoding-type` = "string-array",
                                           `encoding-version` = "0.2.0"))
        } else if (inherits(col_obj, "H5D")) {
          vals <- col_obj$read()
          if (is.character(vals)) {
            .zarr_write_strings(dir = col_dir, strings = vals,
                                compressor = compressor,
                                attrs = list(`encoding-type` = "string-array",
                                             `encoding-version` = "0.2.0"))
          } else if (is.logical(vals)) {
            .zarr_write_numeric(dir = col_dir, data = as.integer(vals),
                                dtype = "|b1", compressor = compressor)
          } else if (is.integer(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<i4", compressor = compressor)
          } else if (is.numeric(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<f8", compressor = compressor)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream var column ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- X (expression matrix: genesxcells CSC -> cellsxgenes CSR) ---
  # Try v5 path first, then v4
  x_path <- NULL
  for (p in c(paste0(assay_path, "/layers/data"),
              paste0(assay_path, "/data"))) {
    if (h5$exists(p)) { x_path <- p; break }
  }
  if (!is.null(x_path)) {
    if (verbose) message("Streaming X (", x_path, ")...")
    x_obj <- h5[[x_path]]
    .stream_h5seurat_matrix_to_zarr(x_obj, file.path(zarr_path, "X"),
                                     compressor, transpose = TRUE)
  }

  # --- layers (counts -> raw/X if separate from data) ---
  counts_path <- NULL
  for (p in c(paste0(assay_path, "/layers/counts"),
              paste0(assay_path, "/counts"))) {
    if (h5$exists(p)) { counts_path <- p; break }
  }

  # If no data layer exists but counts does, also write counts as X
  # so that anndata.read_zarr() puts data in adata.X (not just adata.raw.X)
  if (is.null(x_path) && !is.null(counts_path)) {
    if (verbose) message("Streaming X (from counts: ", counts_path, ")...")
    counts_obj <- h5[[counts_path]]
    .stream_h5seurat_matrix_to_zarr(counts_obj, file.path(zarr_path, "X"),
                                     compressor, transpose = TRUE)
  }

  .zarr_create_group(file.path(zarr_path, "layers"))
  if (!is.null(counts_path) && !identical(counts_path, x_path)) {
    if (verbose) message("Streaming raw/X (", counts_path, ")...")
    .zarr_create_group(file.path(zarr_path, "raw"))
    counts_obj <- h5[[counts_path]]
    .stream_h5seurat_matrix_to_zarr(counts_obj, file.path(zarr_path, "raw", "X"),
                                     compressor, transpose = TRUE)
    # Write raw/var with just _index
    .zarr_create_group(file.path(zarr_path, "raw", "var"), attrs = list(
      `encoding-type` = "dataframe", `encoding-version` = "0.2.0",
      `_index` = "_index", `column-order` = character(0)
    ))
    if (!is.null(feature_names)) {
      .zarr_write_strings(
        dir = file.path(zarr_path, "raw", "var", "_index"),
        strings = feature_names, compressor = compressor,
        attrs = list(`encoding-type` = "string-array",
                     `encoding-version` = "0.2.0")
      )
    }
  }

  # --- obsm (reductions -> obsm/X_{name}) ---
  .zarr_create_group(file.path(zarr_path, "obsm"))
  if (h5$exists("reductions")) {
    reduc_grp <- h5[["reductions"]]
    for (rn in names(reduc_grp)) {
      if (verbose) message("  Streaming obsm: X_", rn)
      tryCatch({
        r_obj <- reduc_grp[[rn]]
        if (r_obj$exists("cell.embeddings")) {
          emb <- r_obj[["cell.embeddings"]]
          if (inherits(emb, "H5D")) {
            ndims_emb <- length(emb$dims)
            mat <- if (ndims_emb == 1L) {
              matrix(emb[], ncol = 1L)
            } else {
              emb[,]
            }
            # hdf5r reads h5seurat cell.embeddings as [n_cells, n_comp]
            # zarr AnnData also needs [n_cells, n_comp] -- no transpose needed
            # Only fix if dimensions are unexpectedly swapped
            if (nrow(mat) != n_cells && ncol(mat) == n_cells) {
              mat <- t(mat)
            }
            .zarr_write_numeric(
              dir = file.path(zarr_path, "obsm", paste0("X_", rn)),
              data = mat, dtype = "<f8", compressor = compressor
            )
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream reduction ", rn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- obsp (graphs -> obsp/{name}) ---
  .zarr_create_group(file.path(zarr_path, "obsp"))
  if (h5$exists("graphs")) {
    graph_grp <- h5[["graphs"]]
    for (gn in names(graph_grp)) {
      if (verbose) message("  Streaming obsp: ", gn)
      tryCatch({
        g_obj <- graph_grp[[gn]]
        if (inherits(g_obj, "H5Group") &&
            g_obj$exists("data") && g_obj$exists("indices") && g_obj$exists("indptr")) {
          # h5seurat graphs are CSC (dgCMatrix layout)
          data_vals <- g_obj[["data"]][]
          indices <- g_obj[["indices"]][]
          indptr <- g_obj[["indptr"]][]
          shape <- if (g_obj$attr_exists("dims")) {
            hdf5r::h5attr(g_obj, "dims")
          } else {
            c(length(indptr) - 1L,
              if (length(indices) > 0) max(indices) + 1L else 0L)
          }

          .zarr_create_group(file.path(zarr_path, "obsp", gn), attrs = list(
            `encoding-type` = "csc_matrix",
            `encoding-version` = "0.1.0",
            shape = as.list(as.integer(shape))
          ))
          .zarr_parallel_write(
            write_specs = list(
              list(dir = file.path(zarr_path, "obsp", gn, "data"),
                   data = as.numeric(data_vals), dtype = "<f8"),
              list(dir = file.path(zarr_path, "obsp", gn, "indices"),
                   data = as.integer(indices), dtype = "<i4"),
              list(dir = file.path(zarr_path, "obsp", gn, "indptr"),
                   data = as.integer(indptr), dtype = "<i4")
            ),
            compressor = compressor
          )
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream graph ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- uns (misc -> uns) ---
  .zarr_create_group(file.path(zarr_path, "uns"))
  if (h5$exists("misc")) {
    misc_grp <- h5[["misc"]]
    for (un in names(misc_grp)) {
      tryCatch({
        item <- misc_grp[[un]]
        if (inherits(item, "H5D")) {
          vals <- item[]
          if (is.character(vals)) {
            .zarr_write_strings(dir = file.path(zarr_path, "uns", un),
                                strings = vals, compressor = compressor)
          } else if (is.numeric(vals) || is.logical(vals)) {
            .zarr_write_numeric(dir = file.path(zarr_path, "uns", un),
                                data = vals, compressor = compressor)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream misc ", un, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  if (verbose) message("Streaming h5seurat -> zarr complete: ", zarr_path)
  invisible(zarr_path)
}

#' Stream h5seurat sparse matrix to zarr
#'
#' h5seurat stores CSC (genes x cells). If transpose=TRUE, reinterprets as
#' CSR (cells x genes) using the zero-copy CSC<->CSR duality.
#'
#' @keywords internal
#' @noRd
.stream_h5seurat_matrix_to_zarr <- function(h5_obj, zarr_mat_path,
                                             compressor, transpose = TRUE) {
  if (inherits(h5_obj, "H5Group") &&
      h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
    # Sparse CSC: i=row indices (0-based), p=col pointers, x=values
    data_vals <- h5_obj[["data"]][]
    indices <- h5_obj[["indices"]][]
    indptr <- h5_obj[["indptr"]][]
    shape <- if (h5_obj$attr_exists("dims")) {
      hdf5r::h5attr(h5_obj, "dims")
    } else {
      c(if (length(indices) > 0) max(indices) + 1L else 0L,
        length(indptr) - 1L)
    }

    if (transpose) {
      # CSC of (genes x cells) == CSR of (cells x genes) -- zero-copy
      out_shape <- as.integer(rev(shape))  # [n_cells, n_genes]
      .zarr_create_group(zarr_mat_path, attrs = list(
        `encoding-type` = "csr_matrix",
        `encoding-version` = "0.1.0",
        shape = as.list(out_shape)
      ))
    } else {
      out_shape <- as.integer(shape)
      .zarr_create_group(zarr_mat_path, attrs = list(
        `encoding-type` = "csc_matrix",
        `encoding-version` = "0.1.0",
        shape = as.list(out_shape)
      ))
    }
    .zarr_parallel_write(
      write_specs = list(
        list(dir = file.path(zarr_mat_path, "data"),
             data = as.numeric(data_vals), dtype = "<f8"),
        list(dir = file.path(zarr_mat_path, "indices"),
             data = as.integer(indices), dtype = "<i4"),
        list(dir = file.path(zarr_mat_path, "indptr"),
             data = as.integer(indptr), dtype = "<i4")
      ),
      compressor = compressor
    )
  } else if (inherits(h5_obj, "H5D")) {
    ndims <- length(h5_obj$dims)
    if (ndims == 1L) {
      # 1D dense: likely an embedding or a single feature.
      # Cannot reliably determine (rows, cols) from a 1D vector alone,
      # so skip rather than write a wrong-shape matrix.
      warning("Skipping 1D dense dataset in H5SeuratToZarr: ",
              "cannot determine matrix dimensions from 1D data")
    } else {
      mat <- h5_obj[,]
      if (transpose) mat <- t(mat)
      .zarr_write_numeric(
        dir = zarr_mat_path, data = mat, dtype = "<f8",
        compressor = compressor,
        attrs = list(`encoding-type` = "array", `encoding-version` = "0.2.0")
      )
    }
  }
}

#' Convert a zarr (AnnData) store to h5seurat format
#'
#' Converts an AnnData zarr store to h5seurat format, translating AnnData
#' conventions (categories/codes 0-based, cells x genes CSR, obsm arrays)
#' to h5seurat layout (levels/values 1-based, genes x cells CSC,
#' reductions/cell.embeddings). By default uses streaming that avoids a
#' Seurat intermediate.
#'
#' @param source Path to input .zarr directory
#' @param dest Path for output .h5seurat file
#' @param assay Assay name (default "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as h5seurat.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToH5Seurat <- function(source, dest, assay = "RNA", overwrite = FALSE,
                            gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!dir.exists(source))
    stop("Zarr store not found: ", source, call. = FALSE)
  if (file.exists(dest) || dir.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_zarr_to_h5seurat(source, dest, assay = assay, gzip = gzip,
                              verbose = verbose)
  } else {
    obj <- readZarr(source, verbose = verbose)
    writeH5Seurat(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_zarr_to_h5seurat <- function(zarr_path, h5seurat_path, assay = "RNA",
                                      gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for ZarrToH5Seurat", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for ZarrToH5Seurat", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5seurat_path, mode = "w")
  on.exit(h5$close_all())

  # Read cell and feature names from zarr
  cell_names <- .zarr_read_anndata_index(zarr_path, "obs")
  feature_names <- .zarr_read_anndata_index(zarr_path, "var")
  n_cells <- length(cell_names)
  n_features <- length(feature_names)

  # Root attributes
  h5$create_attr(attr_name = "active.assay", robj = assay,
                 dtype = CachedGuessDType(assay), space = ScalarSpace())
  h5$create_attr(attr_name = "project", robj = "scConvert",
                 dtype = CachedGuessDType("scConvert"), space = ScalarSpace())
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("SeuratObject")),
    error = function(e) "5.0.0"
  )
  h5$create_attr(attr_name = "version", robj = pkg_ver,
                 dtype = CachedGuessDType(pkg_ver), space = ScalarSpace())

  # --- cell.names ---
  if (!is.null(cell_names)) {
    h5$create_dataset("cell.names", robj = cell_names,
                      dtype = CachedUtf8Type(),
                      chunk_dims = length(cell_names), gzip_level = gzip)
  }

  # --- meta.data (obs -> meta.data with levels/values factors) ---
  if (verbose) message("Streaming meta.data (obs)...")
  obs_attrs <- .zarr_read_attrs(zarr_path, "obs")
  obs_cols <- obs_attrs[["column-order"]]
  if (is.null(obs_cols)) {
    obs_cols <- setdiff(.zarr_list_children(zarr_path, "obs"),
                        c("_index", "index"))
  }

  md_grp <- h5$create_group("meta.data")
  col_written <- character(0)
  for (col in obs_cols) {
    col_path <- file.path("obs", col)
    if (.zarr_node_type(zarr_path, col_path) == "missing") next
    tryCatch({
      node_type <- .zarr_node_type(zarr_path, col_path)
      if (node_type == "group") {
        attrs <- .zarr_read_attrs(zarr_path, col_path)
        enc <- attrs[["encoding-type"]] %||% ""
        if (enc == "categorical") {
          codes <- .zarr_read_numeric(zarr_path, file.path(col_path, "codes"))
          # Try reading categories as strings first
          cats_path <- file.path(col_path, "categories")
          cats_meta <- .zarr_read_json(file.path(zarr_path, cats_path, ".zarray"))
          cats <- if (!is.null(cats_meta$dtype) && cats_meta$dtype == "|O") {
            .zarr_read_strings(zarr_path, cats_path)
          } else {
            tryCatch(.zarr_read_strings(zarr_path, cats_path),
                     error = function(e) {
                       as.character(.zarr_read_numeric(zarr_path, cats_path))
                     })
          }
          # Write as h5seurat factor: levels + values (1-based)
          fac_grp <- md_grp$create_group(col)
          fac_grp$create_dataset("levels", robj = cats,
                                 dtype = CachedUtf8Type(),
                                 chunk_dims = length(cats), gzip_level = gzip)
          values_1based <- as.integer(codes) + 1L
          fac_grp$create_dataset("values", robj = values_1based,
                                 chunk_dims = length(values_1based),
                                 gzip_level = gzip)
          col_written <- c(col_written, col)
        }
      } else if (node_type == "array") {
        meta <- .zarr_read_json(file.path(zarr_path, col_path, ".zarray"))
        if (!is.null(meta$dtype) && meta$dtype == "|O") {
          vals <- .zarr_read_strings(zarr_path, col_path)
          md_grp$create_dataset(col, robj = vals, dtype = CachedUtf8Type(),
                                chunk_dims = length(vals), gzip_level = gzip)
        } else {
          vals <- .zarr_read_numeric(zarr_path, col_path)
          if (is.logical(vals)) vals <- as.integer(vals)
          md_grp$create_dataset(col, robj = vals,
                                chunk_dims = length(vals), gzip_level = gzip)
        }
        col_written <- c(col_written, col)
      }
    }, error = function(e) {
      if (verbose) warning("Could not stream obs column ", col, ": ",
                           e$message, immediate. = TRUE)
    })
  }
  if (length(col_written) > 0) {
    md_grp$create_attr(attr_name = "colnames", robj = col_written,
                       dtype = GuessDType(col_written))
  }

  # --- assays/{assay} ---
  assay_grp <- h5$create_group("assays")
  a_grp <- assay_grp$create_group(assay)
  a_grp$create_attr(attr_name = "key", robj = paste0(tolower(assay), "_"),
                    dtype = CachedGuessDType(paste0(tolower(assay), "_")),
                    space = ScalarSpace())
  a_grp$create_attr(attr_name = "s4class", robj = "SeuratObject::Assay5",
                    dtype = CachedGuessDType("SeuratObject::Assay5"),
                    space = ScalarSpace())

  # features
  if (!is.null(feature_names)) {
    a_grp$create_dataset("features", robj = feature_names,
                         dtype = CachedUtf8Type(),
                         chunk_dims = length(feature_names), gzip_level = gzip)
  }

  # --- X -> layers/data (transpose CSR cellsxgenes -> CSC genesxcells) ---
  layers_grp <- a_grp$create_group("layers")
  if (.zarr_node_type(zarr_path, "X") != "missing") {
    if (verbose) message("Streaming X -> layers/data...")
    .stream_zarr_matrix_to_h5seurat(zarr_path, "X", layers_grp, "data",
                                     gzip, transpose = TRUE)
  }

  # --- raw/X or layers/counts -> layers/counts ---
  counts_src <- NULL
  if (.zarr_node_type(zarr_path, "raw/X") != "missing") {
    counts_src <- "raw/X"
  } else if (.zarr_node_type(zarr_path, "layers") == "group") {
    lc <- .zarr_list_children(zarr_path, "layers")
    if ("counts" %in% lc) counts_src <- "layers/counts"
  }
  if (!is.null(counts_src)) {
    if (verbose) message("Streaming ", counts_src, " -> layers/counts...")
    .stream_zarr_matrix_to_h5seurat(zarr_path, counts_src, layers_grp,
                                     "counts", gzip, transpose = TRUE)
  }

  # --- reductions (obsm -> reductions) ---
  reductions_grp <- h5$create_group("reductions")
  if (.zarr_node_type(zarr_path, "obsm") == "group") {
    for (rn in .zarr_list_children(zarr_path, "obsm")) {
      clean_name <- gsub("^X_", "", rn)
      if (verbose) message("  Streaming reduction: ", clean_name)
      tryCatch({
        r_grp <- reductions_grp$create_group(clean_name)
        r_grp$create_attr(attr_name = "active.assay", robj = assay,
                          dtype = CachedGuessDType(assay), space = ScalarSpace())
        key <- AnnDataReductionKey(clean_name)
        r_grp$create_attr(attr_name = "key", robj = key,
                          dtype = CachedGuessDType(key), space = ScalarSpace())
        r_grp$create_attr(attr_name = "global", robj = 0L,
                          dtype = GuessDType(0L), space = ScalarSpace())
        r_grp$create_group("misc")

        node_type <- .zarr_node_type(zarr_path, file.path("obsm", rn))
        if (node_type == "array") {
          emb <- .zarr_read_numeric(zarr_path, file.path("obsm", rn))
          # zarr stores [n_cells, n_comp] in C order.
          # hdf5r writes R matrices in Fortran order, so writing [n_cells, n_comp]
          # from R produces the same HDF5 layout as writeH5Seurat (which writes
          # Embeddings() directly -- also [n_cells, n_comp] in R).
          # No transpose needed.
          if (is.matrix(emb) && nrow(emb) != n_cells && ncol(emb) == n_cells) {
            emb <- t(emb)  # fix only if dimensions are swapped
          }
          r_grp$create_dataset("cell.embeddings", robj = emb,
                               chunk_dims = dim(emb), gzip_level = gzip)
        } else if (node_type == "group") {
          # Sparse embedding -- read and densify
          mat <- .zarr_read_anndata_matrix(zarr_path, file.path("obsm", rn),
                                            transpose = FALSE)
          mat <- as.matrix(mat)
          if (nrow(mat) != n_cells && ncol(mat) == n_cells) mat <- t(mat)
          r_grp$create_dataset("cell.embeddings", robj = mat,
                               chunk_dims = dim(mat), gzip_level = gzip)
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream reduction ", clean_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- graphs (obsp -> graphs) ---
  graphs_grp <- h5$create_group("graphs")
  if (.zarr_node_type(zarr_path, "obsp") == "group") {
    for (gn in .zarr_list_children(zarr_path, "obsp")) {
      if (verbose) message("  Streaming graph: ", gn)
      tryCatch({
        gp <- file.path("obsp", gn)
        attrs <- .zarr_read_attrs(zarr_path, gp)
        enc <- attrs[["encoding-type"]] %||% ""

        if (enc %in% c("csr_matrix", "csc_matrix")) {
          data_vals <- .zarr_read_numeric(zarr_path, file.path(gp, "data"))
          indices <- .zarr_read_numeric(zarr_path, file.path(gp, "indices"))
          indptr <- .zarr_read_numeric(zarr_path, file.path(gp, "indptr"))
          shape <- attrs[["shape"]]
          if (is.list(shape)) shape <- unlist(shape)

          if (enc == "csr_matrix") {
            # CSR -> CSC: need actual conversion via dgCMatrix
            mat <- Matrix::sparseMatrix(
              i = rep(seq_len(shape[1]), diff(as.integer(indptr))),
              j = as.integer(indices) + 1L,
              x = as.numeric(data_vals),
              dims = as.integer(shape)
            )
            mat <- as(mat, "dgCMatrix")
            g_grp <- graphs_grp$create_group(gn)
            g_grp$create_dataset("data", robj = mat@x,
                                 chunk_dims = min(length(mat@x), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indices", robj = mat@i,
                                 chunk_dims = min(length(mat@i), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indptr", robj = mat@p,
                                 chunk_dims = min(length(mat@p), 65536L),
                                 gzip_level = gzip)
            g_grp$create_attr(attr_name = "dims", robj = dim(mat),
                              dtype = GuessDType(dim(mat)))
          } else {
            # CSC -- direct copy (same as h5seurat native format)
            g_grp <- graphs_grp$create_group(gn)
            g_grp$create_dataset("data", robj = as.numeric(data_vals),
                                 chunk_dims = min(length(data_vals), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indices", robj = as.integer(indices),
                                 chunk_dims = min(length(indices), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indptr", robj = as.integer(indptr),
                                 chunk_dims = min(length(indptr), 65536L),
                                 gzip_level = gzip)
            g_grp$create_attr(attr_name = "dims", robj = as.integer(shape),
                              dtype = GuessDType(as.integer(shape)))
          }
          g_grp$create_attr(attr_name = "assay.used", robj = assay,
                            dtype = CachedGuessDType(assay), space = ScalarSpace())
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream graph ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- misc (uns -> misc) ---
  misc_grp <- h5$create_group("misc")
  if (.zarr_node_type(zarr_path, "uns") == "group") {
    for (un in .zarr_list_children(zarr_path, "uns")) {
      tryCatch({
        node <- .zarr_node_type(zarr_path, file.path("uns", un))
        if (node == "array") {
          meta <- .zarr_read_json(file.path(zarr_path, "uns", un, ".zarray"))
          if (!is.null(meta$dtype) && meta$dtype == "|O") {
            vals <- .zarr_read_strings(zarr_path, file.path("uns", un))
            misc_grp$create_dataset(un, robj = vals, dtype = CachedUtf8Type(),
                                    chunk_dims = length(vals), gzip_level = gzip)
          } else {
            vals <- .zarr_read_numeric(zarr_path, file.path("uns", un))
            misc_grp$create_dataset(un, robj = vals,
                                    chunk_dims = length(vals), gzip_level = gzip)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream uns ", un, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- Required empty groups ---
  for (g in c("tools", "commands", "images", "neighbors")) {
    h5$create_group(g)
  }

  # --- active.ident (default: all cells = single identity) ---
  ident_grp <- h5$create_group("active.ident")
  ident_grp$create_dataset("levels", robj = "scConvert",
                           dtype = CachedUtf8Type(),
                           chunk_dims = 1L, gzip_level = gzip)
  ident_grp$create_dataset("values", robj = rep(1L, n_cells),
                           chunk_dims = n_cells, gzip_level = gzip)

  if (verbose) message("Streaming zarr -> h5seurat complete: ", h5seurat_path)
  invisible(h5seurat_path)
}

#' Stream zarr sparse or dense matrix to h5seurat format
#'
#' AnnData stores cells x genes (CSR). h5seurat needs genes x cells (CSC).
#' CSR of (cells x genes) == CSC of (genes x cells) -- zero-copy reinterpretation.
#'
#' @keywords internal
#' @noRd
.stream_zarr_matrix_to_h5seurat <- function(store_path, zarr_mat_path,
                                             h5_parent, h5_name, gzip,
                                             transpose = TRUE) {
  node_type <- .zarr_node_type(store_path, zarr_mat_path)
  if (node_type == "missing") return()

  chunk_size <- 65536L

  if (node_type == "group") {
    attrs <- .zarr_read_attrs(store_path, zarr_mat_path)
    enc <- attrs[["encoding-type"]] %||% ""

    if (enc %in% c("csr_matrix", "csc_matrix")) {
      data_vals <- .zarr_read_numeric(store_path,
                                       file.path(zarr_mat_path, "data"))
      indices <- .zarr_read_numeric(store_path,
                                     file.path(zarr_mat_path, "indices"))
      indptr <- .zarr_read_numeric(store_path,
                                    file.path(zarr_mat_path, "indptr"))
      shape <- attrs[["shape"]]
      if (is.list(shape)) shape <- unlist(shape)

      if (transpose && enc == "csr_matrix") {
        # CSR(cells x genes) == CSC(genes x cells) -- zero-copy
        dims <- as.integer(rev(shape))
      } else if (transpose && enc == "csc_matrix") {
        # CSC already, but we need to transpose -> build dgCMatrix + reinterpret
        mat <- Matrix::sparseMatrix(
          i = as.integer(indices) + 1L, p = as.integer(indptr),
          x = as.numeric(data_vals), dims = as.integer(shape), repr = "C"
        )
        mat <- Matrix::t(mat)
        data_vals <- mat@x; indices <- mat@i; indptr <- mat@p
        dims <- dim(mat)
      } else {
        dims <- as.integer(shape)
      }

      grp <- h5_parent$create_group(h5_name)
      grp$create_dataset("data", robj = as.numeric(data_vals),
                         chunk_dims = min(length(data_vals), chunk_size),
                         gzip_level = gzip)
      grp$create_dataset("indices", robj = as.integer(indices),
                         chunk_dims = min(length(indices), chunk_size),
                         gzip_level = gzip)
      grp$create_dataset("indptr", robj = as.integer(indptr),
                         chunk_dims = min(length(indptr), chunk_size),
                         gzip_level = gzip)
      grp$create_attr(attr_name = "dims", robj = as.integer(dims),
                      dtype = GuessDType(as.integer(dims)))
    }
  } else if (node_type == "array") {
    mat <- .zarr_read_numeric(store_path, zarr_mat_path)
    if (is.matrix(mat) && transpose) mat <- t(mat)
    h5_parent$create_dataset(h5_name, robj = mat,
                             chunk_dims = if (is.matrix(mat)) dim(mat) else length(mat),
                             gzip_level = gzip)
  }
}

