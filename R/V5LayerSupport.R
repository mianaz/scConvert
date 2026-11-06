#' @include zzz.R
#' @importFrom hdf5r H5Group H5File h5attr
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#'
NULL

#' V5 Layer Support Functions
#'
#' Internal functions for reading and writing V5 layer structures in h5Seurat files
#'
#' @name V5LayerSupport
#' @rdname V5LayerSupport
#'
#' @keywords internal
#'
NULL

#' Write V5 Layer Data
#'
#' Write matrix data in V5 layer format (sparse or dense)
#'
#' @param h5_group The HDF5 group to write to
#' @param layer_name Name of the layer (e.g., "counts", "data", "scale.data")
#' @param matrix_data The matrix data to write
#' @param features Character vector of feature names
#' @param verbose Show progress messages
#'
#' @return Invisible NULL
#'
#' @keywords internal
#'
WriteV5Layer <- function(h5_group, layer_name, matrix_data, features = NULL, verbose = FALSE) {
  if (verbose) message(paste0("Writing V5 layer '", layer_name, "'"))

  if (!h5_group$exists("layers")) {
    h5_group$create_group("layers")
  }
  layers_group <- h5_group[["layers"]]

  if (layers_group$exists(layer_name)) {
    if (verbose) message(paste0("Removing existing layer '", layer_name, "'"))
    layers_group$link_delete(layer_name)
  }

  is_sparse <- inherits(matrix_data, "sparseMatrix") || inherits(matrix_data, "dgCMatrix")

  if (is_sparse) {
    if (!inherits(matrix_data, "dgCMatrix")) {
      matrix_data <- as(matrix_data, "dgCMatrix")
    }

    layer_group <- layers_group$create_group(layer_name)

    # Write sparse matrix in CSC format (V5 standard) with chunking + compression
    gzip <- GetCompressionLevel()
    chunk_cap <- 65536L
    for (.ds_info in list(
      list(name = "data",    vals = matrix_data@x),
      list(name = "indices", vals = as.integer(matrix_data@i)),
      list(name = "indptr",  vals = as.integer(matrix_data@p))
    )) {
      .n <- length(.ds_info$vals)
      if (gzip > 0L && .n > 0L) {
        layer_group$create_dataset(
          name = .ds_info$name, robj = .ds_info$vals,
          dtype = GuessDType(x = .ds_info$vals),
          gzip_level = gzip, chunk_dims = min(.n, chunk_cap)
        )
      } else {
        layer_group$create_dataset(
          name = .ds_info$name, robj = .ds_info$vals,
          dtype = GuessDType(x = .ds_info$vals)
        )
      }
    }

    layer_group$create_attr("shape", c(nrow(matrix_data), ncol(matrix_data)))
    layer_group$create_attr("encoding-type", "csc_matrix")
    layer_group$create_attr("encoding-version", "0.1.0")

    if (verbose) message(paste0("  Wrote sparse matrix: ", length(matrix_data@x), " non-zero values, ",
            "shape ", nrow(matrix_data), " x ", ncol(matrix_data)))
  } else {
    layers_group[[layer_name]] <- as.matrix(matrix_data)

    if (verbose) message(paste0("  Wrote dense matrix: shape ", nrow(matrix_data), " x ", ncol(matrix_data)))
  }

  if (!is.null(features) && !h5_group$exists("meta.data")) {
    h5_group$create_group("meta.data")
  }
  if (!is.null(features) && !h5_group[["meta.data"]]$exists("_index")) {
    h5_group[["meta.data"]][["_index"]] <- features
    if (verbose) message(paste0("  Stored ", length(features), " feature names"))
  }

  return(invisible(NULL))
}

#' Read V5 Layer Data
#'
#' Read matrix data from V5 layer format
#'
#' @param h5_group The HDF5 assay group
#' @param layer_name Name of the layer to read
#' @param features Expected feature names (for dimension validation)
#' @param cells Expected cell names (for dimension validation)
#' @param verbose Show progress messages
#'
#' @return A matrix object, or NULL if layer doesn't exist
#'
#' @keywords internal
#'
ReadV5Layer <- function(h5_group, layer_name, features = NULL, cells = NULL, verbose = FALSE) {
  if (!h5_group$exists("layers")) {
    if (verbose) message("No layers group found")
    return(NULL)
  }

  layers_group <- h5_group[["layers"]]

  if (!layers_group$exists(layer_name)) {
    if (verbose) message(paste0("Layer '", layer_name, "' not found"))
    return(NULL)
  }

  layer_obj <- layers_group[[layer_name]]

  if (inherits(layer_obj, "H5Group")) {
    if (verbose) message(paste0("Reading sparse V5 layer '", layer_name, "'"))

    if (exists("ReadSparseMatrix", envir = asNamespace("scConvert"))) {
      mat <- ReadSparseMatrix(layer_obj, verbose = verbose)
    } else {
      if (layer_obj$exists("data") && layer_obj$exists("indices") && layer_obj$exists("indptr")) {
        data_vals <- layer_obj[["data"]][]
        indices <- layer_obj[["indices"]][] + 1L
        indptr <- layer_obj[["indptr"]][]

        if (layer_obj$attr_exists("shape")) {
          shape <- h5attr(x = layer_obj, which = "shape")
          nrows <- shape[1]
          ncols <- shape[2]
        } else {
          ncols <- length(indptr) - 1L
          nrows <- max(indices)
        }

        sparse_mat <- Matrix::sparseMatrix(
          i = indices,
          p = indptr,
          x = data_vals,
          dims = c(nrows, ncols),
          index1 = TRUE
        )

        mat <- as.matrix(sparse_mat)
      } else {
        warning("Unknown sparse matrix format in layer '", layer_name, "'")
        return(NULL)
      }
    }
  } else {
    if (verbose) message(paste0("Reading dense V5 layer '", layer_name, "'"))
    mat <- as.matrix(layer_obj[,])
  }

  if (!is.null(features) && !is.null(cells)) {
    expected_nrows <- length(features)
    expected_ncols <- length(cells)

    if (nrow(mat) != expected_nrows || ncol(mat) != expected_ncols) {
      if (nrow(mat) == expected_ncols && ncol(mat) == expected_nrows) {
        if (verbose) message("Transposing matrix to match expected dimensions")
        mat <- t(mat)
      } else if (nrow(mat) < expected_nrows && ncol(mat) == expected_ncols) {
        if (verbose) message("Expanding matrix to full feature space")
        full_mat <- matrix(0, nrow = expected_nrows, ncol = expected_ncols)
        full_mat[1:nrow(mat), ] <- mat
        mat <- full_mat
      }
    }

    rownames(mat) <- features
    colnames(mat) <- cells
  }

  return(mat)
}


