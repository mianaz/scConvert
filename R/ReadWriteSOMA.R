#' @include zzz.R
#' @importFrom Matrix t
#' @importFrom Seurat CreateSeuratObject CreateAssayObject CreateDimReducObject
#'   DefaultAssay
#' @importFrom SeuratObject Cells
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SOMA (TileDB-SOMA) I/O
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read a TileDB-SOMA experiment as a Seurat object
#'
#' Reads a TileDB-SOMA experiment (local or remote) and returns a Seurat object.
#' Supports cell/feature filtering via value queries, multiple measurements
#' (assays), embeddings (obsm), and neighbor graphs (obsp).
#'
#' TileDB-SOMA is the cloud-native format underlying CELLxGENE Census (61M+
#' cells). This function requires the \pkg{tiledbsoma} R package.
#'
#' @param uri URI to a SOMA experiment. Can be a local path or a cloud URI
#'   (e.g., \code{tiledb://...} or \code{s3://...}).
#' @param measurement Name of the measurement to use as the default assay
#'   (default: \code{"RNA"}).
#' @param obs_query Optional value filter string for cells, passed to
#'   \code{SOMADataFrame$read(value_filter = ...)}. For example,
#'   \code{"cell_type == 'T cell'"}.
#' @param var_query Optional value filter string for features.
#' @param obs_column_names Optional character vector of obs column names to
#'   read. If \code{NULL} (default), all columns are read.
#' @param var_column_names Optional character vector of var column names to
#'   read. If \code{NULL} (default), all columns are read.
#' @param verbose Show progress messages
#'
#' @return A \code{Seurat} object
#'
#' @details
#' The reader performs the following steps:
#' \enumerate{
#'   \item Opens the SOMA experiment at \code{uri}
#'   \item Reads obs (cell metadata) with optional filtering
#'   \item Reads var (feature metadata) from the specified measurement
#'   \item Reads the X matrix (\code{data} layer) as a sparse matrix
#'   \item Reads any obsm entries as dimensional reductions
#'   \item Reads any obsp entries as neighbor graphs
#'   \item Reads additional measurements as extra Seurat assays
#' }
#'
#' Cell names are determined from the following columns in priority order:
#' \code{obs_id}, \code{_index}, \code{index}, \code{barcode}, \code{cell_id}.
#' If none are present, synthetic names (\code{cell_1}, \code{cell_2}, ...)
#' are generated.
#'
#' @examples
#' \dontrun{
#' # Read a local SOMA experiment
#' obj <- readSOMA("path/to/experiment")
#'
#' # Read from CELLxGENE Census with filtering
#' obj <- readSOMA(
#'   uri = "s3://cellxgene-census/soma/...",
#'   obs_query = "cell_type == 'T cell' & tissue == 'blood'"
#' )
#'
#' # Read a specific measurement
#' obj <- readSOMA("path/to/experiment", measurement = "ATAC")
#' }
#'
#' @seealso \code{\link{writeSOMA}}, \code{\link{scConvert}}
#'
#' @export
#'
readSOMA <- function(uri,
                     measurement = "RNA",
                     obs_query = NULL,
                     var_query = NULL,
                     obs_column_names = NULL,
                     var_column_names = NULL,
                     verbose = TRUE) {
  if (!requireNamespace("tiledbsoma", quietly = TRUE)) {
    stop(
      "The 'tiledbsoma' package is required for SOMA support.\n",
      "Install with: install.packages('tiledbsoma')",
      call. = FALSE
    )
  }

  # Open experiment -----------------------------------------------------------
  if (verbose) message("Opening SOMA experiment: ", uri)
  exp <- tryCatch(
    tiledbsoma::SOMAExperimentOpen(uri),
    error = function(e) {
      stop("Failed to open SOMA experiment at '", uri, "': ", e$message,
           call. = FALSE)
    }
  )
  on.exit(tryCatch(exp$close(), error = function(e) NULL), add = TRUE)

  # Read obs (cell metadata) --------------------------------------------------
  if (verbose) message("Reading obs (cell metadata)...")
  obs_df <- .soma_read_obs(exp, obs_query = obs_query,
                           column_names = obs_column_names)
  n_cells <- nrow(obs_df)
  if (n_cells == 0L) {
    stop("No cells found in SOMA experiment (after filtering).", call. = FALSE)
  }
  if (verbose) message("  ", n_cells, " cells")

  # Get measurement -----------------------------------------------------------
  ms <- tryCatch(
    exp$ms$get(measurement),
    error = function(e) {
      available <- tryCatch(exp$ms$names(), error = function(e2) character(0))
      stop("Measurement '", measurement, "' not found. Available: ",
           paste(available, collapse = ", "), call. = FALSE)
    }
  )

  # Read var (feature metadata) -----------------------------------------------
  if (verbose) message("Reading var (feature metadata)...")
  var_df <- .soma_read_var(ms, var_query = var_query,
                           column_names = var_column_names)
  n_features <- nrow(var_df)
  if (n_features == 0L) {
    stop("No features found in measurement '", measurement,
         "' (after filtering).", call. = FALSE)
  }
  if (verbose) message("  ", n_features, " features")

  feature_names <- rownames(var_df)

  # Read X (expression data) --------------------------------------------------
  if (verbose) message("Reading X matrix...")
  mat <- .soma_read_x(ms, n_cells = n_cells, n_features = n_features,
                       cell_names = rownames(obs_df),
                       feature_names = feature_names)
  if (verbose) message("  Matrix: ", nrow(mat), " features x ", ncol(mat), " cells")

  # Drop soma_joinid from obs/var before creating Seurat (internal column)
  obs_meta <- obs_df[, !names(obs_df) %in% "soma_joinid", drop = FALSE]
  var_meta <- var_df[, !names(var_df) %in% "soma_joinid", drop = FALSE]

  # Create Seurat object ------------------------------------------------------
  if (verbose) message("Creating Seurat object...")
  obj <- tryCatch(
    Seurat::CreateSeuratObject(
      counts = mat, meta.data = obs_meta,
      assay = measurement, min.cells = 0, min.features = 0
    ),
    error = function(e) {
      stop("Failed to create Seurat object: ", e$message, call. = FALSE)
    }
  )

  # Store var metadata in the assay's meta.features
  if (ncol(var_meta) > 0) {
    tryCatch({
      assay_obj <- obj[[measurement]]
      # Only keep columns not already present
      existing <- names(assay_obj[[]])
      new_cols <- setdiff(names(var_meta), existing)
      if (length(new_cols) > 0) {
        assay_obj[[new_cols]] <- var_meta[feature_names, new_cols, drop = FALSE]
        obj[[measurement]] <- assay_obj
      }
    }, error = function(e) {
      if (verbose) warning("Could not attach var metadata: ", e$message,
                           call. = FALSE)
    })
  }

  # Read obsm (embeddings) ----------------------------------------------------
  obj <- .soma_read_obsm(ms, obj, measurement = measurement,
                          cell_names = rownames(obs_df), verbose = verbose)

  # Read obsp (graphs) --------------------------------------------------------
  obj <- .soma_read_obsp(ms, obj, measurement = measurement,
                          cell_names = rownames(obs_df), verbose = verbose)

  # Read additional measurements as extra assays ------------------------------
  obj <- .soma_read_extra_measurements(exp, obj, primary = measurement,
                                        n_cells = n_cells,
                                        cell_names = rownames(obs_df),
                                        verbose = verbose)

  if (verbose) message("Done.")
  obj
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# writeSOMA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write a Seurat object to TileDB-SOMA format
#'
#' Writes a Seurat object to a TileDB-SOMA experiment. Uses the
#' \code{tiledbsoma::write_soma()} function from the \pkg{tiledbsoma} package,
#' which handles the full conversion of assays, metadata, embeddings, and
#' graphs.
#'
#' @param object A Seurat object
#' @param uri Output URI (local path or cloud URI such as \code{tiledb://...}
#'   or \code{s3://...})
#' @param measurement Name for the primary measurement. If \code{NULL}
#'   (default), uses \code{DefaultAssay(object)}.
#' @param overwrite If \code{TRUE}, remove any existing SOMA at \code{uri}
#'   before writing. Default is \code{FALSE}.
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{uri}
#'
#' @examples
#' \dontrun{
#' library(scConvert)
#' library(Seurat)
#'
#' # Write a Seurat object to SOMA
#' writeSOMA(seurat_obj, uri = "output_soma")
#'
#' # Overwrite existing
#' writeSOMA(seurat_obj, uri = "output_soma", overwrite = TRUE)
#' }
#'
#' @seealso \code{\link{readSOMA}}, \code{\link{scConvert}}
#'
#' @export
#'
writeSOMA <- function(object, uri, measurement = NULL, overwrite = FALSE,
                      verbose = TRUE) {
  if (!requireNamespace("tiledbsoma", quietly = TRUE)) {
    stop(
      "The 'tiledbsoma' package is required for SOMA support.\n",
      "Install with: install.packages('tiledbsoma')",
      call. = FALSE
    )
  }
  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  if (!.soma_has_write_soma()) {
    stop(
      "tiledbsoma::write_soma() is not available. ",
      "Please update the tiledbsoma package to a version that includes ",
      "write_soma() support.",
      call. = FALSE
    )
  }

  if ((dir.exists(uri) || file.exists(uri)) && !overwrite) {
    stop("Output exists: ", uri, ". Use overwrite = TRUE.", call. = FALSE)
  }

  if (is.null(measurement)) {
    measurement <- Seurat::DefaultAssay(object)
  }

  # Atomic write: tiledbsoma::write_soma creates a directory tree. Write it
  # under a sibling temp name and rename on success so a mid-write crash
  # leaves no plausible-but-corrupt store at `uri`.
  parent <- dirname(uri)
  if (!nzchar(parent)) parent <- "."
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  tmp_uri <- tempfile(pattern = "scconvert-soma-", tmpdir = parent)

  if (verbose) message("Writing Seurat object to SOMA: ", uri)

  tryCatch(
    tiledbsoma::write_soma(object, uri = tmp_uri),
    error = function(e) {
      unlink(tmp_uri, recursive = TRUE, force = TRUE)
      stop("Failed to write SOMA: ", e$message, call. = FALSE)
    }
  )

  if (dir.exists(uri) || file.exists(uri)) {
    unlink(uri, recursive = TRUE, force = TRUE)
  }
  if (!file.rename(tmp_uri, uri)) {
    unlink(tmp_uri, recursive = TRUE, force = TRUE)
    stop("Failed to rename temporary SOMA store to final path: ", uri,
         call. = FALSE)
  }

  if (verbose) message("Done: ", uri)
  invisible(uri)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Direct conversion helpers: SOMA <-> other formats
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert h5ad to SOMA format
#'
#' Reads an h5ad file and writes it to TileDB-SOMA format. Since SOMA is not
#' an HDF5 format, this conversion routes through a Seurat intermediate.
#'
#' @param source Path to input .h5ad file
#' @param dest Output URI for the SOMA experiment
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readH5AD}}, \code{\link{writeSOMA}}, \code{\link{SOMAToH5AD}}
#'
#' @export
#'
H5ADToSOMA <- function(source, dest, overwrite = FALSE, verbose = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if ((dir.exists(dest) || file.exists(dest)) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting h5ad -> SOMA")
  obj <- readH5AD(source, verbose = verbose)
  writeSOMA(obj, uri = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#' Convert SOMA to h5ad format
#'
#' Reads a SOMA experiment and writes it to h5ad format. Since SOMA is not
#' an HDF5 format, this conversion routes through a Seurat intermediate.
#'
#' @param source URI of the input SOMA experiment
#' @param dest Path for the output .h5ad file
#' @param measurement Name of the measurement to export (default: \code{"RNA"})
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readSOMA}}, \code{\link{writeH5AD}}, \code{\link{H5ADToSOMA}}
#'
#' @export
#'
SOMAToH5AD <- function(source, dest, measurement = "RNA", overwrite = FALSE,
                       verbose = TRUE) {
  if (file.exists(dest) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting SOMA -> h5ad")
  obj <- readSOMA(source, measurement = measurement, verbose = verbose)
  writeH5AD(obj, filename = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#' Convert h5Seurat to SOMA format
#'
#' Reads an h5Seurat file and writes it to TileDB-SOMA format. Routes through
#' a Seurat intermediate.
#'
#' @param source Path to input .h5seurat file
#' @param dest Output URI for the SOMA experiment
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readH5Seurat}}, \code{\link{writeSOMA}},
#'   \code{\link{SOMAToH5Seurat}}
#'
#' @export
#'
H5SeuratToSOMA <- function(source, dest, overwrite = FALSE, verbose = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if ((dir.exists(dest) || file.exists(dest)) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting h5Seurat -> SOMA")
  obj <- readH5Seurat(file = source, verbose = verbose)
  writeSOMA(obj, uri = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#' Convert SOMA to h5Seurat format
#'
#' Reads a SOMA experiment and writes it to h5Seurat format. Routes through
#' a Seurat intermediate.
#'
#' @param source URI of the input SOMA experiment
#' @param dest Path for the output .h5seurat file
#' @param measurement Name of the measurement to export (default: \code{"RNA"})
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readSOMA}}, \code{\link{writeH5Seurat}},
#'   \code{\link{H5SeuratToSOMA}}
#'
#' @export
#'
SOMAToH5Seurat <- function(source, dest, measurement = "RNA",
                           overwrite = FALSE, verbose = TRUE) {
  if (file.exists(dest) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting SOMA -> h5Seurat")
  obj <- readSOMA(source, measurement = measurement, verbose = verbose)
  writeH5Seurat(object = obj, filename = dest, overwrite = overwrite,
                verbose = verbose)
  invisible(dest)
}


#' Convert H5MU (MuData) to SOMA format
#'
#' Reads an h5mu file and writes it to TileDB-SOMA format. Routes through
#' a Seurat intermediate.
#'
#' @param source Path to input .h5mu file
#' @param dest Output URI for the SOMA experiment
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readH5MU}}, \code{\link{writeSOMA}},
#'   \code{\link{SOMAToH5MU}}
#'
#' @export
#'
H5MUToSOMA <- function(source, dest, overwrite = FALSE, verbose = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if ((dir.exists(dest) || file.exists(dest)) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting h5mu -> SOMA")
  obj <- readH5MU(file = source, verbose = verbose)
  writeSOMA(obj, uri = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#' Convert SOMA to H5MU (MuData) format
#'
#' Reads a SOMA experiment and writes it to h5mu format. Routes through
#' a Seurat intermediate. Best suited for SOMA experiments with multiple
#' measurements (assays).
#'
#' @param source URI of the input SOMA experiment
#' @param dest Path for the output .h5mu file
#' @param measurement Name of the primary measurement (default: \code{"RNA"})
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readSOMA}}, \code{\link{writeH5MU}},
#'   \code{\link{H5MUToSOMA}}
#'
#' @export
#'
SOMAToH5MU <- function(source, dest, measurement = "RNA", overwrite = FALSE,
                       verbose = TRUE) {
  if (file.exists(dest) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting SOMA -> h5mu")
  obj <- readSOMA(source, measurement = measurement, verbose = verbose)
  writeH5MU(object = obj, filename = dest, overwrite = overwrite,
            verbose = verbose)
  invisible(dest)
}


#' Convert Loom to SOMA format
#'
#' Reads a Loom file and writes it to TileDB-SOMA format. Routes through
#' a Seurat intermediate.
#'
#' @param source Path to input .loom file
#' @param dest Output URI for the SOMA experiment
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readLoom}}, \code{\link{writeSOMA}},
#'   \code{\link{SOMAToLoom}}
#'
#' @export
#'
LoomToSOMA <- function(source, dest, overwrite = FALSE, verbose = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if ((dir.exists(dest) || file.exists(dest)) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting Loom -> SOMA")
  obj <- readLoom(file = source, verbose = verbose)
  writeSOMA(obj, uri = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#' Convert SOMA to Loom format
#'
#' Reads a SOMA experiment and writes it to Loom format. Routes through
#' a Seurat intermediate.
#'
#' @param source URI of the input SOMA experiment
#' @param dest Path for the output .loom file
#' @param measurement Name of the measurement to export (default: \code{"RNA"})
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readSOMA}}, \code{\link{writeLoom}},
#'   \code{\link{LoomToSOMA}}
#'
#' @export
#'
SOMAToLoom <- function(source, dest, measurement = "RNA", overwrite = FALSE,
                       verbose = TRUE) {
  if (file.exists(dest) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting SOMA -> Loom")
  obj <- readSOMA(source, measurement = measurement, verbose = verbose)
  writeLoom(object = obj, filename = dest, overwrite = overwrite,
            verbose = verbose)
  invisible(dest)
}


#' Convert Zarr to SOMA format
#'
#' Reads an AnnData Zarr store and writes it to TileDB-SOMA format. Routes
#' through a Seurat intermediate.
#'
#' @param source Path to input .zarr directory
#' @param dest Output URI for the SOMA experiment
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readZarr}}, \code{\link{writeSOMA}},
#'   \code{\link{SOMAToZarr}}
#'
#' @export
#'
ZarrToSOMA <- function(source, dest, overwrite = FALSE, verbose = TRUE) {
  if (!dir.exists(source))
    stop("Input directory not found: ", source, call. = FALSE)
  if ((dir.exists(dest) || file.exists(dest)) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting Zarr -> SOMA")
  obj <- readZarr(source, verbose = verbose)
  writeSOMA(obj, uri = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#' Convert SOMA to Zarr format
#'
#' Reads a SOMA experiment and writes it to AnnData Zarr format. Routes
#' through a Seurat intermediate.
#'
#' @param source URI of the input SOMA experiment
#' @param dest Path for the output .zarr directory
#' @param measurement Name of the measurement to export (default: \code{"RNA"})
#' @param overwrite If \code{TRUE}, overwrite existing output
#' @param verbose Show progress messages
#'
#' @return Invisibly returns \code{dest}
#'
#' @seealso \code{\link{readSOMA}}, \code{\link{writeZarr}},
#'   \code{\link{ZarrToSOMA}}
#'
#' @export
#'
SOMAToZarr <- function(source, dest, measurement = "RNA", overwrite = FALSE,
                       verbose = TRUE) {
  if ((dir.exists(dest) || file.exists(dest)) && !overwrite)
    stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)

  if (verbose) message("Converting SOMA -> Zarr")
  obj <- readSOMA(source, measurement = measurement, verbose = verbose)
  writeZarr(obj, filename = dest, overwrite = overwrite, verbose = verbose)
  invisible(dest)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers for readSOMA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read obs from a SOMA experiment
#'
#' @param exp A SOMAExperiment object
#' @param obs_query Optional value filter string
#' @param column_names Optional column names to read
#'
#' @return A data.frame with cell metadata, rownames set to cell identifiers
#'
#' @keywords internal
#'
.soma_read_obs <- function(exp, obs_query = NULL, column_names = NULL) {
  read_args <- list()
  if (!is.null(obs_query)) {
    read_args$value_filter <- obs_query
  }
  if (!is.null(column_names)) {
    # Always include soma_joinid and potential index columns
    idx_candidates <- c("soma_joinid", "obs_id", "_index", "index",
                        "barcode", "cell_id")
    column_names <- unique(c(intersect(idx_candidates, column_names),
                             column_names))
    read_args$column_names <- column_names
  }

  obs_iter <- tryCatch(
    do.call(exp$obs$read, read_args),
    error = function(e) {
      stop("Failed to read obs: ", e$message, call. = FALSE)
    }
  )

  obs_tbl <- tryCatch(
    obs_iter$concat(),
    error = function(e) {
      stop("Failed to concatenate obs iterator: ", e$message, call. = FALSE)
    }
  )

  obs_df <- as.data.frame(obs_tbl)
  if (nrow(obs_df) == 0L) return(obs_df)

  # Set rownames from the best available index column
  idx_col <- intersect(
    c("obs_id", "_index", "index", "barcode", "cell_id"),
    names(obs_df)
  )
  if (length(idx_col) > 0) {
    row_ids <- as.character(obs_df[[idx_col[1]]])
    # Handle duplicates
    if (anyDuplicated(row_ids)) {
      row_ids <- make.unique(row_ids)
    }
    rownames(obs_df) <- row_ids
  } else if ("soma_joinid" %in% names(obs_df)) {
    rownames(obs_df) <- paste0("cell_", obs_df$soma_joinid)
  } else {
    rownames(obs_df) <- paste0("cell_", seq_len(nrow(obs_df)))
  }

  obs_df
}


#' Read var from a SOMA measurement
#'
#' @param ms A SOMAMeasurement object
#' @param var_query Optional value filter string
#' @param column_names Optional column names to read
#'
#' @return A data.frame with feature metadata, rownames set to feature names
#'
#' @keywords internal
#'
.soma_read_var <- function(ms, var_query = NULL, column_names = NULL) {
  read_args <- list()
  if (!is.null(var_query)) {
    read_args$value_filter <- var_query
  }
  if (!is.null(column_names)) {
    idx_candidates <- c("soma_joinid", "var_id", "_index", "index",
                        "gene_id", "feature_id", "feature_name")
    column_names <- unique(c(intersect(idx_candidates, column_names),
                             column_names))
    read_args$column_names <- column_names
  }

  var_iter <- tryCatch(
    do.call(ms$var$read, read_args),
    error = function(e) {
      stop("Failed to read var: ", e$message, call. = FALSE)
    }
  )

  var_tbl <- tryCatch(
    var_iter$concat(),
    error = function(e) {
      stop("Failed to concatenate var iterator: ", e$message, call. = FALSE)
    }
  )

  var_df <- as.data.frame(var_tbl)
  if (nrow(var_df) == 0L) return(var_df)

  # Set rownames from the best available index column
  idx_col <- intersect(
    c("var_id", "_index", "index", "gene_id", "feature_id", "feature_name"),
    names(var_df)
  )
  if (length(idx_col) > 0) {
    feat_ids <- as.character(var_df[[idx_col[1]]])
    if (anyDuplicated(feat_ids)) {
      feat_ids <- make.unique(feat_ids)
    }
    rownames(var_df) <- feat_ids
  } else if ("soma_joinid" %in% names(var_df)) {
    rownames(var_df) <- paste0("gene_", var_df$soma_joinid)
  } else {
    rownames(var_df) <- paste0("gene_", seq_len(nrow(var_df)))
  }

  var_df
}


#' Read X matrix from a SOMA measurement
#'
#' @param ms A SOMAMeasurement object
#' @param n_cells Expected number of cells
#' @param n_features Expected number of features
#' @param cell_names Character vector of cell names
#' @param feature_names Character vector of feature names
#'
#' @return A sparse matrix (features x cells)
#'
#' @keywords internal
#'
.soma_read_x <- function(ms, n_cells, n_features, cell_names, feature_names) {
  # Try "data" layer first, then "counts", then "raw"
  x_layer <- NULL
  for (layer_name in c("data", "counts", "raw")) {
    x_layer <- tryCatch(ms$X$get(layer_name), error = function(e) NULL)
    if (!is.null(x_layer)) break
  }
  if (is.null(x_layer)) {
    # Try to get whatever is available
    available <- tryCatch(ms$X$names(), error = function(e) character(0))
    if (length(available) > 0) {
      x_layer <- ms$X$get(available[1])
    } else {
      stop("No X layers found in measurement. ",
           "Expected at least one of: data, counts, raw", call. = FALSE)
    }
  }

  mat <- tryCatch({
    iter <- x_layer$read()
    sp <- iter$sparse_matrix(zero_based = TRUE)
    result <- sp$concat()
    # tiledbsoma 2.x returns matrixZeroBasedView R6; convert to standard sparse
    if (inherits(result, "R6")) {
      result$get_one_based_matrix()
    } else {
      result
    }
  }, error = function(e) {
    stop("Failed to read X matrix: ", e$message, call. = FALSE)
  })

  # Ensure features x cells orientation (SOMA stores cells x features)
  mat_rows <- nrow(mat)
  mat_cols <- ncol(mat)

  if (mat_rows == n_cells && mat_cols == n_features) {
    # cells x features -> transpose to features x cells
    mat <- Matrix::t(mat)
  } else if (mat_rows == n_features && mat_cols == n_cells) {
    # Already features x cells
  } else {
    # Dimensions don't match exactly; try heuristic
    if (mat_rows > mat_cols && n_features > n_cells) {
      # Likely features x cells already
    } else if (mat_rows < mat_cols && n_features < n_cells) {
      # Likely features x cells already
    } else {
      # Ambiguous; assume cells x features (SOMA convention) and transpose
      mat <- Matrix::t(mat)
    }
  }

  # Coerce to dgCMatrix for Seurat compatibility
  if (!inherits(mat, "dgCMatrix")) {
    mat <- tryCatch(
      methods::as(mat, "dgCMatrix"),
      error = function(e) {
        tryCatch(
          methods::as(mat, "CsparseMatrix"),
          error = function(e2) mat
        )
      }
    )
  }

  # Set dimnames
  if (nrow(mat) == length(feature_names)) {
    rownames(mat) <- feature_names
  }
  if (ncol(mat) == length(cell_names)) {
    colnames(mat) <- cell_names
  }

  mat
}


#' Read obsm embeddings from a SOMA measurement
#'
#' @param ms A SOMAMeasurement object
#' @param obj A Seurat object to add reductions to
#' @param measurement Assay name
#' @param cell_names Cell names for rownames
#' @param verbose Verbosity
#'
#' @return The Seurat object with embeddings added
#'
#' @keywords internal
#'
.soma_read_obsm <- function(ms, obj, measurement, cell_names, verbose = TRUE) {
  obsm_coll <- tryCatch(ms$obsm, error = function(e) NULL)
  if (is.null(obsm_coll)) return(obj)

  obsm_names <- tryCatch(obsm_coll$names(), error = function(e) character(0))
  if (length(obsm_names) == 0L) return(obj)

  for (emb_name in obsm_names) {
    tryCatch({
      if (verbose) message("  Reading obsm: ", emb_name)
      emb_arr <- obsm_coll$get(emb_name)
      rd <- emb_arr$read()
      # tiledbsoma 2.x: use sparse_matrix path for SOMASparseNDArray
      emb_result <- if ("sparse_matrix" %in% ls(rd)) {
        m <- rd$sparse_matrix(zero_based = TRUE)$concat()
        if (inherits(m, "R6")) m$get_one_based_matrix() else m
      } else {
        rd$concat()
      }
      emb_mat <- as.matrix(emb_result)

      # Remove soma_joinid column if present
      if ("soma_joinid" %in% colnames(emb_mat)) {
        emb_mat <- emb_mat[, colnames(emb_mat) != "soma_joinid", drop = FALSE]
      }

      if (nrow(emb_mat) != length(cell_names)) {
        if (verbose) warning("Skipping obsm '", emb_name,
                             "': row count mismatch (", nrow(emb_mat),
                             " vs ", length(cell_names), " cells)",
                             call. = FALSE)
        next
      }
      rownames(emb_mat) <- cell_names

      # Clean up reduction name (X_pca -> pca, X_umap -> umap)
      clean_name <- gsub("^X_", "", emb_name)
      key <- paste0(tolower(clean_name), "_")

      # Set column names if missing
      if (is.null(colnames(emb_mat))) {
        colnames(emb_mat) <- paste0(gsub("_$", "", key), "_",
                                    seq_len(ncol(emb_mat)))
      }

      obj[[clean_name]] <- Seurat::CreateDimReducObject(
        embeddings = emb_mat, key = key, assay = measurement
      )
    }, error = function(e) {
      if (verbose) warning("Could not read obsm '", emb_name, "': ",
                           e$message, call. = FALSE)
    })
  }

  obj
}


#' Read obsp graphs from a SOMA measurement
#'
#' @param ms A SOMAMeasurement object
#' @param obj A Seurat object to add graphs to
#' @param measurement Assay name
#' @param cell_names Cell names for row/column names
#' @param verbose Verbosity
#'
#' @return The Seurat object with graphs added
#'
#' @keywords internal
#'
.soma_read_obsp <- function(ms, obj, measurement, cell_names, verbose = TRUE) {
  obsp_coll <- tryCatch(ms$obsp, error = function(e) NULL)
  if (is.null(obsp_coll)) return(obj)

  obsp_names <- tryCatch(obsp_coll$names(), error = function(e) character(0))
  if (length(obsp_names) == 0L) return(obj)

  for (graph_name in obsp_names) {
    tryCatch({
      if (verbose) message("  Reading obsp: ", graph_name)
      graph_arr <- obsp_coll$get(graph_name)
      graph_result <- graph_arr$read()$sparse_matrix(zero_based = TRUE)$concat()
      # tiledbsoma 2.x returns matrixZeroBasedView R6
      graph_data <- if (inherits(graph_result, "R6")) {
        graph_result$get_one_based_matrix()
      } else {
        graph_result
      }

      # Coerce to dgCMatrix
      if (!inherits(graph_data, "dgCMatrix")) {
        graph_data <- tryCatch(
          methods::as(graph_data, "dgCMatrix"),
          error = function(e) graph_data
        )
      }

      if (nrow(graph_data) != length(cell_names)) {
        if (verbose) warning("Skipping obsp '", graph_name,
                             "': dimension mismatch", call. = FALSE)
        next
      }
      rownames(graph_data) <- colnames(graph_data) <- cell_names

      graph_obj <- Seurat::as.Graph(graph_data)
      Seurat::DefaultAssay(graph_obj) <- measurement
      full_name <- paste0(measurement, "_", graph_name)
      obj[[full_name]] <- graph_obj
    }, error = function(e) {
      if (verbose) warning("Could not read obsp '", graph_name, "': ",
                           e$message, call. = FALSE)
    })
  }

  obj
}


#' Read additional measurements as extra assays
#'
#' @param exp A SOMAExperiment object
#' @param obj A Seurat object to add assays to
#' @param primary Name of the primary measurement (already loaded)
#' @param n_cells Number of cells in the primary measurement
#' @param cell_names Cell names
#' @param verbose Verbosity
#'
#' @return The Seurat object with extra assays added
#'
#' @keywords internal
#'
.soma_read_extra_measurements <- function(exp, obj, primary, n_cells,
                                           cell_names, verbose = TRUE) {
  all_ms <- tryCatch(exp$ms$names(), error = function(e) character(0))
  extra_ms <- setdiff(all_ms, primary)
  if (length(extra_ms) == 0L) return(obj)

  for (ms_name in extra_ms) {
    tryCatch({
      if (verbose) message("  Reading measurement: ", ms_name)
      other_ms <- exp$ms$get(ms_name)

      # Read var
      other_var <- as.data.frame(other_ms$var$read()$concat())
      idx <- intersect(c("var_id", "_index", "index", "gene_id",
                         "feature_id", "feature_name"), names(other_var))
      feat <- if (length(idx) > 0) {
        as.character(other_var[[idx[1]]])
      } else {
        paste0("feat_", seq_len(nrow(other_var)))
      }

      # Read X
      x_layer <- NULL
      for (ln in c("data", "counts", "raw")) {
        x_layer <- tryCatch(other_ms$X$get(ln), error = function(e) NULL)
        if (!is.null(x_layer)) break
      }
      if (is.null(x_layer)) {
        available <- tryCatch(other_ms$X$names(), error = function(e) character(0))
        if (length(available) > 0) x_layer <- other_ms$X$get(available[1])
      }
      if (is.null(x_layer)) {
        if (verbose) warning("No X layer in measurement '", ms_name, "'",
                             call. = FALSE)
        next
      }

      other_mat <- x_layer$read()$sparse_matrix(zero_based = TRUE)$concat()

      # Orient to features x cells
      if (nrow(other_mat) == n_cells && ncol(other_mat) == length(feat)) {
        other_mat <- Matrix::t(other_mat)
      }

      # Coerce
      if (!inherits(other_mat, "dgCMatrix")) {
        other_mat <- tryCatch(
          methods::as(other_mat, "dgCMatrix"),
          error = function(e) other_mat
        )
      }

      if (nrow(other_mat) == length(feat)) rownames(other_mat) <- feat
      if (ncol(other_mat) == length(cell_names)) {
        colnames(other_mat) <- cell_names
      }

      obj[[ms_name]] <- Seurat::CreateAssayObject(counts = other_mat)
    }, error = function(e) {
      if (verbose) warning("Could not read measurement '", ms_name, "': ",
                           e$message, call. = FALSE)
    })
  }

  obj
}


#' Check whether tiledbsoma::write_soma is available
#'
#' @return Logical
#'
#' @keywords internal
#'
.soma_has_write_soma <- function() {
  tryCatch(
    is.function(get("write_soma", envir = asNamespace("tiledbsoma"),
                    mode = "function")),
    error = function(e) FALSE
  )
}
