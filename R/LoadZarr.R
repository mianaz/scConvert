#' @include zzz.R
#' @include AnnDataEncoding.R
#' @importFrom Matrix sparseMatrix t
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom SeuratObject Cells
NULL

#' Apply a user-supplied selection filter to a vector of available names
#'
#' \itemize{
#'   \item \code{NULL} returns \code{available} unchanged (keep all).
#'   \item \code{character(0)} returns \code{character(0)} (drop all).
#'   \item Otherwise returns the intersection.
#' }
#'
#' @keywords internal
.zarr_select_keep <- function(filter, available) {
  if (is.null(filter)) return(available)
  if (!is.character(filter)) {
    stop("filter must be NULL or a character vector", call. = FALSE)
  }
  if (length(filter) == 0L) return(character(0))
  intersect(filter, available)
}

#' Load an AnnData Zarr store as a Seurat object
#'
#' Direct conversion from AnnData Zarr format (.zarr directory) to Seurat object.
#' Supports zarr v2 stores as produced by Python anndata. This enables reading
#' of cloud-native AnnData stores from CELLxGENE, Human Cell Atlas, and
#' SpatialData workflows.
#'
#' \code{file} may also be a remote URL: \code{s3://bucket/key.zarr},
#' \code{gs://bucket/key.zarr}. Anonymous (public) buckets only; private
#' buckets requiring SigV4 signing are not supported. When \code{cache} is
#' \code{TRUE}, the downloaded store is kept under
#' \code{tools::R_user_dir("scConvert", "cache")} and reused on subsequent
#' calls with the same URL.
#'
#' @param file Path to a local .zarr directory, or a remote
#'   \code{s3://}/\code{gs://} URL.
#' @param assay.name Name for the primary assay (default: "RNA")
#' @param verbose Show progress messages
#' @param cache For remote URLs: \code{TRUE} to persist the download under
#'   \code{tools::R_user_dir("scConvert", "cache")}; \code{FALSE} to use a
#'   tempdir that is discarded with the R session. Ignored for local paths.
#' @param layers Character vector of \code{/layers/*} groups to read,
#'   or \code{NULL} (default) for all, or \code{character(0)} to skip all
#'   layers entirely. Skipping unused layers avoids their chunk fetches
#'   on remote stores.
#' @param obs_idx,var_idx Integer indices (1-based) selecting cells / features
#'   to keep. \code{NULL} (default) reads all. Index slicing pushes down to
#'   chunk fetches: only chunks of the cell / feature index, obs columns,
#'   obsm embeddings, var columns, and sparse \code{indptr}-resolved data
#'   blocks containing the requested indices are pulled.
#' @param obsm,obsp,varm,varp,uns Same semantics as \code{layers} for the
#'   corresponding AnnData groups: \code{NULL} keeps everything,
#'   \code{character()} drops the entire group, otherwise reads only the
#'   named items.
#' @param include_x If \code{FALSE}, skip reading the main expression
#'   matrix entirely. Useful for metadata-only reads from a remote store.
#'   Default \code{TRUE}. (If \code{FALSE}, a zero-row placeholder assay
#'   is constructed so the returned Seurat object still has the right
#'   number of cells.)
#' @param ... Additional arguments (currently unused)
#'
#' @return A \code{Seurat} object
#'
#' @export
#'
readZarr <- function(file, assay.name = "RNA", verbose = TRUE,
                     cache = TRUE,
                     obs_idx = NULL, var_idx = NULL,
                     layers = NULL, obsm = NULL, obsp = NULL,
                     varm = NULL, varp = NULL, uns = NULL,
                     include_x = TRUE, ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required for readZarr. ",
         "Install with: install.packages('jsonlite')", call. = FALSE)
  }

  is_remote <- .is_remote_zarr_url(file)
  if (is_remote) {
    if (verbose) message("Opening remote Zarr store (lazy fetch): ", file)
    store <- .zarr_make_store(file, cache = cache)
    # Subsequent helpers receive the store object and reuse the same manifest
    # + chunk cache; no download-then-read mirror.
    file <- store
  } else {
    if (!dir.exists(file)) {
      stop("Zarr store not found: ", file, call. = FALSE)
    }
  }

  # Validate zarr store
  version <- .zarr_store_version(file)
  if (verbose && !is_remote) {
    message("Loading Zarr store (v", version, "): ", file)
  }

  # Read root attributes to verify it's an AnnData store
  root_attrs <- .zarr_read_attrs(file)

  # Resolve obs_idx / var_idx into integer vectors (1-based). These are
  # threaded through every downstream chunk-aware read so that the lazy
  # HTTP store only fetches blocks containing the requested indices.
  n_obs_full <- .zarr_get_n_obs(file)
  n_var_full <- .zarr_get_n_vars(file)
  if (!is.null(obs_idx)) {
    if (is.logical(obs_idx)) obs_idx <- which(obs_idx)
    obs_idx <- as.integer(obs_idx)
    if (any(obs_idx < 1L | obs_idx > n_obs_full)) {
      stop("obs_idx out of range (1..", n_obs_full, ")", call. = FALSE)
    }
  }
  if (!is.null(var_idx)) {
    if (is.logical(var_idx)) var_idx <- which(var_idx)
    var_idx <- as.integer(var_idx)
    if (any(var_idx < 1L | var_idx > n_var_full)) {
      stop("var_idx out of range (1..", n_var_full, ")", call. = FALSE)
    }
  }

  # 1. Read cell names (sliced)
  if (verbose) message("Reading cell names...")
  cell.names <- .zarr_read_anndata_index(file, "obs", slice_idx = obs_idx)

  # 2. Read feature names (sliced)
  if (verbose) message("Reading feature names...")
  feature.names <- .zarr_read_anndata_index(file, "var", slice_idx = var_idx)

  # Generate fallback names if needed (using sliced length, not full).
  n_obs_kept <- if (is.null(obs_idx)) n_obs_full else length(obs_idx)
  n_var_kept <- if (is.null(var_idx)) n_var_full else length(var_idx)
  if (is.null(cell.names)) cell.names <- paste0("Cell", seq_len(n_obs_kept))
  if (is.null(feature.names)) feature.names <- paste0("Gene", seq_len(n_var_kept))

  # 3. Read expression matrix (or skip if include_x = FALSE)
  x_path <- "X"
  if (include_x) {
    if (verbose) message("Reading expression matrix...")
    if (.zarr_node_type(file, x_path) == "missing") {
      stop("No expression matrix (X) found in zarr store", call. = FALSE)
    }
    expr_matrix <- .zarr_read_anndata_matrix(file, x_path,
                                              obs_idx = obs_idx,
                                              var_idx = var_idx)

    # Handle dimension mismatches
    if (nrow(expr_matrix) != length(feature.names)) {
      warning("Adjusting feature names to match matrix dimensions", immediate. = TRUE)
      feature.names <- feature.names[seq_len(nrow(expr_matrix))]
    }
    if (ncol(expr_matrix) != length(cell.names)) {
      warning("Adjusting cell names to match matrix dimensions", immediate. = TRUE)
      cell.names <- cell.names[seq_len(ncol(expr_matrix))]
    }
    rownames(expr_matrix) <- feature.names
    colnames(expr_matrix) <- cell.names
  } else {
    if (verbose) message("Skipping expression matrix (include_x = FALSE)")
    expr_matrix <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(length(feature.names), length(cell.names)),
      dimnames = list(feature.names, cell.names)
    )
  }

  # 4. Create Seurat object
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    project = "Zarr",
    assay = assay.name,
    min.cells = 0,
    min.features = 0
  )

  # 5. Add layers if present
  layers_path <- "layers"
  layer_names <- .zarr_select_keep(layers,
    if (.zarr_node_type(file, layers_path) == "group") {
      .zarr_list_children(file, layers_path)
    } else character(0))
  if (length(layer_names) > 0L) {
    if (verbose) message("Adding layers...")
    for (layer_name in layer_names) {
      if (verbose) message("  Adding layer: ", layer_name)
      tryCatch({
        layer_matrix <- .zarr_read_anndata_matrix(
          file, file.path(layers_path, layer_name),
          obs_idx = obs_idx, var_idx = var_idx
        )
        if (nrow(layer_matrix) == nrow(expr_matrix) &&
            ncol(layer_matrix) == ncol(expr_matrix)) {
          rownames(layer_matrix) <- feature.names
          colnames(layer_matrix) <- cell.names
          seurat_slot <- AnnDataLayerToSeurat(layer_name)
          tryCatch({
            seurat_obj[[assay.name]] <- SetAssayData(
              object = seurat_obj[[assay.name]],
              layer = seurat_slot,
              new.data = layer_matrix
            )
          }, error = function(e) {
            if (verbose) warning("Could not add layer ", layer_name, ": ",
                                 e$message, immediate. = TRUE)
          })
        }
      }, error = function(e) {
        if (verbose) warning("Could not read layer ", layer_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 6. Add cell metadata
  obs_path <- "obs"
  if (.zarr_node_type(file, obs_path) == "group") {
    if (verbose) message("Adding cell metadata...")
    obs_attrs <- .zarr_read_attrs(file, obs_path)
    index_col <- obs_attrs[["_index"]] %||% "_index"
    column_order <- obs_attrs[["column-order"]]

    obs_cols <- if (!is.null(column_order)) {
      column_order
    } else {
      children <- .zarr_list_children(file, obs_path)
      setdiff(children, c("_index", "index", "__categories"))
    }

    for (col in obs_cols) {
      col_path <- file.path(obs_path, col)
      if (.zarr_node_type(file, col_path) == "missing") next
      if (verbose) message("  Adding metadata: ", col)
      tryCatch({
        meta_values <- .zarr_read_anndata_column(file, col_path,
                                                  slice_idx = obs_idx)
        if (!is.null(meta_values) && length(meta_values) == ncol(seurat_obj)) {
          seurat_obj[[col]] <- meta_values
        }
      }, error = function(e) {
        if (verbose) warning("Could not add metadata '", col, "': ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 7. Add dimensional reductions
  obsm_path <- "obsm"
  reduc_names <- .zarr_select_keep(obsm,
    if (.zarr_node_type(file, obsm_path) == "group") {
      .zarr_list_children(file, obsm_path)
    } else character(0))
  if (length(reduc_names) > 0L) {
    if (verbose) message("Adding dimensional reductions...")
    for (reduc_name in reduc_names) {
      clean_name <- gsub("^X_", "", reduc_name)
      if (verbose) message("  Adding reduction: ", clean_name)

      tryCatch({
        # obsm is (n_obs, n_components) — slice rows with obs_idx only.
        obsm_slice <- if (is.null(obs_idx)) NULL else list(d1 = obs_idx)
        embeddings <- .zarr_read_numeric(file, file.path(obsm_path, reduc_name),
                                          slice = obsm_slice)

        # Handle row-major vs column-major dimension ordering
        if (nrow(embeddings) != length(cell.names) &&
            ncol(embeddings) == length(cell.names)) {
          embeddings <- t(embeddings)
        }

        if (nrow(embeddings) != length(cell.names)) {
          warning("Skipping reduction ", clean_name, " - dimension mismatch",
                  immediate. = TRUE)
          next
        }

        rownames(embeddings) <- cell.names
        key <- AnnDataReductionKey(clean_name)
        colnames(embeddings) <- paste0(
          gsub("_$", "", key), "_", seq_len(ncol(embeddings))
        )

        reduc_obj <- CreateDimReducObject(
          embeddings = embeddings,
          key = key,
          assay = assay.name
        )
        seurat_obj[[clean_name]] <- reduc_obj
      }, error = function(e) {
        if (verbose) warning("Could not add reduction ", clean_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 8. Add feature metadata
  var_path <- "var"
  if (.zarr_node_type(file, var_path) == "group") {
    if (verbose) message("Adding feature metadata...")
    var_attrs <- .zarr_read_attrs(file, var_path)
    index_col <- var_attrs[["_index"]] %||% "_index"
    column_order <- var_attrs[["column-order"]]

    var_cols <- if (!is.null(column_order)) {
      column_order
    } else {
      children <- .zarr_list_children(file, var_path)
      setdiff(children, c("_index", "index", "__categories"))
    }

    for (col in var_cols) {
      col_path <- file.path(var_path, col)
      if (.zarr_node_type(file, col_path) == "missing") next
      if (verbose) message("  Adding feature metadata: ", col)

      tryCatch({
        meta_values <- .zarr_read_anndata_column(file, col_path,
                                                  slice_idx = var_idx)
        if (is.null(meta_values)) next

        # Special handling for highly_variable
        if (col == "highly_variable") {
          if (is.factor(meta_values)) {
            meta_values <- as.character(meta_values) == "True"
          } else if (is.numeric(meta_values)) {
            meta_values <- as.logical(meta_values)
          }
        }

        if (length(meta_values) == nrow(seurat_obj)) {
          names(meta_values) <- feature.names
          seurat_obj[[assay.name]][[col]] <- meta_values

          if (col == "highly_variable" && is.logical(meta_values)) {
            VariableFeatures(seurat_obj) <- feature.names[meta_values]
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not add feature metadata ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 9. Add neighbor graphs from obsp
  obsp_path <- "obsp"
  graph_names <- .zarr_select_keep(obsp,
    if (.zarr_node_type(file, obsp_path) == "group") {
      .zarr_list_children(file, obsp_path)
    } else character(0))
  if (length(graph_names) > 0L) {
    if (verbose) message("Adding neighbor graphs...")
    for (graph_name in graph_names) {
      if (verbose) message("  Adding graph: ", graph_name)
      tryCatch({
        graph_matrix <- .zarr_read_anndata_matrix(
          file, file.path(obsp_path, graph_name), transpose = FALSE,
          obs_idx = obs_idx, var_idx = obs_idx
        )
        if (nrow(graph_matrix) == ncol(graph_matrix) &&
            nrow(graph_matrix) == ncol(seurat_obj)) {
          rownames(graph_matrix) <- cell.names
          colnames(graph_matrix) <- cell.names
          seurat_graph_name <- switch(graph_name,
            "connectivities" = paste0(assay.name, "_snn"),
            "distances" = paste0(assay.name, "_nn"),
            graph_name
          )
          seurat_obj@graphs[[seurat_graph_name]] <- as.Graph(graph_matrix)
        }
      }, error = function(e) {
        if (verbose) warning("Could not add graph ", graph_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 9b. Add varp (pairwise variable annotations) to misc
  varp_path <- "varp"
  varp_names <- .zarr_select_keep(varp,
    if (.zarr_node_type(file, varp_path) == "group") {
      .zarr_list_children(file, varp_path)
    } else character(0))
  if (length(varp_names) > 0L) {
    if (verbose) message("Adding pairwise variable annotations (varp)...")
    varp_list <- list()
    for (varp_name in varp_names) {
      tryCatch({
        varp_mat <- .zarr_read_anndata_matrix(
          file, file.path(varp_path, varp_name), transpose = FALSE,
          obs_idx = var_idx, var_idx = var_idx
        )
        if (nrow(varp_mat) == length(feature.names) && ncol(varp_mat) == length(feature.names)) {
          rownames(varp_mat) <- feature.names
          colnames(varp_mat) <- feature.names
        }
        varp_list[[varp_name]] <- varp_mat
        if (verbose) message("  Added varp: ", varp_name)
      }, error = function(e) {
        if (verbose) warning("Could not read varp ", varp_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
    if (length(varp_list) > 0) {
      seurat_obj@misc[["__varp__"]] <- varp_list
    }
  }

  # 10. Add uns (unstructured) data to misc
  uns_path <- "uns"
  uns_children <- .zarr_select_keep(uns,
    if (.zarr_node_type(file, uns_path) == "group") {
      .zarr_list_children(file, uns_path)
    } else character(0))
  if (length(uns_children) > 0L) {
    if (verbose) message("Adding unstructured data...")
    for (item in uns_children) {
      item_path <- file.path(uns_path, item)
      tryCatch({
        node_type <- .zarr_node_type(file, item_path)
        if (node_type == "array") {
          # Try reading as numeric, fall back to string
          .uns_store <- .zarr_make_store(file)
          meta <- .zarr_read_array_meta(.uns_store, item_path)
          if (meta$is_string) {
            seurat_obj@misc[[item]] <- .zarr_read_strings(file, item_path)
          } else {
            seurat_obj@misc[[item]] <- .zarr_read_numeric(file, item_path)
          }
        } else if (node_type == "group") {
          if (verbose) message("  Storing complex uns item: ", item)
          seurat_obj@misc[[paste0(item, "_present")]] <- TRUE
        }
      }, error = function(e) {
        if (verbose) warning("Could not add uns item ", item, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  if (verbose) {
    message("\nSuccessfully loaded Zarr store")
    message("  Cells: ", ncol(seurat_obj))
    message("  Features: ", nrow(seurat_obj))
    message("  Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
    if (length(seurat_obj@reductions) > 0) {
      message("  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
    }
    if (length(seurat_obj@graphs) > 0) {
      message("  Graphs: ", paste(names(seurat_obj@graphs), collapse = ", "))
    }
    message("  Metadata columns: ", ncol(seurat_obj@meta.data))
  }

  return(seurat_obj)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal AnnData Zarr Reading Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read AnnData index from obs or var group
#'
#' @param store_path Path to zarr store
#' @param group "obs" or "var"
#'
#' @return Character vector of index values, or NULL
#'
#' @keywords internal
#'
.zarr_read_anndata_index <- function(store_path, group, slice_idx = NULL) {
  attrs <- .zarr_read_attrs(store_path, group)
  index_col <- attrs[["_index"]] %||% "_index"

  index_path <- file.path(group, index_col)
  if (.zarr_node_type(store_path, index_path) != "missing") {
    return(as.character(.zarr_read_strings(store_path, index_path,
                                            slice_idx = slice_idx)))
  }

  # Try "index" as fallback
  index_path <- file.path(group, "index")
  if (.zarr_node_type(store_path, index_path) != "missing") {
    return(as.character(.zarr_read_strings(store_path, index_path,
                                            slice_idx = slice_idx)))
  }

  NULL
}

#' Read an AnnData matrix (dense or sparse) from zarr store
#'
#' Dispatches based on encoding-type attribute.
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path to matrix node
#' @param transpose If TRUE (default), transpose cells x genes to genes x cells
#'
#' @return A matrix or dgCMatrix
#'
#' @keywords internal
#'
.zarr_read_anndata_matrix <- function(store_path, rel_path, transpose = TRUE,
                                       obs_idx = NULL, var_idx = NULL) {
  node_type <- .zarr_node_type(store_path, rel_path)

  if (node_type == "group") {
    attrs <- .zarr_read_attrs(store_path, rel_path)
    encoding_type <- attrs[["encoding-type"]] %||% ""

    if (encoding_type %in% c("csr_matrix", "csc_matrix")) {
      shape <- attrs[["shape"]]
      if (is.list(shape)) shape <- unlist(shape)
      shape <- as.integer(shape)

      if (is.null(obs_idx) && is.null(var_idx)) {
        # Full read (existing path)
        data <- .zarr_read_numeric(store_path, file.path(rel_path, "data"))
        indices <- .zarr_read_numeric(store_path, file.path(rel_path, "indices"))
        indptr <- .zarr_read_numeric(store_path, file.path(rel_path, "indptr"))
        if (encoding_type == "csc_matrix") {
          mat <- sparseMatrix(
            i = as.integer(indices) + 1L,
            p = as.integer(indptr),
            x = as.numeric(data),
            dims = shape, repr = "C"
          )
          if (transpose) mat <- t(mat)
          return(mat)
        }
        return(ReconstructSparseCSR(
          data = data, indices = indices, indptr = indptr,
          shape = shape, transpose = transpose
        ))
      }

      # Sliced read via indptr pushdown
      return(.zarr_read_sparse_sliced(store_path, rel_path, encoding_type,
                                       shape, transpose, obs_idx, var_idx))
    }

    # Group but not sparse - check for sub-arrays
    stop("Unsupported matrix encoding: ", encoding_type, call. = FALSE)
  }

  if (node_type == "array") {
    # Dense matrix. AnnData convention: rows = cells, cols = features.
    # obs_idx selects rows; var_idx selects cols (pre-transpose).
    slice <- NULL
    if (!is.null(obs_idx) || !is.null(var_idx)) {
      slice <- list(d1 = obs_idx, d2 = var_idx)
    }
    mat <- .zarr_read_numeric(store_path, rel_path, slice = slice)
    if (!is.matrix(mat)) {
      stop("Expected 2D array for expression matrix", call. = FALSE)
    }
    if (transpose) mat <- t(mat)
    return(mat)
  }

  stop("Matrix node not found at: ", rel_path, call. = FALSE)
}

#' Sparse-matrix slice via indptr pushdown
#'
#' AnnData stores sparse X as rows = cells / cols = features in CSR, or the
#' reverse in CSC. For CSR we use indptr to locate the data/indices ranges
#' for the requested rows (cells = obs_idx) and only fetch those chunks of
#' \code{data}/\code{indices}; columns (var_idx) are then filtered in
#' memory. For CSC the same trick works on var_idx; obs_idx falls back to
#' a memory subset.
#'
#' @keywords internal
.zarr_read_sparse_sliced <- function(store_path, rel_path, encoding_type,
                                      shape, transpose, obs_idx, var_idx) {
  store <- .zarr_make_store(store_path)

  # Always read the full indptr (small: 4*(n_outer+1) bytes uncompressed).
  indptr <- as.integer(.zarr_read_numeric(store,
                                          file.path(rel_path, "indptr")))

  if (encoding_type == "csr_matrix") {
    # outer = rows = cells; inner = cols = features
    n_obs <- shape[1]; n_var <- shape[2]
    if (is.logical(obs_idx)) obs_idx <- which(obs_idx)
    if (is.null(obs_idx)) obs_idx <- seq_len(n_obs)
    obs_idx <- as.integer(obs_idx)
    if (any(obs_idx < 1L | obs_idx > n_obs)) {
      stop("obs_idx out of range", call. = FALSE)
    }

    # For each selected row, the on-disk range is indptr[i]:indptr[i+1]-1
    # (indptr is 0-based offsets into the flat data/indices arrays).
    starts <- indptr[obs_idx]
    ends   <- indptr[obs_idx + 1L] - 1L
    nz_per_row <- ends - starts + 1L
    flat_idx <- unlist(mapply(seq.int, starts + 1L, ends + 1L,
                              SIMPLIFY = FALSE))  # 1-based

    if (length(flat_idx) == 0L) {
      sub <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                          dims = c(length(obs_idx), n_var), repr = "R")
    } else {
      sub_data    <- .zarr_read_numeric(store, file.path(rel_path, "data"),
                                         slice = list(d1 = flat_idx))
      sub_indices <- .zarr_read_numeric(store, file.path(rel_path, "indices"),
                                         slice = list(d1 = flat_idx))
      new_indptr <- as.integer(c(0L, cumsum(nz_per_row)))
      sub <- ReconstructSparseCSR(
        data = as.numeric(sub_data),
        indices = as.integer(sub_indices),
        indptr = new_indptr,
        shape = as.integer(c(length(obs_idx), n_var)),
        transpose = FALSE
      )
    }
    # Column filter in memory (CSR cannot push columns down efficiently).
    if (!is.null(var_idx)) {
      if (is.logical(var_idx)) var_idx <- which(var_idx)
      sub <- sub[, as.integer(var_idx), drop = FALSE]
    }
    if (transpose) sub <- t(sub)
    return(sub)
  }

  # CSC
  n_obs <- shape[1]; n_var <- shape[2]
  if (is.logical(var_idx)) var_idx <- which(var_idx)
  if (is.null(var_idx)) var_idx <- seq_len(n_var)
  var_idx <- as.integer(var_idx)
  if (any(var_idx < 1L | var_idx > n_var)) {
    stop("var_idx out of range", call. = FALSE)
  }

  starts <- indptr[var_idx]
  ends   <- indptr[var_idx + 1L] - 1L
  nz_per_col <- ends - starts + 1L
  flat_idx <- unlist(mapply(seq.int, starts + 1L, ends + 1L,
                            SIMPLIFY = FALSE))

  if (length(flat_idx) == 0L) {
    sub <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                        dims = c(n_obs, length(var_idx)), repr = "C")
  } else {
    sub_data    <- .zarr_read_numeric(store, file.path(rel_path, "data"),
                                       slice = list(d1 = flat_idx))
    sub_indices <- .zarr_read_numeric(store, file.path(rel_path, "indices"),
                                       slice = list(d1 = flat_idx))
    new_indptr <- as.integer(c(0L, cumsum(nz_per_col)))
    sub <- sparseMatrix(
      i = as.integer(sub_indices) + 1L,
      p = new_indptr,
      x = as.numeric(sub_data),
      dims = c(n_obs, length(var_idx)), repr = "C"
    )
  }
  if (!is.null(obs_idx)) {
    if (is.logical(obs_idx)) obs_idx <- which(obs_idx)
    sub <- sub[as.integer(obs_idx), , drop = FALSE]
  }
  if (transpose) sub <- t(sub)
  sub
}

#' Read an AnnData DataFrame column from zarr store
#'
#' Handles categorical, string-array, and numeric encodings.
#'
#' @param store_path Path to zarr store
#' @param col_path Full relative path to column (e.g., "obs/cell_type")
#'
#' @return Vector (factor, character, numeric, or logical)
#'
#' @keywords internal
#'
.zarr_read_anndata_column <- function(store_path, col_path, slice_idx = NULL) {
  node_type <- .zarr_node_type(store_path, col_path)
  slice1d <- if (is.null(slice_idx)) NULL else list(d1 = slice_idx)

  if (node_type == "group") {
    # Could be categorical encoding: group with codes/ and categories/ sub-arrays
    attrs <- .zarr_read_attrs(store_path, col_path)
    encoding_type <- attrs[["encoding-type"]] %||% ""

    if (encoding_type == "categorical") {
      store <- .zarr_make_store(store_path)
      codes <- .zarr_read_numeric(store, file.path(col_path, "codes"),
                                   slice = slice1d)
      # Categories are global (not row-aligned), always read in full.
      cats_path <- file.path(col_path, "categories")
      cats_meta <- .zarr_read_array_meta(store, cats_path)
      if (cats_meta$is_string) {
        categories <- .zarr_read_strings(store_path, cats_path)
      } else {
        categories <- as.character(.zarr_read_numeric(store_path, cats_path))
      }
      return(DecodeCategorical(as.integer(codes), as.character(categories)))
    }

    return(NULL)
  }

  if (node_type == "array") {
    store <- .zarr_make_store(store_path)
    meta  <- .zarr_read_array_meta(store, col_path)

    if (meta$is_string) {
      return(.zarr_read_strings(store_path, col_path, slice_idx = slice_idx))
    }

    values <- .zarr_read_numeric(store_path, col_path, slice = slice1d)

    if (meta$is_bool) {
      values <- as.logical(values)
    }

    return(values)
  }

  NULL
}

#' Get number of observations from zarr store
#'
#' @keywords internal
#'
.zarr_get_n_obs <- function(store_path) {
  # Try X shape
  x_type <- .zarr_node_type(store_path, "X")
  if (x_type == "group") {
    attrs <- .zarr_read_attrs(store_path, "X")
    shape <- attrs[["shape"]]
    if (!is.null(shape)) return(shape[[1]])
  }
  if (x_type == "array") {
    store <- .zarr_make_store(store_path)
    meta <- .zarr_read_array_meta(store, "X")
    return(meta$shape[[1]])
  }
  stop("Cannot determine number of observations", call. = FALSE)
}

#' Get number of variables from zarr store
#'
#' @keywords internal
#'
.zarr_get_n_vars <- function(store_path) {
  x_type <- .zarr_node_type(store_path, "X")
  if (x_type == "group") {
    attrs <- .zarr_read_attrs(store_path, "X")
    shape <- attrs[["shape"]]
    if (!is.null(shape)) return(shape[[2]])
  }
  if (x_type == "array") {
    store <- .zarr_make_store(store_path)
    meta <- .zarr_read_array_meta(store, "X")
    return(meta$shape[[2]])
  }
  stop("Cannot determine number of variables", call. = FALSE)
}
