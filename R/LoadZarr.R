#' @include zzz.R
#' @include AnnDataEncoding.R
#' @importFrom Matrix sparseMatrix t
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom SeuratObject Cells
NULL

#' Load an AnnData Zarr store as a Seurat object
#'
#' Direct conversion from AnnData Zarr format (.zarr directory) to Seurat object.
#' Supports zarr v2 stores as produced by Python anndata. This enables reading
#' of cloud-native AnnData stores from CELLxGENE, Human Cell Atlas, and
#' SpatialData workflows.
#'
#' @param file Path to .zarr directory
#' @param assay.name Name for the primary assay (default: "RNA")
#' @param verbose Show progress messages
#' @param ... Additional arguments (currently unused)
#'
#' @return A \code{Seurat} object
#'
#' @export
#'
readZarr <- function(file, assay.name = "RNA", verbose = TRUE, ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required for readZarr. ",
         "Install with: install.packages('jsonlite')", call. = FALSE)
  }

  if (!dir.exists(file)) {
    stop("Zarr store not found: ", file, call. = FALSE)
  }

  # Validate zarr store
  version <- .zarr_store_version(file)
  if (verbose) message("Loading Zarr store (v", version, "): ", file)

  # Read root attributes to verify it's an AnnData store
  root_attrs <- .zarr_read_attrs(file)

  # 1. Read cell names
  if (verbose) message("Reading cell names...")
  cell.names <- .zarr_read_anndata_index(file, "obs")

  # 2. Read feature names
  if (verbose) message("Reading feature names...")
  feature.names <- .zarr_read_anndata_index(file, "var")

  # Generate fallback names if needed
  if (is.null(cell.names)) {
    cell.names <- paste0("Cell", seq_len(.zarr_get_n_obs(file)))
  }
  if (is.null(feature.names)) {
    feature.names <- paste0("Gene", seq_len(.zarr_get_n_vars(file)))
  }

  # 3. Read expression matrix
  if (verbose) message("Reading expression matrix...")
  x_path <- "X"
  if (.zarr_node_type(file, x_path) == "missing") {
    stop("No expression matrix (X) found in zarr store", call. = FALSE)
  }

  expr_matrix <- .zarr_read_anndata_matrix(file, x_path)

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
  if (.zarr_node_type(file, layers_path) == "group") {
    if (verbose) message("Adding layers...")
    layer_names <- .zarr_list_children(file, layers_path)
    for (layer_name in layer_names) {
      if (verbose) message("  Adding layer: ", layer_name)
      tryCatch({
        layer_matrix <- .zarr_read_anndata_matrix(
          file, file.path(layers_path, layer_name)
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
        meta_values <- .zarr_read_anndata_column(file, col_path)
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
  if (.zarr_node_type(file, obsm_path) == "group") {
    if (verbose) message("Adding dimensional reductions...")
    reduc_names <- .zarr_list_children(file, obsm_path)

    for (reduc_name in reduc_names) {
      clean_name <- gsub("^X_", "", reduc_name)
      if (verbose) message("  Adding reduction: ", clean_name)

      tryCatch({
        embeddings <- .zarr_read_numeric(file, file.path(obsm_path, reduc_name))

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
        meta_values <- .zarr_read_anndata_column(file, col_path)
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
  if (.zarr_node_type(file, obsp_path) == "group") {
    if (verbose) message("Adding neighbor graphs...")
    graph_names <- .zarr_list_children(file, obsp_path)

    for (graph_name in graph_names) {
      if (verbose) message("  Adding graph: ", graph_name)
      tryCatch({
        graph_matrix <- .zarr_read_anndata_matrix(
          file, file.path(obsp_path, graph_name), transpose = FALSE
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
  if (.zarr_node_type(file, varp_path) == "group") {
    if (verbose) message("Adding pairwise variable annotations (varp)...")
    varp_list <- list()
    varp_names <- .zarr_list_children(file, varp_path)
    for (varp_name in varp_names) {
      tryCatch({
        varp_mat <- .zarr_read_anndata_matrix(
          file, file.path(varp_path, varp_name), transpose = FALSE
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
  if (.zarr_node_type(file, uns_path) == "group") {
    if (verbose) message("Adding unstructured data...")
    uns_children <- .zarr_list_children(file, uns_path)

    for (item in uns_children) {
      item_path <- file.path(uns_path, item)
      tryCatch({
        node_type <- .zarr_node_type(file, item_path)
        if (node_type == "array") {
          # Try reading as numeric, fall back to string
          meta <- .zarr_read_json(file.path(file, item_path, ".zarray"))
          if (!is.null(meta$dtype) && meta$dtype == "|O") {
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
.zarr_read_anndata_index <- function(store_path, group) {
  attrs <- .zarr_read_attrs(store_path, group)
  index_col <- attrs[["_index"]] %||% "_index"

  index_path <- file.path(group, index_col)
  if (.zarr_node_type(store_path, index_path) != "missing") {
    return(as.character(.zarr_read_strings(store_path, index_path)))
  }

  # Try "index" as fallback
  index_path <- file.path(group, "index")
  if (.zarr_node_type(store_path, index_path) != "missing") {
    return(as.character(.zarr_read_strings(store_path, index_path)))
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
.zarr_read_anndata_matrix <- function(store_path, rel_path, transpose = TRUE) {
  node_type <- .zarr_node_type(store_path, rel_path)

  if (node_type == "group") {
    attrs <- .zarr_read_attrs(store_path, rel_path)
    encoding_type <- attrs[["encoding-type"]] %||% ""

    if (encoding_type %in% c("csr_matrix", "csc_matrix")) {
      # Sparse matrix
      data <- .zarr_read_numeric(store_path, file.path(rel_path, "data"))
      indices <- .zarr_read_numeric(store_path, file.path(rel_path, "indices"))
      indptr <- .zarr_read_numeric(store_path, file.path(rel_path, "indptr"))
      shape <- attrs[["shape"]]
      if (is.list(shape)) shape <- unlist(shape)

      if (encoding_type == "csc_matrix") {
        # CSC: indptr = column pointers, indices = row indices (0-based)
        # Directly construct dgCMatrix (R's native CSC format)
        mat <- sparseMatrix(
          i = as.integer(indices) + 1L,
          p = as.integer(indptr),
          x = as.numeric(data),
          dims = as.integer(shape),
          repr = "C"
        )
        if (transpose) mat <- t(mat)
        return(mat)
      }

      return(ReconstructSparseCSR(
        data = data,
        indices = indices,
        indptr = indptr,
        shape = as.integer(shape),
        transpose = transpose
      ))
    }

    # Group but not sparse - check for sub-arrays
    stop("Unsupported matrix encoding: ", encoding_type, call. = FALSE)
  }

  if (node_type == "array") {
    # Dense matrix
    mat <- .zarr_read_numeric(store_path, rel_path)
    if (!is.matrix(mat)) {
      stop("Expected 2D array for expression matrix", call. = FALSE)
    }
    if (transpose) mat <- t(mat)
    return(mat)
  }

  stop("Matrix node not found at: ", rel_path, call. = FALSE)
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
.zarr_read_anndata_column <- function(store_path, col_path) {
  node_type <- .zarr_node_type(store_path, col_path)

  if (node_type == "group") {
    # Could be categorical encoding: group with codes/ and categories/ sub-arrays
    attrs <- .zarr_read_attrs(store_path, col_path)
    encoding_type <- attrs[["encoding-type"]] %||% ""

    if (encoding_type == "categorical") {
      codes <- .zarr_read_numeric(store_path, file.path(col_path, "codes"))
      # Categories may be strings or numeric
      cats_path <- file.path(col_path, "categories")
      cats_meta <- .zarr_read_json(file.path(store_path, cats_path, ".zarray"))
      if (!is.null(cats_meta$dtype) && cats_meta$dtype == "|O") {
        categories <- .zarr_read_strings(store_path, cats_path)
      } else {
        categories <- as.character(.zarr_read_numeric(store_path, cats_path))
      }
      return(DecodeCategorical(as.integer(codes), as.character(categories)))
    }

    return(NULL)
  }

  if (node_type == "array") {
    meta <- .zarr_read_json(file.path(store_path, col_path, ".zarray"))

    if (!is.null(meta$dtype) && meta$dtype == "|O") {
      # String array
      return(.zarr_read_strings(store_path, col_path))
    }

    # Numeric or boolean array
    values <- .zarr_read_numeric(store_path, col_path)

    # Check for boolean encoding
    if (!is.null(meta$dtype) && grepl("^[|<>]?b1$", meta$dtype)) {
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
    meta <- .zarr_read_json(file.path(store_path, "X", ".zarray"))
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
    meta <- .zarr_read_json(file.path(store_path, "X", ".zarray"))
    return(meta$shape[[2]])
  }
  stop("Cannot determine number of variables", call. = FALSE)
}
