#' Load an AnnData H5AD file as a Seurat object
#'
#' Direct conversion from H5AD format to Seurat object without intermediate h5Seurat.
#' Supports optional BPCells on-disk matrix loading for large datasets that exceed
#' available memory. When compiled C routines are available, uses a fast native
#' reader (typically 2-3x faster than the pure-R path).
#'
#' @param file Path to H5AD file
#' @param assay.name Name for the primary assay (default: "RNA")
#' @param use.bpcells If not NULL, a directory path where BPCells will store the
#'   expression matrix on disk. Requires the BPCells package. The resulting Seurat
#'   object will reference the on-disk matrix instead of loading it into memory,
#'   enabling analysis of datasets larger than available RAM.
#' @param components Character vector of h5ad components to load. Default loads
#'   everything. Use \code{c("X")} for matrix-only (fastest), or any subset of
#'   \code{c("X", "obs", "var", "obsm", "obsp", "varp", "layers", "uns")}.
#'   The file path is stored in \code{misc[[".__h5ad_path__"]]} for deferred
#'   loading via \code{\link{scLoadMeta}}.
#' @param use.c Use compiled C reader when available (default: TRUE). Set to
#'   FALSE to force the pure-R hdf5r path.
#' @param verbose Show progress messages
#'
#' @return A \code{Seurat} object. If \code{use.bpcells} is set, the count matrix
#'   is stored on disk in BPCells format and the object uses minimal memory.
#'
#' @importFrom hdf5r H5File h5attr h5attr_names
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom SeuratObject Cells
#'
#' @export
#'
readH5AD <- function(file, assay.name = "RNA", use.bpcells = NULL,
                     components = NULL, use.c = TRUE, verbose = TRUE) {
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  # Normalize components parameter
  all_components <- c("X", "obs", "var", "obsm", "obsp", "varp", "layers", "uns")
  if (is.null(components)) {
    components <- all_components
  } else {
    components <- match.arg(components, all_components, several.ok = TRUE)
    if (!"X" %in% components) components <- c("X", components)
  }

  # Try C reader (Optimization 1: ~15% faster than hdf5r path for sparse data)
  c_available <- use.c && is.null(use.bpcells) && is.loaded("C_read_h5ad", PACKAGE = "scConvert")
  if (c_available) {
    c_result <- tryCatch(
      .readH5AD_c(file, assay.name = assay.name, components = components,
                   verbose = verbose),
      error = function(e) {
        if (verbose) message("C reader failed (", conditionMessage(e), "), using R reader")
        NULL
      }
    )
    if (!is.null(c_result)) return(c_result)
  }

  h5ad <- H5File$new(file, mode = "r")
  on.exit(h5ad$close_all())

  if (verbose) {
    message("Loading H5AD file: ", file)
  }

  # Helper function to read H5AD sparse or dense matrix
  ReadH5ADMatrix <- function(h5_obj, transpose = TRUE) {
    if (inherits(h5_obj, "H5Group")) {
      # Sparse matrix (CSR or CSC format in h5ad)
      if (h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
        data_vals <- h5_obj[["data"]][]
        indices <- h5_obj[["indices"]][]   # 0-based
        indptr <- h5_obj[["indptr"]][]

        # Detect encoding type: CSR vs CSC
        encoding <- tryCatch(h5attr(h5_obj, "encoding-type"), error = function(e) "csr_matrix")
        is_csc <- identical(encoding, "csc_matrix")

        # Get dimensions from shape attribute
        if (h5_obj$attr_exists("shape")) {
          shape <- h5attr(h5_obj, "shape")
          n_rows <- shape[1]
          n_cols <- shape[2]
        } else if (is_csc) {
          n_cols <- length(indptr) - 1L
          n_rows <- if (length(indices) > 0) max(indices) + 1L else 0L
        } else {
          n_rows <- length(indptr) - 1L
          n_cols <- if (length(indices) > 0) max(indices) + 1L else 0L
        }

        # Helper: sort row indices within each column for valid dgCMatrix
        # dgCMatrix requires @i to be increasing within each column (defined by @p).
        # scipy CSR/CSC matrices may have unsorted indices (e.g. scanpy pbmc3k raw/X).
        .sort_dgc_indices <- function(i, p, x) {
          needs_sort <- FALSE
          n_col <- length(p) - 1L
          for (ci in seq_len(n_col)) {
            start <- p[ci] + 1L
            end <- p[ci + 1L]
            if (end > start && is.unsorted(i[start:end])) {
              needs_sort <- TRUE
              break
            }
          }
          if (!needs_sort) return(list(i = i, x = x))
          # Sort indices within each column
          for (ci in seq_len(n_col)) {
            start <- p[ci] + 1L
            end <- p[ci + 1L]
            if (end > start) {
              seg <- start:end
              ord <- order(i[seg])
              i[seg] <- i[seg][ord]
              x[seg] <- x[seg][ord]
            }
          }
          list(i = i, x = x)
        }

        if (is_csc) {
          # CSC format: indptr = column pointers, indices = row indices
          # dgCMatrix is natively CSC, so construct directly
          sorted <- .sort_dgc_indices(as.integer(indices), as.integer(indptr),
                                       as.numeric(data_vals))
          mat <- new("dgCMatrix",
            i = sorted$i,
            p = as.integer(indptr),
            x = sorted$x,
            Dim = c(as.integer(n_rows), as.integer(n_cols))
          )
          if (transpose) {
            mat <- Matrix::t(mat)
          }
        } else if (transpose) {
          # CSR->CSC reinterpretation: CSR of (n_rows x n_cols) == CSC of (n_cols x n_rows)
          # h5ad CSR indices become dgCMatrix @i (row indices in transposed view).
          # scipy CSR may have unsorted column indices within rows, so we must
          # sort after reinterpretation to produce a valid dgCMatrix.
          i_int <- as.integer(indices)
          p_int <- as.integer(indptr)
          x_num <- as.numeric(data_vals)
          sorted <- .sort_dgc_indices(i_int, p_int, x_num)
          mat <- new("dgCMatrix",
            i = sorted$i,
            p = p_int,
            x = sorted$x,
            Dim = c(as.integer(n_cols), as.integer(n_rows))
          )
        } else {
          # Keep original CSR orientation (e.g. obsp graphs)
          indices_1based <- indices + 1L
          row_indices <- rep(seq_len(n_rows), diff(indptr))
          mat <- sparseMatrix(
            i = row_indices,
            j = indices_1based,
            x = data_vals,
            dims = c(n_rows, n_cols),
            index1 = TRUE
          )
        }
        return(mat)
      }
    } else if (inherits(h5_obj, "H5D")) {
      # Dense matrix
      mat <- h5_obj[,]
      if (transpose) {
        mat <- t(mat)
      }
      return(mat)
    }
    stop("Unknown matrix format", call. = FALSE)
  }

  # 1. Read cell names
  if (verbose) message("Reading cell names...")
  cell.names <- NULL
  obs_compound_df <- NULL
  if (h5ad$exists("obs")) {
    obs_obj <- h5ad[["obs"]]
    if (inherits(obs_obj, "H5Group")) {
      # Modern h5ad format: obs is a group with _index dataset or _index attribute
      if (obs_obj$exists("_index")) {
        cell.names <- as.character(obs_obj[["_index"]][])
      } else if (obs_obj$exists("index")) {
        cell.names <- as.character(obs_obj[["index"]][])
      } else if (obs_obj$attr_exists("_index")) {
        # AnnData convention: _index attribute names the index column
        idx_col <- h5attr(obs_obj, "_index")
        if (obs_obj$exists(idx_col)) {
          cell.names <- as.character(obs_obj[[idx_col]][])
        }
      }
    } else if (inherits(obs_obj, "H5D")) {
      # Legacy h5ad format: obs is a compound HDF5 dataset
      obs_compound_df <- obs_obj$read()
      if ("index" %in% names(obs_compound_df)) {
        cell.names <- as.character(obs_compound_df$index)
      } else if ("_index" %in% names(obs_compound_df)) {
        cell.names <- as.character(obs_compound_df[["_index"]])
      } else if (!is.null(rownames(obs_compound_df))) {
        cell.names <- rownames(obs_compound_df)
      }
    }
  }

  # Generate cell names if not found
  if (is.null(cell.names)) {
    n_cells <- if (h5ad$exists("X")) {
      if (inherits(h5ad[["X"]], "H5Group") && h5ad[["X"]]$attr_exists("shape")) {
        h5attr(h5ad[["X"]], "shape")[1]
      } else if (inherits(h5ad[["X"]], "H5D")) {
        h5ad[["X"]]$dims[1]
      }
    } else {
      stop("Cannot determine number of cells", call. = FALSE)
    }
    cell.names <- paste0("Cell", seq_len(n_cells))
  }

  # 2. Read feature names
  if (verbose) message("Reading feature names...")
  feature.names <- NULL
  var_compound_df <- NULL
  if (h5ad$exists("var")) {
    var_obj <- h5ad[["var"]]
    if (inherits(var_obj, "H5Group")) {
      if (var_obj$exists("_index")) {
        feature.names <- as.character(var_obj[["_index"]][])
      } else if (var_obj$exists("index")) {
        feature.names <- as.character(var_obj[["index"]][])
      } else if (var_obj$attr_exists("_index")) {
        idx_col <- h5attr(var_obj, "_index")
        if (var_obj$exists(idx_col)) {
          feature.names <- as.character(var_obj[[idx_col]][])
        }
      }
    } else if (inherits(var_obj, "H5D")) {
      # Legacy compound dataset
      var_compound_df <- var_obj$read()
      if ("index" %in% names(var_compound_df)) {
        feature.names <- as.character(var_compound_df$index)
      } else if ("_index" %in% names(var_compound_df)) {
        feature.names <- as.character(var_compound_df[["_index"]])
      }
    }
  }

  # Generate feature names if not found
  if (is.null(feature.names)) {
    n_features <- if (h5ad$exists("X")) {
      if (inherits(h5ad[["X"]], "H5Group") && h5ad[["X"]]$attr_exists("shape")) {
        h5attr(h5ad[["X"]], "shape")[2]
      } else if (inherits(h5ad[["X"]], "H5D")) {
        h5ad[["X"]]$dims[2]
      }
    } else {
      stop("Cannot determine number of features", call. = FALSE)
    }
    feature.names <- paste0("Gene", seq_len(n_features))
  }

  # Deduplicate feature names (some datasets have duplicates, e.g. squidpy four_i)
  if (anyDuplicated(feature.names)) {
    feature.names <- make.unique(feature.names)
  }

  # 3. Read main expression matrix
  if (verbose) message("Reading expression matrix...")
  if (!h5ad$exists("X")) {
    stop("No expression matrix (X) found in h5ad file", call. = FALSE)
  }

  use_bpcells <- !is.null(use.bpcells)
  if (use_bpcells) {
    if (!requireNamespace("BPCells", quietly = TRUE)) {
      stop("BPCells package is required for on-disk loading. ",
           "Install with: remotes::install_github('bnprks/BPCells')", call. = FALSE)
    }

    # Prefer raw/X (raw counts) over X (often normalized) -- matches non-BPCells path
    has_raw <- h5ad$exists("raw") && h5ad[["raw"]]$exists("X")
    bp_group <- if (has_raw) "raw/X" else "X"
    if (verbose) {
      if (has_raw) {
        message("Loading raw counts (raw/X) via BPCells (on-disk)...")
      } else {
        message("Loading expression matrix (X) via BPCells (on-disk)...")
      }
    }

    # Close hdf5r handle before BPCells opens the file (avoid lock conflicts)
    h5ad$close_all()

    if (isTRUE(use.bpcells)) {
      # use.bpcells = TRUE: read directly from h5ad (backed by HDF5)
      bp_ns <- asNamespace("BPCells")
      expr_matrix <- get("open_matrix_anndata_hdf5", envir = bp_ns)(path = file, group = bp_group)
    } else {
      # use.bpcells = directory path: write to BPCells dir for faster repeated access
      bp_ns <- asNamespace("BPCells")
      bpcells_dir <- file.path(use.bpcells, "counts")
      dir.create(bpcells_dir, recursive = TRUE, showWarnings = FALSE)
      bpcells_mat <- get("open_matrix_anndata_hdf5", envir = bp_ns)(path = file, group = bp_group)
      get("write_matrix_dir", envir = bp_ns)(mat = bpcells_mat, dir = bpcells_dir, overwrite = TRUE)
      expr_matrix <- get("open_matrix_dir", envir = bp_ns)(dir = bpcells_dir)
    }
    # Reopen h5ad for metadata reading
    h5ad <- H5File$new(file, mode = "r")
  } else {
    expr_matrix <- ReadH5ADMatrix(h5ad[["X"]], transpose = TRUE)
  }

  # Handle dimension mismatches: check if matrix needs transposing
  # This can happen with dense matrices where R/Python storage conventions differ
  n_features <- length(feature.names)
  n_cells <- length(cell.names)
  if (nrow(expr_matrix) == n_cells && ncol(expr_matrix) == n_features &&
      nrow(expr_matrix) != n_features) {
    # Matrix is in (cells x genes) orientation, transpose to (genes x cells)
    expr_matrix <- t(expr_matrix)
  }

  if (nrow(expr_matrix) != n_features) {
    warning(sprintf(
      "Feature count mismatch: matrix has %d rows but %d feature names found. Truncating feature names to match matrix.",
      nrow(expr_matrix), n_features
    ))
    feature.names <- feature.names[seq_len(nrow(expr_matrix))]
  }
  if (ncol(expr_matrix) != n_cells) {
    warning(sprintf(
      "Cell count mismatch: matrix has %d columns but %d cell names found. Truncating cell names to match matrix.",
      ncol(expr_matrix), n_cells
    ))
    cell.names <- cell.names[seq_len(ncol(expr_matrix))]
  }

  rownames(expr_matrix) <- feature.names
  colnames(expr_matrix) <- cell.names

  # 4. Create Seurat object
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(
    counts = expr_matrix,
    project = "H5AD",
    assay = assay.name,
    min.cells = 0,
    min.features = 0
  )

  # 5. Add raw counts if present
  # In scanpy convention: X = normalized/processed, raw/X = raw counts
  # When raw/X exists, X should become the "data" layer and raw/X the "counts" layer
  if (h5ad$exists("raw") && h5ad[["raw"]]$exists("X")) {
    if (use_bpcells) {
      # In BPCells mode, skip raw/X to preserve on-disk matrix.
      # Loading raw/X in-memory would defeat the purpose of on-disk mode.
      if (verbose) message("Skipping raw counts (BPCells on-disk mode preserves X as counts)")
    } else {
      if (verbose) message("Adding raw counts...")

      raw_features <- NULL
      if (h5ad[["raw"]]$exists("var")) {
        raw_var <- h5ad[["raw/var"]]
        if (raw_var$exists("_index")) {
          raw_features <- as.character(raw_var[["_index"]][])
        } else if (raw_var$exists("index")) {
          raw_features <- as.character(raw_var[["index"]][])
        }
      }

      if (!is.null(raw_features)) {
        raw_matrix <- ReadH5ADMatrix(h5ad[["raw/X"]], transpose = TRUE)

        # Handle dimension mismatches (dense matrices may need additional transpose)
        n_raw_features <- length(raw_features)
        if (nrow(raw_matrix) == n_cells && ncol(raw_matrix) == n_raw_features &&
            nrow(raw_matrix) != n_raw_features) {
          raw_matrix <- t(raw_matrix)
        }

        # Match dimensions
        raw_features <- raw_features[seq_len(min(length(raw_features), nrow(raw_matrix)))]
        rownames(raw_matrix) <- raw_features
        colnames(raw_matrix) <- cell.names

        # Find common features
        common_features <- intersect(feature.names, raw_features)
        if (length(common_features) > 0) {
          raw_subset <- raw_matrix[common_features, , drop = FALSE]
          # raw/X -> counts layer (actual raw counts)
          seurat_obj[[assay.name]] <- SetAssayData(
            object = seurat_obj[[assay.name]],
            layer = "counts",
            new.data = raw_subset
          )
          # X (already loaded as counts in step 4) -> data layer (normalized)
          # expr_matrix contains the X values which are normalized when raw exists
          x_subset <- expr_matrix[common_features, , drop = FALSE]
          seurat_obj[[assay.name]] <- SetAssayData(
            object = seurat_obj[[assay.name]],
            layer = "data",
            new.data = x_subset
          )
          if (verbose) message("  Set raw/X as counts, X as data (normalized)")
        }
      }
    }
  }

  # 6. Add layers if present
  if (h5ad$exists("layers")) {
    if (use_bpcells) {
      if (verbose) message("Skipping additional layers (BPCells on-disk mode)")
    } else {
      if (verbose) message("Adding layers...")
      layer_names <- names(h5ad[["layers"]])

      for (layer_name in layer_names) {
        if (verbose) message("  Adding layer: ", layer_name)

        layer_matrix <- ReadH5ADMatrix(h5ad[["layers"]][[layer_name]], transpose = TRUE)

        # Ensure dimensions match
        if (nrow(layer_matrix) == nrow(expr_matrix) && ncol(layer_matrix) == ncol(expr_matrix)) {
          rownames(layer_matrix) <- feature.names
          colnames(layer_matrix) <- cell.names

          # Map layer names to Seurat slots
          seurat_slot <- switch(layer_name,
            "counts" = "counts",
            "data" = "data",
            "log_normalized" = "data",
            "scale.data" = "scale.data",
            "scaled" = "scale.data",
            layer_name
          )

          tryCatch({
            seurat_obj[[assay.name]] <- SetAssayData(
              object = seurat_obj[[assay.name]],
              layer = seurat_slot,
              new.data = layer_matrix
            )
          }, error = function(e) {
            if (verbose) message("Could not add layer ", layer_name, ": ", e$message)
          })
        }
      }
    }
  }

  # 7. Add cell metadata
  if (h5ad$exists("obs")) {
    if (verbose) message("Adding cell metadata...")

    if (!is.null(obs_compound_df)) {
      # Legacy compound dataset: obs was already read as a data.frame
      meta_cols <- setdiff(names(obs_compound_df), c("_index", "index"))
      for (col in meta_cols) {
        tryCatch({
          seurat_obj[[col]] <- obs_compound_df[[col]]
        }, error = function(e) {
          if (verbose) message("Could not add metadata column '", col, "': ", e$message)
        })
      }
    } else {
      obs_group <- h5ad[["obs"]]
      # Exclude index column (may be _index, index, or named by _index attribute)
      obs_exclude <- c("_index", "index", "__categories")
      if (obs_group$attr_exists("_index")) {
        obs_exclude <- c(obs_exclude, h5attr(obs_group, "_index"))
      }
      obs_cols <- setdiff(names(obs_group), obs_exclude)

      # Cache __categories group reference if it exists (legacy format)
      has_legacy_cats <- obs_group$exists("__categories")
      legacy_cats <- if (has_legacy_cats) obs_group[["__categories"]] else NULL
      legacy_cat_names <- if (has_legacy_cats) names(legacy_cats) else character(0)

      # Pre-classify columns into groups (modern categoricals) vs datasets
      # using ls() to get all member types in a single HDF5 call
      obs_ls <- obs_group$ls()
      obs_types <- setNames(as.character(obs_ls$obj_type), obs_ls$name)
      cat_cols <- obs_cols[obs_cols %in% names(obs_types) & obs_types[obs_cols] == "H5I_GROUP"]
      plain_cols <- obs_cols[obs_cols %in% names(obs_types) & obs_types[obs_cols] == "H5I_DATASET"]

      # Batch-read all columns into a local data frame first, then assign once
      # (avoids per-column Seurat validation overhead which is O(n_cols * n_cells))
      obs_batch <- list()

      for (col in plain_cols) {
        tryCatch({
          col_obj <- obs_group[[col]]
          if (has_legacy_cats && col %in% legacy_cat_names) {
            codes <- col_obj$read()
            categories <- as.character(legacy_cats[[col]]$read())
            obs_batch[[col]] <- DecodeCategorical(codes, categories)
          } else {
            meta_values <- col_obj$read()
            if (!is.null(meta_values)) {
              obs_batch[[col]] <- meta_values
            }
          }
        }, error = function(e) {
          if (verbose) message("Could not add metadata column '", col, "': ", e$message)
        })
      }

      # Read modern categoricals (groups with codes/categories)
      for (col in cat_cols) {
        tryCatch({
          col_obj <- obs_group[[col]]
          encoding_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
          if (encoding_type == "categorical" && col_obj$exists("categories") && col_obj$exists("codes")) {
            codes <- col_obj[["codes"]]$read()
            categories <- as.character(col_obj[["categories"]]$read())
            obs_batch[[col]] <- DecodeCategorical(codes, categories)
          }
        }, error = function(e) {
          if (verbose) message("Could not add metadata column '", col, "': ", e$message)
        })
      }

      # Single bulk assignment (one Seurat validation pass instead of N)
      if (length(obs_batch) > 0) {
        batch_df <- data.frame(row.names = colnames(seurat_obj))
        for (col in names(obs_batch)) {
          batch_df[[col]] <- obs_batch[[col]]
        }
        seurat_obj <- AddMetaData(seurat_obj, metadata = batch_df)
      }
    }
  }

  # 8. Add dimensional reductions
  if (h5ad$exists("obsm")) {
    if (verbose) message("Adding dimensional reductions...")

    obsm_obj <- h5ad[["obsm"]]

    if (inherits(obsm_obj, "H5D")) {
      # Legacy compound dataset: obsm is a structured array with named fields
      obsm_df <- obsm_obj$read()
      n_cells <- length(cell.names)
      for (reduc_name in names(obsm_df)) {
        clean_name <- gsub("^X_", "", reduc_name)
        if (verbose) message("  Adding reduction: ", clean_name)
        tryCatch({
          vals <- obsm_df[[reduc_name]]
          n_dims <- length(vals) %/% n_cells
          if (n_dims < 1 || length(vals) != n_cells * n_dims) next
          embeddings <- matrix(vals, nrow = n_cells, ncol = n_dims, byrow = FALSE)
          rownames(embeddings) <- cell.names
          key <- AnnDataReductionKey(clean_name)
          colnames(embeddings) <- paste0(gsub("_$", "", key), "_", seq_len(n_dims))
          reduc_obj <- CreateDimReducObject(embeddings = embeddings, key = key, assay = assay.name)
          seurat_obj[[clean_name]] <- reduc_obj
        }, error = function(e) {
          if (verbose) message("Could not add reduction ", clean_name, ": ", e$message)
        })
      }
    } else {
      # Modern format: obsm is an H5Group with named datasets
      for (reduc_name in names(obsm_obj)) {
        clean_name <- gsub("^X_", "", reduc_name)
        # Skip 'spatial' -- handled separately in step 12
        if (clean_name == "spatial") next
        if (verbose) message("  Adding reduction: ", clean_name)

        tryCatch({
          reduc_item <- obsm_obj[[reduc_name]]

          # Handle both dense datasets and sparse groups
          if (inherits(reduc_item, "H5Group")) {
            embeddings <- as.matrix(ReadH5ADMatrix(reduc_item, transpose = FALSE))
          } else {
            embeddings <- reduc_item[,]
          }

          # HDF5 stores in row-major (C) order, R reads in column-major (Fortran) order
          # h5ad obsm shape is (n_obs, n_dims) but hdf5r may return (n_dims, n_obs)
          if (nrow(embeddings) != length(cell.names) && ncol(embeddings) == length(cell.names)) {
            embeddings <- t(embeddings)
          }

          # Check dimensions
          if (nrow(embeddings) != length(cell.names)) {
            if (verbose) message("Skipping reduction ", clean_name, " - dimension mismatch")
            next
          }

          rownames(embeddings) <- cell.names
          key <- AnnDataReductionKey(clean_name)
          colnames(embeddings) <- paste0(gsub("_$", "", key), "_", seq_len(ncol(embeddings)))

          reduc_obj <- CreateDimReducObject(
            embeddings = embeddings,
            key = key,
            assay = assay.name
          )

          seurat_obj[[clean_name]] <- reduc_obj
        }, error = function(e) {
          if (verbose) message("Could not add reduction ", clean_name, ": ", e$message)
        })
      }
    }
  }

  # 9. Add feature metadata
  if (h5ad$exists("var")) {
    if (verbose) message("Adding feature metadata...")

    if (!is.null(var_compound_df)) {
      # Legacy compound dataset: var was already read as a data.frame
      var_cols <- setdiff(names(var_compound_df), c("_index", "index"))
      for (col in var_cols) {
        if (verbose) message("  Adding feature metadata: ", col)
        tryCatch({
          meta_values <- var_compound_df[[col]]

          # Special handling for highly_variable
          if (col == "highly_variable") {
            if (is.numeric(meta_values)) meta_values <- as.logical(meta_values)
          }

          if (length(meta_values) == nrow(seurat_obj)) {
            names(meta_values) <- feature.names
            seurat_obj[[assay.name]][[col]] <- meta_values
            if (col == "highly_variable" && is.logical(meta_values)) {
              VariableFeatures(seurat_obj) <- feature.names[meta_values]
            }
          }
        }, error = function(e) {
          if (verbose) message("Could not add feature metadata ", col, ": ", e$message)
        })
      }
    } else {
    var_group <- h5ad[["var"]]
    # Exclude index column (may be _index, index, or named by _index attribute)
    var_exclude <- c("_index", "index", "__categories")
    if (var_group$attr_exists("_index")) {
      var_exclude <- c(var_exclude, h5attr(var_group, "_index"))
    }
    var_cols <- setdiff(names(var_group), var_exclude)

    # Cache __categories if present (legacy format)
    has_var_cats <- var_group$exists("__categories")
    var_cats <- if (has_var_cats) var_group[["__categories"]] else NULL
    var_cat_names <- if (has_var_cats) names(var_cats) else character(0)

    for (col in var_cols) {
      if (verbose) message("  Adding feature metadata: ", col)

      tryCatch({
        meta_values <- NULL
        col_obj <- var_group[[col]]

        if (inherits(col_obj, "H5Group")) {
          # Modern h5ad categorical format
          encoding_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
          if (encoding_type == "categorical" && col_obj$exists("categories") && col_obj$exists("codes")) {
            codes <- col_obj[["codes"]]$read()
            categories <- as.character(col_obj[["categories"]]$read())
            meta_values <- DecodeCategorical(codes, categories)
          }
        } else if (inherits(col_obj, "H5D")) {
          # Check legacy categorical format
          if (has_var_cats && col %in% var_cat_names) {
            codes <- col_obj$read()
            categories <- as.character(var_cats[[col]]$read())
            meta_values <- DecodeCategorical(codes, categories)
          } else {
            meta_values <- col_obj$read()
          }
        }

        if (is.null(meta_values)) next

        # Special handling for highly_variable: convert to logical
        if (col == "highly_variable") {
          if (is.factor(meta_values)) {
            # Categorical "True"/"False" strings
            meta_values <- as.character(meta_values) == "True"
          } else if (is.numeric(meta_values)) {
            meta_values <- as.logical(meta_values)
          }
        }

        # Ensure length matches and name with feature names for Seurat v5 compatibility
        if (length(meta_values) == nrow(seurat_obj)) {
          names(meta_values) <- feature.names
          seurat_obj[[assay.name]][[col]] <- meta_values

          # Set variable features if highly_variable column exists
          if (col == "highly_variable" && is.logical(meta_values)) {
            VariableFeatures(seurat_obj) <- feature.names[meta_values]
          }
        }
      }, error = function(e) {
        if (verbose) message("Could not add feature metadata ", col, ": ", e$message)
      })
    }
    } # end else (non-compound var)
  }

  # 10. Add neighbor graphs from obsp
  if (h5ad$exists("obsp")) {
    if (verbose) message("Adding neighbor graphs...")

    for (graph_name in names(h5ad[["obsp"]])) {
      if (verbose) message("  Adding graph: ", graph_name)

      # obsp stores CSR; read as-is then transpose to get correct orientation
      graph_matrix <- Matrix::t(ReadH5ADMatrix(h5ad[["obsp"]][[graph_name]], transpose = FALSE))

      # Ensure square matrix with correct dimensions
      if (nrow(graph_matrix) == ncol(graph_matrix) && nrow(graph_matrix) == ncol(seurat_obj)) {
        rownames(graph_matrix) <- cell.names
        colnames(graph_matrix) <- cell.names

        # Map to Seurat graph names
        seurat_graph_name <- switch(graph_name,
          "connectivities" = paste0(assay.name, "_snn"),
          "distances" = paste0(assay.name, "_nn"),
          graph_name
        )

        graph_obj <- as.Graph(graph_matrix)
        DefaultAssay(graph_obj) <- assay.name
        seurat_obj@graphs[[seurat_graph_name]] <- graph_obj
      }
    }
  }

  # 10b. Add varp (pairwise variable annotations) to misc
  if (h5ad$exists("varp")) {
    if (verbose) message("Adding pairwise variable annotations (varp)...")
    varp_list <- list()
    for (varp_name in names(h5ad[["varp"]])) {
      tryCatch({
        varp_mat <- ReadH5ADMatrix(h5ad[["varp"]][[varp_name]], transpose = FALSE)
        if (nrow(varp_mat) == length(feature.names) && ncol(varp_mat) == length(feature.names)) {
          rownames(varp_mat) <- feature.names
          colnames(varp_mat) <- feature.names
        }
        varp_list[[varp_name]] <- varp_mat
        if (verbose) message("  Added varp: ", varp_name,
                             " (", nrow(varp_mat), " x ", ncol(varp_mat), ")")
      }, error = function(e) {
        if (verbose) message("  Could not read varp item ", varp_name, ": ", e$message)
      })
    }
    if (length(varp_list) > 0) {
      seurat_obj@misc[["__varp__"]] <- varp_list
    }
  }

  # 11. Add uns (unstructured) data to misc
  if (h5ad$exists("uns")) {
    if (verbose) message("Adding unstructured data...")
    uns_group <- h5ad[["uns"]]

    .read_uns_group <- function(grp) {
      result <- list()
      for (item in names(grp)) {
        tryCatch({
          if (inherits(grp[[item]], "H5D")) {
            result[[item]] <- grp[[item]]$read()
          } else if (inherits(grp[[item]], "H5Group")) {
            result[[item]] <- .read_uns_group(grp[[item]])
          }
        }, error = function(e) {
          if (verbose) message("Could not read uns item '", item, "': ", e$message)
        })
      }
      result
    }

    for (item in names(uns_group)) {
      tryCatch({
        if (inherits(uns_group[[item]], "H5D")) {
          seurat_obj@misc[[item]] <- uns_group[[item]]$read()
        } else if (inherits(uns_group[[item]], "H5Group")) {
          if (verbose) message("  Reading complex uns item: ", item)
          seurat_obj@misc[[item]] <- .read_uns_group(uns_group[[item]])
        }
      }, error = function(e) {
        if (verbose) message("Could not add uns item ", item, ": ", e$message)
      })
    }
  }

  # 12. Add spatial data support
  if (h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
    seurat_obj <- H5ADSpatialToSeurat(
      h5ad_file = h5ad,
      seurat_obj = seurat_obj,
      assay_name = assay.name,
      verbose = verbose
    )
  }

  # Store source path for deferred loading (Optimization 4)
  if (!setequal(components, all_components)) {
    seurat_obj@misc[[".__h5ad_path__"]] <- normalizePath(file)
    seurat_obj@misc[[".__h5ad_loaded__"]] <- components
  }

  if (verbose) {
    message("\nSuccessfully loaded H5AD file")
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

#' Load deferred h5ad components into a Seurat object
#'
#' When \code{readH5AD} was called with a subset of \code{components},
#' this function loads the remaining components from the stored file path.
#'
#' @param object A Seurat object created by \code{readH5AD} with partial components
#' @param components Character vector of additional components to load, or NULL for all
#' @param verbose Show progress messages
#'
#' @return The Seurat object with additional components loaded
#'
#' @export
#'
scLoadMeta <- function(object, components = NULL, verbose = TRUE) {
  path <- object@misc[[".__h5ad_path__"]]
  if (is.null(path)) {
    stop("No h5ad source path stored -- object was not created with partial components",
         call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Source h5ad file no longer exists: ", path, call. = FALSE)
  }
  already <- object@misc[[".__h5ad_loaded__"]] %||% character(0)
  all_components <- c("X", "obs", "var", "obsm", "obsp", "varp", "layers", "uns")
  if (is.null(components)) {
    components <- setdiff(all_components, already)
  }
  if (length(components) == 0) return(object)

  # Re-read with the missing components using the full reader
  assay.name <- names(object@assays)[1]
  full <- readH5AD(path, assay.name = assay.name, components = all_components,
                   verbose = verbose)

  # Merge missing components
  if ("obs" %in% components && !"obs" %in% already) {
    meta_cols <- setdiff(names(full@meta.data), names(object@meta.data))
    if (length(meta_cols) > 0) {
      object <- AddMetaData(object, metadata = full@meta.data[, meta_cols, drop = FALSE])
    }
  }
  if ("obsm" %in% components && !"obsm" %in% already) {
    for (r in names(full@reductions)) {
      if (!r %in% names(object@reductions)) {
        object@reductions[[r]] <- full@reductions[[r]]
      }
    }
  }
  if ("obsp" %in% components && !"obsp" %in% already) {
    for (g in names(full@graphs)) {
      if (!g %in% names(object@graphs)) {
        object@graphs[[g]] <- full@graphs[[g]]
      }
    }
  }
  object@misc[[".__h5ad_loaded__"]] <- union(already, components)
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# C-accelerated h5ad reader (Optimization 1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Fast Seurat object constructor (bypass CreateSeuratObject validation)
#'
#' Constructs a Seurat object directly from a pre-validated sparse matrix,
#' skipping the expensive validation in CreateSeuratObject (cell/feature
#' filtering, name checks, etc.). Use only when data is known-good (e.g.,
#' freshly read from a validated HDF5 file).
#'
#' @param counts dgCMatrix with rownames (features) and colnames (cells)
#' @param assay.name Name for the assay
#' @return Seurat object
#' @keywords internal
#'
.fast_create_seurat <- function(counts, assay.name = "RNA") {
  CreateSeuratObject(counts = counts, project = "H5AD",
                     assay = assay.name, min.cells = 0, min.features = 0)
}

.readH5AD_c <- function(file, assay.name = "RNA", components = NULL, verbose = TRUE) {
  if (verbose) message("Loading H5AD file (C reader): ", file)

  # Call C reader for requested components
  result <- .Call(C_read_h5ad, file, components)
  if (is.null(result)) {
    if (verbose) message("C reader failed, falling back to R reader")
    return(readH5AD(file, assay.name = assay.name, components = components,
                    use.c = FALSE, verbose = verbose))
  }

  # 1. Construct expression matrix from C result
  mat_data <- result[["X"]]
  if (is.null(mat_data)) {
    stop("C reader returned no X matrix", call. = FALSE)
  }

  is_dense <- !is.null(attr(mat_data, "dense"))
  if (is_dense) {
    # Dense matrix returned as vector with dim attribute
    expr_matrix <- mat_data
  } else {
    # Sparse: construct dgCMatrix from components
    expr_matrix <- new("dgCMatrix",
      i = mat_data$i,
      p = mat_data$p,
      x = mat_data$x,
      Dim = c(mat_data$nrow, mat_data$ncol)
    )
  }

  # Set dimnames
  cell.names <- mat_data$colnames
  feature.names <- mat_data$rownames
  if (is.null(cell.names)) cell.names <- paste0("Cell", seq_len(ncol(expr_matrix)))
  if (is.null(feature.names)) feature.names <- paste0("Gene", seq_len(nrow(expr_matrix)))
  if (anyDuplicated(feature.names)) feature.names <- make.unique(feature.names)
  rownames(expr_matrix) <- feature.names
  colnames(expr_matrix) <- cell.names

  # 2. Create Seurat object (fast path bypasses validation overhead)
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- .fast_create_seurat(expr_matrix, assay.name = assay.name)

  # 3. Add obs metadata
  if ("obs" %in% components && !is.null(result[["obs"]])) {
    if (verbose) message("Adding cell metadata...")
    obs_data <- result[["obs"]]
    # Remove _index (already used as cell names)
    obs_data[["_index"]] <- NULL
    if (length(obs_data) > 0) {
      batch_df <- data.frame(row.names = colnames(seurat_obj))
      for (col in names(obs_data)) {
        batch_df[[col]] <- obs_data[[col]]
      }
      seurat_obj <- AddMetaData(seurat_obj, metadata = batch_df)
    }
  }

  # 4. Add var metadata
  if ("var" %in% components && !is.null(result[["var"]])) {
    if (verbose) message("Adding feature metadata...")
    var_data <- result[["var"]]
    var_data[["_index"]] <- NULL
    for (col in names(var_data)) {
      tryCatch({
        meta_values <- var_data[[col]]
        if (col == "highly_variable") {
          if (is.factor(meta_values)) meta_values <- as.character(meta_values) == "True"
          else if (is.numeric(meta_values)) meta_values <- as.logical(meta_values)
        }
        if (length(meta_values) == nrow(seurat_obj)) {
          names(meta_values) <- feature.names
          seurat_obj[[assay.name]][[col]] <- meta_values
          if (col == "highly_variable" && is.logical(meta_values)) {
            VariableFeatures(seurat_obj) <- feature.names[meta_values]
          }
        }
      }, error = function(e) {
        if (verbose) message("Could not add feature metadata ", col, ": ", e$message)
      })
    }
  }

  # 5. Add obsm reductions
  if ("obsm" %in% components && !is.null(result[["obsm"]])) {
    if (verbose) message("Adding dimensional reductions...")
    for (reduc_name in names(result[["obsm"]])) {
      clean_name <- gsub("^X_", "", reduc_name)
      if (clean_name == "spatial") next
      tryCatch({
        embeddings <- result[["obsm"]][[reduc_name]]
        if (!is.matrix(embeddings)) embeddings <- as.matrix(embeddings)
        if (nrow(embeddings) != length(cell.names) && ncol(embeddings) == length(cell.names)) {
          embeddings <- t(embeddings)
        }
        if (nrow(embeddings) != length(cell.names)) next
        rownames(embeddings) <- cell.names
        key <- AnnDataReductionKey(clean_name)
        colnames(embeddings) <- paste0(gsub("_$", "", key), "_", seq_len(ncol(embeddings)))
        seurat_obj[[clean_name]] <- CreateDimReducObject(
          embeddings = embeddings, key = key, assay = assay.name
        )
        if (verbose) message("  Added reduction: ", clean_name)
      }, error = function(e) {
        if (verbose) message("Could not add reduction ", clean_name, ": ", e$message)
      })
    }
  }

  # 6. Add obsp graphs
  if ("obsp" %in% components && !is.null(result[["obsp"]])) {
    if (verbose) message("Adding neighbor graphs...")
    for (graph_name in names(result[["obsp"]])) {
      tryCatch({
        gd <- result[["obsp"]][[graph_name]]
        graph_matrix <- new("dgCMatrix",
          i = gd$i, p = gd$p, x = gd$x,
          Dim = c(gd$nrow, gd$ncol)
        )
        # obsp CSR was reinterpreted as CSC = transpose, so transpose back
        graph_matrix <- Matrix::t(graph_matrix)
        if (nrow(graph_matrix) == ncol(graph_matrix) && nrow(graph_matrix) == ncol(seurat_obj)) {
          rownames(graph_matrix) <- cell.names
          colnames(graph_matrix) <- cell.names
          seurat_graph_name <- switch(graph_name,
            "connectivities" = paste0(assay.name, "_snn"),
            "distances" = paste0(assay.name, "_nn"),
            graph_name
          )
          graph_obj <- as.Graph(graph_matrix)
          DefaultAssay(graph_obj) <- assay.name
          seurat_obj@graphs[[seurat_graph_name]] <- graph_obj
          if (verbose) message("  Added graph: ", graph_name)
        }
      }, error = function(e) {
        if (verbose) message("Could not add graph ", graph_name, ": ", e$message)
      })
    }
  }

  # 7. Handle remaining components via hdf5r fallback (varp, layers, uns, spatial, raw)
  # These are less performance-critical, so use the existing R reader
  needs_hdf5r <- any(c("varp", "layers", "uns") %in% components)
  if (needs_hdf5r) {
    h5ad <- H5File$new(file, mode = "r")
    on.exit(h5ad$close_all())

    # raw/X handling: set counts and data layers
    if ("layers" %in% components || "X" %in% components) {
      if (h5ad$exists("raw") && h5ad[["raw"]]$exists("X")) {
        if (verbose) message("Adding raw counts...")
        raw_features <- NULL
        if (h5ad[["raw"]]$exists("var")) {
          raw_var <- h5ad[["raw/var"]]
          if (raw_var$exists("_index")) raw_features <- as.character(raw_var[["_index"]][])
          else if (raw_var$exists("index")) raw_features <- as.character(raw_var[["index"]][])
        }
        if (!is.null(raw_features)) {
          # Use the same ReadH5ADMatrix helper pattern
          raw_obj <- h5ad[["raw/X"]]
          if (inherits(raw_obj, "H5Group") && raw_obj$exists("data")) {
            data_vals <- raw_obj[["data"]][]
            indices <- raw_obj[["indices"]][]
            indptr <- raw_obj[["indptr"]][]
            encoding <- tryCatch(h5attr(raw_obj, "encoding-type"), error = function(e) "csr_matrix")
            if (identical(encoding, "csc_matrix")) {
              raw_matrix <- new("dgCMatrix",
                i = as.integer(indices), p = as.integer(indptr), x = as.numeric(data_vals),
                Dim = c(as.integer(length(raw_features)), as.integer(length(cell.names)))
              )
            } else {
              # CSR of (n_cells x n_features) == CSC of (n_features x n_cells)
              raw_matrix <- new("dgCMatrix",
                i = as.integer(indices), p = as.integer(indptr), x = as.numeric(data_vals),
                Dim = c(as.integer(length(raw_features)), as.integer(length(cell.names)))
              )
            }
          } else if (inherits(raw_obj, "H5D")) {
            raw_matrix <- t(raw_obj[,])
          }
          if (exists("raw_matrix")) {
            common <- intersect(feature.names, raw_features)
            if (length(common) > 0) {
              rownames(raw_matrix) <- raw_features
              colnames(raw_matrix) <- cell.names
              seurat_obj[[assay.name]] <- SetAssayData(
                seurat_obj[[assay.name]], layer = "counts", new.data = raw_matrix[common, , drop = FALSE]
              )
              seurat_obj[[assay.name]] <- SetAssayData(
                seurat_obj[[assay.name]], layer = "data", new.data = expr_matrix[common, , drop = FALSE]
              )
            }
          }
        }
      }
    }

    # Layers
    if ("layers" %in% components && h5ad$exists("layers")) {
      if (verbose) message("Adding layers...")
      for (layer_name in names(h5ad[["layers"]])) {
        tryCatch({
          layer_obj <- h5ad[["layers"]][[layer_name]]
          if (inherits(layer_obj, "H5Group") && layer_obj$exists("data")) {
            ld <- layer_obj[["data"]][]; li <- layer_obj[["indices"]][]; lp <- layer_obj[["indptr"]][]
            # CSR of (cells x features) == CSC of (features x cells)
            layer_matrix <- new("dgCMatrix", i = as.integer(li), p = as.integer(lp),
                                x = as.numeric(ld), Dim = c(as.integer(length(feature.names)),
                                                             as.integer(length(cell.names))))
          } else if (inherits(layer_obj, "H5D")) {
            layer_matrix <- t(layer_obj[,])
          } else next
          if (nrow(layer_matrix) == nrow(expr_matrix) && ncol(layer_matrix) == ncol(expr_matrix)) {
            rownames(layer_matrix) <- feature.names
            colnames(layer_matrix) <- cell.names
            seurat_slot <- switch(layer_name,
              "counts" = "counts", "data" = "data",
              "log_normalized" = "data", "scale.data" = "scale.data",
              "scaled" = "scale.data", layer_name)
            seurat_obj[[assay.name]] <- SetAssayData(seurat_obj[[assay.name]],
                                                     layer = seurat_slot, new.data = layer_matrix)
          }
        }, error = function(e) {
          if (verbose) message("Could not add layer ", layer_name, ": ", e$message)
        })
      }
    }

    # varp
    if ("varp" %in% components && h5ad$exists("varp")) {
      if (verbose) message("Adding varp...")
      varp_list <- list()
      for (vn in names(h5ad[["varp"]])) {
        tryCatch({
          vo <- h5ad[["varp"]][[vn]]
          if (inherits(vo, "H5Group") && vo$exists("data")) {
            vd <- vo[["data"]][]; vi <- vo[["indices"]][]; vp <- vo[["indptr"]][]
            vm <- new("dgCMatrix", i = as.integer(vi), p = as.integer(vp),
                      x = as.numeric(vd), Dim = rep(as.integer(length(feature.names)), 2))
          } else if (inherits(vo, "H5D")) {
            vm <- vo[,]
          }
          if (exists("vm")) { varp_list[[vn]] <- vm; rm(vm) }
        }, error = function(e) NULL)
      }
      if (length(varp_list) > 0) seurat_obj@misc[["__varp__"]] <- varp_list
    }

    # uns
    if ("uns" %in% components && h5ad$exists("uns")) {
      if (verbose) message("Adding unstructured data...")
      .read_uns_group <- function(grp) {
        result <- list()
        for (item in names(grp)) {
          tryCatch({
            if (inherits(grp[[item]], "H5D")) result[[item]] <- grp[[item]]$read()
            else if (inherits(grp[[item]], "H5Group")) result[[item]] <- .read_uns_group(grp[[item]])
          }, error = function(e) NULL)
        }
        result
      }
      for (item in names(h5ad[["uns"]])) {
        tryCatch({
          if (inherits(h5ad[["uns"]][[item]], "H5D")) seurat_obj@misc[[item]] <- h5ad[["uns"]][[item]]$read()
          else if (inherits(h5ad[["uns"]][[item]], "H5Group"))
            seurat_obj@misc[[item]] <- .read_uns_group(h5ad[["uns"]][[item]])
        }, error = function(e) NULL)
      }
    }

    # Spatial
    if ("obsm" %in% components && h5ad$exists("obsm") && "spatial" %in% names(h5ad[["obsm"]])) {
      seurat_obj <- H5ADSpatialToSeurat(h5ad_file = h5ad, seurat_obj = seurat_obj,
                                         assay_name = assay.name, verbose = verbose)
    }
  }

  # Store source path for deferred loading
  all_components <- c("X", "obs", "var", "obsm", "obsp", "varp", "layers", "uns")
  if (!setequal(components, all_components)) {
    seurat_obj@misc[[".__h5ad_path__"]] <- normalizePath(file)
    seurat_obj@misc[[".__h5ad_loaded__"]] <- components
  }

  if (verbose) {
    message("\nSuccessfully loaded H5AD file")
    message("  Cells: ", ncol(seurat_obj))
    message("  Features: ", nrow(seurat_obj))
  }

  return(seurat_obj)
}
