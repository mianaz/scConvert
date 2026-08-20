#' Construct a typed scConvert_data_error condition.
#' @keywords internal
#' @noRd
.scconvert_data_error <- function(msg, ...) {
  structure(
    class = c("scConvert_data_error", "error", "condition"),
    list(message = msg, call = NULL, ...)
  )
}

#' Check that a sparse h5ad group has internally consistent
#' `data`/`indices` lengths. Returns silently on a non-sparse group;
#' raises a scConvert_data_error condition on mismatch. Used to catch
#' truncated payloads up-front so the C-level sparse reader never sees
#' a partial buffer.
#' @keywords internal
#' @noRd
.h5ad_check_sparse_lengths <- function(grp, label) {
  if (!inherits(grp, "H5Group")) return(invisible(NULL))
  if (!(grp$exists("data") && grp$exists("indices") && grp$exists("indptr"))) {
    return(invisible(NULL))
  }
  n_data    <- grp[["data"]]$dims[1]
  n_indices <- grp[["indices"]]$dims[1]
  if (is.null(n_data) || is.null(n_indices)) return(invisible(NULL))
  if (n_data != n_indices) {
    stop(.scconvert_data_error(sprintf(
      "Malformed h5ad: /%s/data has %d entries but /%s/indices has %d",
      label, n_data, label, n_indices)))
  }
  invisible(TRUE)
}

#' Structural validation for h5ad inputs.
#'
#' Two-tier check:
#'   * Required skeleton: `/X` (or `/raw/X`), `/obs`, `/var` must exist.
#'     Missing any of these raises `scConvert_data_error` immediately;
#'     the reader cannot recover.
#'   * Sparse-payload consistency: every sparse group reachable under
#'     `/X`, `/raw/X`, `/layers/*`, `/obsp/*`, and `/varp/*` is checked
#'     for `len(data) == len(indices)`. A truncated payload here would
#'     otherwise segfault the C-level sparse constructor; raising the
#'     typed condition lets the harness catch and report it.
#'
#' Auxiliary slot integrity (per-component shape vs `n_cells` /
#' `n_features` for `obsm`/`varm`, etc.) is enforced by the readers at
#' load time with skip-and-warn semantics, since a single corrupt
#' embedding should not prevent the rest of the object from loading.
#'
#' @keywords internal
#' @noRd
.h5ad_validate_minimum_structure <- function(file) {
  h <- tryCatch(hdf5r::H5File$new(file, mode = "r"),
                error = function(e) {
                  stop(.scconvert_data_error(sprintf(
                    "Cannot open as HDF5: %s", conditionMessage(e))))
                })
  on.exit(tryCatch(h$close_all(), error = function(e) NULL), add = TRUE)
  has_X    <- h$exists("X") || (h$exists("raw") && h[["raw"]]$exists("X"))
  has_obs  <- h$exists("obs")
  has_var  <- h$exists("var")
  missing <- character()
  if (!has_X)   missing <- c(missing, "/X")
  if (!has_obs) missing <- c(missing, "/obs")
  if (!has_var) missing <- c(missing, "/var")
  if (length(missing) > 0L) {
    stop(.scconvert_data_error(sprintf(
      "Malformed h5ad: missing required group(s): %s",
      paste(missing, collapse = ", ")), missing = missing))
  }

  # Walk every sparse-encoded matrix and check data/indices length parity.
  if (h$exists("X")) .h5ad_check_sparse_lengths(h[["X"]], "X")
  if (h$exists("raw") && h[["raw"]]$exists("X")) {
    .h5ad_check_sparse_lengths(h[["raw"]][["X"]], "raw/X")
  }
  for (parent in c("layers", "obsm", "obsp", "varm", "varp")) {
    if (!h$exists(parent)) next
    pg <- h[[parent]]
    if (!inherits(pg, "H5Group")) next
    for (nm in names(pg)) {
      child <- tryCatch(pg[[nm]], error = function(e) NULL)
      if (is.null(child)) next
      .h5ad_check_sparse_lengths(child, sprintf("%s/%s", parent, nm))
    }
  }
  invisible(TRUE)
}

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
#' @param reductions Which \code{obsm} entries to load as dimensional
#'   reductions. \code{NULL} (default) loads all of them under cleaned names
#'   (leading \code{X_} stripped). Pass a character vector of obsm keys to
#'   load a subset, optionally named to control the Seurat reduction names:
#'   \code{reductions = c(scvi = "X_scVI", umap = "X_umap")} loads only those
#'   two keys as reductions \code{scvi} and \code{umap}. Name collisions
#'   (e.g. \code{X_pca} and \code{pca} both present) are resolved by keeping
#'   the raw obsm key for later claimants, with a warning, instead of the
#'   previous silent overwrite.
#'
#' @return A \code{Seurat} object. If \code{use.bpcells} is set, the count matrix
#'   is stored on disk in BPCells format and the object uses minimal memory.
#'
#' @section Post-read verification:
#' Before returning, the loaded object is verified against the file: dims
#' must match the file's obs/var counts (orientation check), cell names must
#' equal the obs index in order, feature names must equal the var index
#' modulo Seurat's documented underscore-to-dash replacement, and no
#' duplicate barcodes/features may survive. Violations raise
#' \code{scConvert_data_error}. Duplicate names in the file are made unique
#' with a \code{scConvert_names_warning}; a counts layer left holding
#' non-integer values raises a \code{scConvert_counts_warning}. What the
#' reader did (which slot became the counts layer, where X went, the
#' scConvert version) is recorded in \code{misc$scConvert_read}.
#'
#' @importFrom hdf5r H5File h5attr h5attr_names
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject SetAssayData CreateDimReducObject
#' @importFrom SeuratObject Cells
#'
#' @export
#'
readH5AD <- function(file, assay.name = "RNA", use.bpcells = NULL,
                     components = NULL, use.c = TRUE, verbose = TRUE,
                     reductions = NULL) {
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  # Structural validation: the AnnData spec requires /X, /obs, and /var.
  # We reject inputs that are missing any of the three rather than silently
  # substituting empty placeholders, which previously let malformed files
  # through and produced a Seurat object with auto-generated cell/feature
  # names. Errors raised here inherit from "scConvert_data_error" so callers
  # can distinguish a malformed-input refusal from a library-side bug.
  .h5ad_validate_minimum_structure(file)

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
                   reductions = reductions, verbose = verbose),
      error = function(e) {
        if (verbose) message("C reader failed (", conditionMessage(e), "), using R reader")
        NULL
      }
    )
    if (!is.null(c_result)) return(c_result)
  }

  h5ad <- H5File$new(file, mode = "r")
  # close_all() throws on hdf5r/HDF5 1.12.1 (CRAN Windows binary) when a
  # leaf-object ID isn't fully released. Read already succeeded by then.
  # Matches the write-path wrap in R/WriteH5AD.R (commit cf362e0).
  on.exit(tryCatch(h5ad$close_all(), error = function(e) NULL))

  if (verbose) {
    message("Loading H5AD file: ", file)
  }

  # Read h5ad sparse/dense matrices via the shared package-level reader
  # (.h5ad_read_matrix in R/VerifyH5AD.R), which this closure used to inline.
  ReadH5ADMatrix <- .h5ad_read_matrix

  # Resolve where the raw counts live before anything is loaded: the
  # /uns/scConvert stamp when present (version-aware branch), else
  # /layers/counts (modern convention), /raw/X (legacy scanpy), or /X.
  counts_source <- .h5ad_resolve_counts_source(h5ad)

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

  # Deduplicate names (some datasets have duplicates, e.g. squidpy four_i).
  # Duplicates are a data-integrity hazard (silent misalignment on any
  # name-based join), so the rename is loud (scConvert_names_warning) and
  # recorded in misc$scConvert_read. The original var index survives in the
  # 'orig_var_index' feature metadata column via .h5ad_preserve_var_identity.
  var_index_original <- feature.names
  dedup_features <- anyDuplicated(feature.names) > 0L
  if (dedup_features) {
    .h5ad_warn_duplicate_names("feature", sum(duplicated(feature.names)))
    feature.names <- make.unique(feature.names)
  }
  dedup_cells <- anyDuplicated(cell.names) > 0L
  if (dedup_cells) {
    .h5ad_warn_duplicate_names("cell", sum(duplicated(cell.names)))
    cell.names <- make.unique(cell.names)
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

    # Load the resolved counts source on disk -- matches the non-BPCells
    # path's convention. Previously only raw/X was considered, so files
    # following the modern layers['counts'] convention silently served
    # log-normalized X values as on-disk "counts".
    bp_group <- counts_source
    if (verbose) {
      if (identical(bp_group, "X")) {
        message("Loading expression matrix (X) via BPCells (on-disk)...")
      } else {
        message("Loading raw counts (", bp_group, ") via BPCells (on-disk)...")
      }
    }

    # Close hdf5r handle before BPCells opens the file (avoid lock conflicts)
    tryCatch(h5ad$close_all(), error = function(e) NULL)

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

  # 5. Apply the counts/data convention for the resolved counts source
  # (.h5ad_apply_counts_convention in R/VerifyH5AD.R):
  #   raw/X         -> counts <- raw/X (common features), data <- X
  #   layers/counts -> counts <- layers/counts,           data <- X
  #   X             -> counts stays X; data <- copy (call NormalizeData()
  #                    when real normalized values are needed)
  # In BPCells mode the chosen source is already the on-disk matrix, so the
  # in-memory relocation is skipped; when that source is not X itself, X was
  # never loaded at all.
  if (use_bpcells) {
    x_mapped_to <- if (identical(counts_source, "X")) "counts" else "not_loaded"
  } else {
    counts_conv <- .h5ad_apply_counts_convention(
      seurat_obj, h5ad, assay.name, expr_matrix, feature.names, cell.names,
      counts_source, verbose = verbose
    )
    seurat_obj <- counts_conv$object
    x_mapped_to <- counts_conv$x_mapped_to
  }

  # 6. Add layers if present
  if (h5ad$exists("layers")) {
    if (use_bpcells) {
      if (verbose) message("Skipping additional layers (BPCells on-disk mode)")
    } else {
      if (verbose) message("Adding layers...")
      layer_names <- names(h5ad[["layers"]])

      for (layer_name in layer_names) {
        # layers/counts was already placed in the counts layer by the
        # counts-convention step (5); re-reading it here would be wasted IO
        # and, worse, its failure mode used to be a swallowed message that
        # left log-normalized X sitting in the counts slot.
        if (identical(layer_name, "counts") &&
            identical(counts_source, "layers/counts")) {
          next
        }
        if (verbose) message("  Adding layer: ", layer_name)

        # hdf5r reads a dense h5py (cells x features) dataset as an R matrix
        # already shaped (features x cells), matching Seurat's expr_matrix
        # convention. ReadH5ADMatrix(..., transpose = TRUE) used to add a
        # second t() on top of that, leaving every dense layer in the wrong
        # orientation and silently failing the shape check that follows.
        # For dense layers, read at the natural hdf5r orientation; for
        # sparse, ReadH5ADMatrix already builds a (features x cells)
        # dgCMatrix.
        layer_obj <- h5ad[["layers"]][[layer_name]]
        layer_matrix <- if (inherits(layer_obj, "H5D")) {
          layer_obj[,]
        } else {
          ReadH5ADMatrix(layer_obj, transpose = TRUE)
        }

        # Layers must match X's (features x cells) shape. A mismatch means
        # the h5ad is malformed (ragged across slots). Previously we silently
        # skipped the layer with no diagnostic; that let bad files through
        # and made fidelity checks pass-by-omission. Refuse instead.
        if (nrow(layer_matrix) != nrow(expr_matrix) ||
            ncol(layer_matrix) != ncol(expr_matrix)) {
          cond <- structure(
            class = c("scConvert_data_error", "error", "condition"),
            list(message = sprintf(
                   "Malformed h5ad: layer '/layers/%s' has shape (%d x %d) but X has shape (%d x %d)",
                   layer_name,
                   nrow(layer_matrix), ncol(layer_matrix),
                   nrow(expr_matrix),  ncol(expr_matrix)),
                 call = NULL,
                 layer = layer_name))
          stop(cond)
        }
        # Label with the object's final dimnames (Seurat may have replaced
        # underscores with dashes in feature names; a file-index label would
        # make SetAssayData fail and silently drop the layer)
        rownames(layer_matrix) <- rownames(seurat_obj)
        colnames(layer_matrix) <- colnames(seurat_obj)

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

      # Read modern categoricals (groups with codes/categories). The stored
      # category order becomes the factor level order verbatim, and the
      # AnnData `ordered` flag round-trips into an ordered factor.
      for (col in cat_cols) {
        tryCatch({
          col_obj <- obs_group[[col]]
          encoding_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
          if (encoding_type == "categorical" && col_obj$exists("categories") && col_obj$exists("codes")) {
            codes <- col_obj[["codes"]]$read()
            categories <- as.character(col_obj[["categories"]]$read())
            is_ordered <- tryCatch(isTRUE(as.logical(h5attr(col_obj, "ordered"))[1]),
                                   error = function(e) FALSE)
            obs_batch[[col]] <- DecodeCategorical(codes, categories,
                                                  ordered = is_ordered)
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
      # Legacy compound dataset: obsm is a structured array with named fields.
      # exclude_spatial = FALSE: the modern spatial pipeline (step 12) only
      # sees H5Group obsm, so a legacy 'spatial' field must stay a reduction.
      obsm_df <- obsm_obj$read()
      obsm_plan <- .h5ad_plan_reductions(names(obsm_df),
                                         reductions = reductions,
                                         exclude_spatial = FALSE)
      n_cells <- length(cell.names)
      for (plan_i in seq_along(obsm_plan)) {
        reduc_name <- obsm_plan[[plan_i]]
        clean_name <- names(obsm_plan)[plan_i]
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
      # Modern format: obsm is an H5Group with named datasets.
      # 'spatial' is excluded from the plan -- handled separately in step 12.
      obsm_plan <- .h5ad_plan_reductions(names(obsm_obj),
                                         reductions = reductions,
                                         exclude_spatial = TRUE)
      for (plan_i in seq_along(obsm_plan)) {
        reduc_name <- obsm_plan[[plan_i]]
        clean_name <- names(obsm_plan)[plan_i]
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

          # Check dimensions. An embedding whose first dim doesn't match
          # n_cells is a malformed obsm entry; skip with a warning rather
          # than silently dropping it so the user is alerted.
          if (nrow(embeddings) != length(cell.names)) {
            warning(sprintf(
              "Skipping obsm/%s: embedding shape (%d x %d) does not match n_cells = %d",
              reduc_name, nrow(embeddings), ncol(embeddings),
              length(cell.names)), call. = FALSE)
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
            # Name by the object's final rownames, not the file var index:
            # Seurat's underscore-to-dash replacement would otherwise leave
            # the names unmatched and the assignment misaligned.
            names(meta_values) <- rownames(seurat_obj)
            seurat_obj[[assay.name]][[col]] <- meta_values
            if (col == "highly_variable" && is.logical(meta_values)) {
              VariableFeatures(seurat_obj) <- rownames(seurat_obj)[meta_values]
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
          # Modern h5ad categorical format; category order and the `ordered`
          # flag both round-trip into the factor.
          encoding_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
          if (encoding_type == "categorical" && col_obj$exists("categories") && col_obj$exists("codes")) {
            codes <- col_obj[["codes"]]$read()
            categories <- as.character(col_obj[["categories"]]$read())
            is_ordered <- tryCatch(isTRUE(as.logical(h5attr(col_obj, "ordered"))[1]),
                                   error = function(e) FALSE)
            meta_values <- DecodeCategorical(codes, categories,
                                             ordered = is_ordered)
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

        # Ensure length matches and name with the object's final rownames for
        # Seurat v5 compatibility (the file var index may differ after
        # Seurat's underscore-to-dash replacement).
        if (length(meta_values) == nrow(seurat_obj)) {
          names(meta_values) <- rownames(seurat_obj)
          seurat_obj[[assay.name]][[col]] <- meta_values

          # Set variable features if highly_variable column exists
          if (col == "highly_variable" && is.logical(meta_values)) {
            VariableFeatures(seurat_obj) <- rownames(seurat_obj)[meta_values]
          }
        }
      }, error = function(e) {
        if (verbose) message("Could not add feature metadata ", col, ": ", e$message)
      })
    }
    } # end else (non-compound var)
  }

  # 9b. Gene identity: if the final rownames differ from the file's var index
  # (make.unique dedup or Seurat's underscore mangling), keep the original
  # identifiers as feature metadata so identity is never silently lost.
  seurat_obj <- .h5ad_preserve_var_identity(seurat_obj, assay.name,
                                            var_index_original)

  # 10. Add neighbor graphs from obsp. Auxiliary slot: a single corrupt
  # graph should warn and skip rather than abort the whole load, matching
  # the varp pattern below.
  if (h5ad$exists("obsp")) {
    if (verbose) message("Adding neighbor graphs...")

    for (graph_name in names(h5ad[["obsp"]])) {
      if (verbose) message("  Adding graph: ", graph_name)

      tryCatch({
        # obsp stores CSR; read as-is then transpose to get correct orientation
        graph_matrix <- Matrix::t(ReadH5ADMatrix(h5ad[["obsp"]][[graph_name]], transpose = FALSE))

        if (nrow(graph_matrix) != ncol(graph_matrix) ||
            nrow(graph_matrix) != ncol(seurat_obj)) {
          warning(sprintf(
            "Skipping obsp/%s: shape (%d x %d) does not match n_cells = %d",
            graph_name, nrow(graph_matrix), ncol(graph_matrix), ncol(seurat_obj)),
            call. = FALSE)
          return(invisible(NULL))
        }
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
      }, error = function(e) {
        warning(sprintf("Skipping obsp/%s: %s", graph_name,
                        conditionMessage(e)), call. = FALSE)
      })
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

  # 12b. FOV rebuild: if any library under uns/spatial/ has segmentation or
  # molecules subgroups, reconstruct an FOV object and attach to the Seurat
  # object keyed by library name. Closes the CosMx/Xenium silent-loss gap.
  seurat_obj <- .rebuild_fovs_from_h5ad(h5ad, seurat_obj, verbose = verbose)

  # 13. Post-read verification (read-back discipline): assert shape /
  # orientation, cell/feature identity and order, and name uniqueness
  # against the file rather than trusting the load; surface a counts layer
  # that ended up holding non-integer values; record what the reader did.
  .h5ad_verify_read(seurat_obj, cell.names, feature.names, file)
  if (!use_bpcells) {
    .h5ad_warn_noninteger_counts(seurat_obj, assay.name, counts_source)
  }
  seurat_obj <- .h5ad_record_provenance(seurat_obj, file, counts_source,
                                        x_mapped_to,
                                        dedup_cells = dedup_cells,
                                        dedup_features = dedup_features)

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
  obj <- CreateSeuratObject(counts = counts, project = "H5AD",
                             assay = assay.name, min.cells = 0, min.features = 0)
  # Seurat 5's CreateSeuratObject(counts=) populates only the counts layer.
  # FeaturePlot / FetchData require "data" and error with "layer 'data' is
  # not found" when it is missing. Copy counts into data so downstream code
  # that assumes Seurat 4 semantics keeps working.
  data_layer <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay.name, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(data_layer) || length(data_layer) == 0L) {
    obj[[assay.name]] <- SetAssayData(
      object = obj[[assay.name]],
      layer = "data",
      new.data = counts
    )
  }
  obj
}

.readH5AD_c <- function(file, assay.name = "RNA", components = NULL,
                        reductions = NULL, verbose = TRUE) {
  if (verbose) message("Loading H5AD file (C reader): ", file)

  # Call C reader for requested components
  result <- .Call(C_read_h5ad, file, components)
  if (is.null(result)) {
    if (verbose) message("C reader failed, falling back to R reader")
    return(readH5AD(file, assay.name = assay.name, components = components,
                    use.c = FALSE, verbose = verbose, reductions = reductions))
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

  # Set dimnames; duplicates are made unique loudly (scConvert_names_warning)
  # and recorded in misc$scConvert_read, mirroring the R path.
  cell.names <- mat_data$colnames
  feature.names <- mat_data$rownames
  if (is.null(cell.names)) cell.names <- paste0("Cell", seq_len(ncol(expr_matrix)))
  if (is.null(feature.names)) feature.names <- paste0("Gene", seq_len(nrow(expr_matrix)))
  var_index_original <- feature.names
  dedup_features <- anyDuplicated(feature.names) > 0L
  if (dedup_features) {
    .h5ad_warn_duplicate_names("feature", sum(duplicated(feature.names)))
    feature.names <- make.unique(feature.names)
  }
  dedup_cells <- anyDuplicated(cell.names) > 0L
  if (dedup_cells) {
    .h5ad_warn_duplicate_names("cell", sum(duplicated(cell.names)))
    cell.names <- make.unique(cell.names)
  }
  rownames(expr_matrix) <- feature.names
  colnames(expr_matrix) <- cell.names

  # 2. Create Seurat object (fast path bypasses validation overhead)
  if (verbose) message("Creating Seurat object...")
  seurat_obj <- .fast_create_seurat(expr_matrix, assay.name = assay.name)

  # 2b. Counts/data convention. The C reader hands back /X only; where the
  # raw counts actually live (layers/counts, raw/X, or X itself) is resolved
  # against the file and the matrices are relocated accordingly. Runs
  # unconditionally -- counts identity is part of the primary-matrix
  # contract, not an optional component.
  h5ad <- H5File$new(file, mode = "r")
  on.exit(tryCatch(h5ad$close_all(), error = function(e) NULL), add = TRUE)
  counts_source <- .h5ad_resolve_counts_source(h5ad)
  counts_conv <- .h5ad_apply_counts_convention(
    seurat_obj, h5ad, assay.name, expr_matrix, feature.names, cell.names,
    counts_source, verbose = verbose
  )
  seurat_obj <- counts_conv$object
  x_mapped_to <- counts_conv$x_mapped_to

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
          names(meta_values) <- rownames(seurat_obj)
          seurat_obj[[assay.name]][[col]] <- meta_values
          if (col == "highly_variable" && is.logical(meta_values)) {
            VariableFeatures(seurat_obj) <- rownames(seurat_obj)[meta_values]
          }
        }
      }, error = function(e) {
        if (verbose) message("Could not add feature metadata ", col, ": ", e$message)
      })
    }
  }

  # 4b. The compiled reader preserves category order but drops the AnnData
  # `ordered` flag; restore it from the file's attributes. Also keep the
  # original var index when feature names were deduplicated or mangled.
  seurat_obj <- .h5ad_restore_ordered_factors(seurat_obj, h5ad, assay.name)
  seurat_obj <- .h5ad_preserve_var_identity(seurat_obj, assay.name,
                                            var_index_original)

  # 5. Add obsm reductions ('spatial' excluded from the plan: the spatial
  # pipeline in the hdf5r fallback below handles it)
  if ("obsm" %in% components && !is.null(result[["obsm"]])) {
    if (verbose) message("Adding dimensional reductions...")
    obsm_plan <- .h5ad_plan_reductions(names(result[["obsm"]]),
                                       reductions = reductions,
                                       exclude_spatial = TRUE)
    for (plan_i in seq_along(obsm_plan)) {
      reduc_name <- obsm_plan[[plan_i]]
      clean_name <- names(obsm_plan)[plan_i]
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

  # 7. Handle remaining components via the already-open hdf5r handle.
  # raw/X and layers/counts were consumed by the counts-convention step
  # (2b), which replaced the hand-rolled raw/X reconstruction that used to
  # live here (and, unlike it, validates encoding, sorts indices, and covers
  # the layers/counts convention).
  needs_hdf5r <- any(c("varp", "layers", "uns") %in% components)
  if (needs_hdf5r) {
    # Layers. Mirror the R-reader behaviour: a shape mismatch between a
    # layer and X is structurally malformed; raise scConvert_data_error
    # instead of silently dropping the layer.
    if ("layers" %in% components && h5ad$exists("layers")) {
      if (verbose) message("Adding layers...")
      for (layer_name in names(h5ad[["layers"]])) {
        # Already placed into the counts layer by the counts-convention step
        if (identical(layer_name, "counts") &&
            identical(counts_source, "layers/counts")) {
          next
        }
        layer_obj <- h5ad[["layers"]][[layer_name]]
        if (inherits(layer_obj, "H5Group") && layer_obj$exists("data")) {
          ld <- layer_obj[["data"]][]; li <- layer_obj[["indices"]][]; lp <- layer_obj[["indptr"]][]
          layer_matrix <- new("dgCMatrix", i = as.integer(li), p = as.integer(lp),
                              x = as.numeric(ld), Dim = c(as.integer(length(feature.names)),
                                                           as.integer(length(cell.names))))
        } else if (inherits(layer_obj, "H5D")) {
          layer_matrix <- t(layer_obj[,])
        } else next
        if (nrow(layer_matrix) != nrow(expr_matrix) ||
            ncol(layer_matrix) != ncol(expr_matrix)) {
          cond <- structure(
            class = c("scConvert_data_error", "error", "condition"),
            list(message = sprintf(
                   "Malformed h5ad: layer '/layers/%s' has shape (%d x %d) but X has shape (%d x %d)",
                   layer_name,
                   nrow(layer_matrix), ncol(layer_matrix),
                   nrow(expr_matrix),  ncol(expr_matrix)),
                 call = NULL, layer = layer_name))
          stop(cond)
        }
        tryCatch({
          rownames(layer_matrix) <- rownames(seurat_obj)
          colnames(layer_matrix) <- colnames(seurat_obj)
          seurat_slot <- switch(layer_name,
            "counts" = "counts", "data" = "data",
            "log_normalized" = "data", "scale.data" = "scale.data",
            "scaled" = "scale.data", layer_name)
          seurat_obj[[assay.name]] <- SetAssayData(seurat_obj[[assay.name]],
                                                   layer = seurat_slot, new.data = layer_matrix)
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

    # FOV rebuild (mirrors the non-components path above)
    seurat_obj <- .rebuild_fovs_from_h5ad(h5ad, seurat_obj, verbose = verbose)
  }

  # Spatial: gated on obsm alone. The handle is open regardless of
  # needs_hdf5r now, so a components subset like c("X", "obs", "obsm") no
  # longer silently drops the spatial image data.
  if ("obsm" %in% components && h5ad$exists("obsm") &&
      "spatial" %in% names(h5ad[["obsm"]])) {
    seurat_obj <- H5ADSpatialToSeurat(h5ad_file = h5ad, seurat_obj = seurat_obj,
                                       assay_name = assay.name, verbose = verbose)
  }

  # Post-read verification + provenance (mirrors the R reader)
  .h5ad_verify_read(seurat_obj, cell.names, feature.names, file)
  .h5ad_warn_noninteger_counts(seurat_obj, assay.name, counts_source)
  seurat_obj <- .h5ad_record_provenance(seurat_obj, file, counts_source,
                                        x_mapped_to,
                                        dedup_cells = dedup_cells,
                                        dedup_features = dedup_features)

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

# Walk uns/spatial/{library}/ groups and rebuild FOV objects for any library
# that has segmentation/ or molecules/ subgroups (written by WriteFOVToH5AD).
# Attaches each rebuilt FOV to the Seurat object keyed by library name,
# matching the library keys emitted by the writer.
.rebuild_fovs_from_h5ad <- function(h5ad, seurat_obj, verbose = TRUE) {
  if (!h5ad$exists("uns")) return(seurat_obj)
  uns_grp <- h5ad[["uns"]]
  if (!uns_grp$exists("spatial")) return(seurat_obj)
  sp_grp <- uns_grp[["spatial"]]
  lib_names <- names(sp_grp)
  if (length(lib_names) == 0) return(seurat_obj)

  # Use the Seurat object's default assay so the attached FOV binds to the
  # existing assay and Seurat does not warn about orphaned image data.
  default_assay <- tryCatch(Seurat::DefaultAssay(seurat_obj),
                             error = function(e) "Spatial")

  for (lib_id in lib_names) {
    lib_grp <- tryCatch(sp_grp[[lib_id]], error = function(e) NULL)
    if (is.null(lib_grp)) next
    if (!inherits(lib_grp, "H5Group")) next
    has_seg <- tryCatch(lib_grp$exists("segmentation"), error = function(e) FALSE)
    has_mol <- tryCatch(lib_grp$exists("molecules"),    error = function(e) FALSE)
    if (!isTRUE(has_seg) && !isTRUE(has_mol)) next

    # If H5ADSpatialToSeurat already attached a Visium-family image at this
    # library key, leave it alone. SeuratSpatialToH5AD emits boundary
    # centroids under segmentation/ as a Visium roundtrip byproduct, so the
    # seg/mol probe above fires on pure Visium inputs — overwriting with a
    # plain FOV would discard the scale.factors needed by SpatialFeaturePlot.
    existing_img <- tryCatch(seurat_obj[[lib_id]], error = function(e) NULL)
    if (inherits(existing_img, c("VisiumV1", "VisiumV2"))) {
      if (verbose) {
        message("  Keeping ", class(existing_img)[1], " image for '", lib_id,
                "'; skipping FOV rebuild")
      }
      next
    }

    fov <- tryCatch(
      ReadFOVFromH5AD(lib_grp, key = lib_id, assay = default_assay),
      error = function(e) {
        if (verbose) {
          message("  FOV rebuild failed for '", lib_id, "': ", e$message)
        }
        NULL
      }
    )
    if (is.null(fov)) next

    suppressWarnings(seurat_obj[[lib_id]] <- fov)
    if (verbose) message("  Rebuilt FOV '", lib_id, "' from uns/spatial")
  }
  seurat_obj
}
