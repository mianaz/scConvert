#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reverse-conversion (h5ad -> Seurat) integrity helpers
#
# readH5AD() historically trusted its own bookkeeping: whatever /X held was
# loaded into the Seurat counts layer, obsm keys were cleaned into reduction
# names with silent last-writer-wins collisions, and callers had no way to
# confirm that the returned object matched the file. The helpers here close
# those gaps:
#
#   * .h5ad_resolve_counts_source()   version/layout-aware detection of where
#                                     the raw counts live (uns stamp,
#                                     /layers/counts, /raw/X, /X)
#   * .h5ad_apply_counts_convention() single place that maps the resolved
#                                     source onto Seurat counts/data layers
#   * .h5ad_warn_noninteger_counts()  loud, classed warning when the counts
#                                     layer ends up holding normalized values
#   * .h5ad_plan_reductions()         obsm -> reduction-name plan with
#                                     collision handling and user remapping
#   * .h5ad_verify_read()             post-read shape/orientation and
#                                     name-identity assertions
#   * .h5ad_record_provenance()       what-the-reader-did record in misc
#   * .h5ad_stamp_provenance()        writer-side /uns/scConvert stamp that
#                                     future reads can branch on
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read an h5ad-encoded matrix (sparse CSR/CSC group or dense dataset)
#'
#' Package-level version of the reader previously defined as a closure inside
#' \code{readH5AD}, so the counts-convention helper and the C-path fallback can
#' share one implementation. \code{transpose = TRUE} returns the Seurat
#' orientation (features x cells) for an h5ad (cells x genes) matrix.
#'
#' @param h5_obj hdf5r H5Group (sparse) or H5D (dense)
#' @param transpose Transpose into features x cells orientation
#' @return dgCMatrix or dense matrix
#' @keywords internal
#' @noRd
.h5ad_read_matrix <- function(h5_obj, transpose = TRUE) {
  if (inherits(h5_obj, "H5Group")) {
    # Sparse matrix (CSR or CSC format in h5ad)
    if (h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
      # Pre-read structural sanity check: indptr[length(indptr)] is the
      # number of stored non-zeros, which must equal both length(data) and
      # length(indices). A truncated /X/data would otherwise read in
      # quietly and the downstream sparseMatrix call may segfault inside
      # the C-level constructor. Compare HDF5-reported dataset sizes
      # before pulling data into memory.
      data_d    <- h5_obj[["data"]]
      indices_d <- h5_obj[["indices"]]
      indptr_d  <- h5_obj[["indptr"]]
      n_data_h5    <- if (!is.null(data_d$dims))    data_d$dims[1]    else NA_integer_
      n_indices_h5 <- if (!is.null(indices_d$dims)) indices_d$dims[1] else NA_integer_
      n_indptr_h5  <- if (!is.null(indptr_d$dims))  indptr_d$dims[1]  else NA_integer_
      if (!is.na(n_data_h5) && !is.na(n_indices_h5) &&
          n_data_h5 != n_indices_h5) {
        cond <- structure(
          class = c("scConvert_data_error", "error", "condition"),
          list(message = sprintf(
                 "Malformed h5ad sparse group: data has %d entries but indices has %d",
                 n_data_h5, n_indices_h5),
               call = NULL))
        stop(cond)
      }
      # We can also cross-check against indptr[-1] once indptr is read.
      data_vals <- h5_obj[["data"]][]
      indices <- h5_obj[["indices"]][]   # 0-based
      indptr <- h5_obj[["indptr"]][]
      if (length(indptr) >= 1L) {
        declared_nnz <- as.integer(indptr[length(indptr)])
        if (declared_nnz != length(data_vals) ||
            declared_nnz != length(indices)) {
          cond <- structure(
            class = c("scConvert_data_error", "error", "condition"),
            list(message = sprintf(
                   "Malformed h5ad sparse group: indptr declares %d nonzeros but data has %d and indices has %d",
                   declared_nnz, length(data_vals), length(indices)),
                 call = NULL))
          stop(cond)
        }
      }

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

# Slots an h5ad counts stamp may legally point at.
.h5ad_counts_locations <- c("layers/counts", "raw/X", "X")

#' Does a nested h5 path ("layers/counts") exist in an open handle?
#' @keywords internal
#' @noRd
.h5ad_slot_exists <- function(h5ad, path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  node <- h5ad
  for (p in parts) {
    ok <- tryCatch(node$exists(p), error = function(e) FALSE)
    if (!isTRUE(ok)) return(FALSE)
    node <- tryCatch(node[[p]], error = function(e) NULL)
    if (is.null(node)) return(FALSE)
  }
  TRUE
}

#' Decide which slot of an h5ad file holds the raw counts
#'
#' Version/layout-aware resolution, in priority order:
#' \enumerate{
#'   \item \code{/uns/scConvert/counts_location} -- stamp written by
#'     scConvert's own h5ad writers recording where counts were placed.
#'     Files written by older scConvert (or other tools) carry no stamp and
#'     fall through to the structural heuristics.
#'   \item \code{/layers/counts} -- modern anndata/scanpy convention
#'     (X = normalized data).
#'   \item \code{/raw/X} -- legacy scanpy convention.
#'   \item \code{/X} -- no dedicated counts slot; X is all there is.
#' }
#'
#' @param h5ad Open hdf5r H5File handle
#' @return character(1): "layers/counts", "raw/X", or "X"
#' @keywords internal
#' @noRd
.h5ad_resolve_counts_source <- function(h5ad) {
  stamped <- tryCatch({
    if (.h5ad_slot_exists(h5ad, "uns/scConvert/counts_location")) {
      as.character(h5ad[["uns"]][["scConvert"]][["counts_location"]]$read())[1]
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (!is.null(stamped) && stamped %in% .h5ad_counts_locations &&
      .h5ad_slot_exists(h5ad, stamped)) {
    return(stamped)
  }
  # A stamp of "none" (data-only export) still resolves to X below: the
  # non-integer counts warning is what surfaces the missing-counts state.
  if (.h5ad_slot_exists(h5ad, "layers/counts")) return("layers/counts")
  if (.h5ad_slot_exists(h5ad, "raw/X")) return("raw/X")
  "X"
}

#' Move X / raw/X / layers/counts into the Seurat layers the convention implies
#'
#' \code{CreateSeuratObject()} is always fed /X, whatever /X holds. This step
#' relocates matrices so the Seurat layers match the resolved counts source:
#' \itemize{
#'   \item \code{"raw/X"}: counts <- raw/X (common features), data <- X
#'   \item \code{"layers/counts"}: counts <- layers/counts, data <- X
#'   \item \code{"X"}: counts stays X; data <- copy of counts
#' }
#' Runs on every load, not gated on \code{components}: counts identity is part
#' of the primary-matrix contract, exactly like the raw/X handling always was.
#' Failure to relocate X out of the counts layer is the silent
#' lognorm-in-counts failure this step exists to prevent, so the
#' layers/counts branch raises \code{scConvert_data_error} instead of
#' downgrading to a message.
#'
#' @param seurat_obj Seurat object freshly built from /X
#' @param h5ad Open hdf5r H5File handle
#' @param assay.name Assay to populate
#' @param expr_matrix The matrix read from /X (features x cells, dimnames set)
#' @param feature.names Final (deduplicated) feature names
#' @param cell.names Final (deduplicated) cell names
#' @param counts_source Result of \code{.h5ad_resolve_counts_source()}
#' @param verbose Emit progress messages
#' @return list(object = <Seurat>, x_mapped_to = "counts"|"data")
#' @keywords internal
#' @noRd
.h5ad_apply_counts_convention <- function(seurat_obj, h5ad, assay.name,
                                          expr_matrix, feature.names, cell.names,
                                          counts_source, verbose = TRUE) {
  x_mapped_to <- "counts"

  # Label relocation matrices with the object's final dimnames, not the file
  # index: Seurat's underscore-to-dash replacement would otherwise make
  # SetAssayData fail with "no feature overlap" (or, in the pre-fix layers
  # loop, silently skip -- leaving log-normalized X in the counts slot for
  # any file with underscore gene names). Positional relabeling is safe:
  # the row/col order is the file order, which .h5ad_verify_read asserts.
  object_features <- rownames(seurat_obj)
  object_cells <- colnames(seurat_obj)
  x_relabeled <- expr_matrix
  if (nrow(x_relabeled) == length(object_features) &&
      ncol(x_relabeled) == length(object_cells)) {
    rownames(x_relabeled) <- object_features
    colnames(x_relabeled) <- object_cells
  }

  if (identical(counts_source, "raw/X")) {
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
      raw_matrix <- .h5ad_read_matrix(h5ad[["raw/X"]], transpose = TRUE)

      # Handle dimension mismatches (dense matrices may need additional transpose)
      n_raw_features <- length(raw_features)
      n_cells <- length(cell.names)
      if (nrow(raw_matrix) == n_cells && ncol(raw_matrix) == n_raw_features &&
          nrow(raw_matrix) != n_raw_features) {
        raw_matrix <- t(raw_matrix)
      }

      # Match dimensions
      raw_features <- raw_features[seq_len(min(length(raw_features), nrow(raw_matrix)))]
      # Same Seurat name transformation the object went through, so the
      # intersection happens in the object's namespace.
      rownames(raw_matrix) <- gsub("_", "-", raw_features, fixed = TRUE)
      colnames(raw_matrix) <- object_cells

      # Find common features
      common_features <- intersect(object_features, rownames(raw_matrix))
      if (length(common_features) > 0) {
        # raw/X -> counts layer (actual raw counts)
        seurat_obj[[assay.name]] <- SetAssayData(
          object = seurat_obj[[assay.name]],
          layer = "counts",
          new.data = raw_matrix[common_features, , drop = FALSE]
        )
        # X (loaded as counts at creation) -> data layer (normalized)
        seurat_obj[[assay.name]] <- SetAssayData(
          object = seurat_obj[[assay.name]],
          layer = "data",
          new.data = x_relabeled[common_features, , drop = FALSE]
        )
        x_mapped_to <- "data"
        if (verbose) message("  Set raw/X as counts, X as data (normalized)")
      }
    }
  } else if (identical(counts_source, "layers/counts")) {
    if (verbose) message("Adding counts from layers/counts...")

    lyr_obj <- h5ad[["layers"]][["counts"]]
    # Dense layers arrive from hdf5r already shaped (features x cells);
    # sparse groups need the transpose into Seurat orientation. Mirrors the
    # generic layers loop in readH5AD().
    counts_matrix <- if (inherits(lyr_obj, "H5D")) {
      lyr_obj[, ]
    } else {
      .h5ad_read_matrix(lyr_obj, transpose = TRUE)
    }

    if (nrow(counts_matrix) != nrow(expr_matrix) ||
        ncol(counts_matrix) != ncol(expr_matrix)) {
      stop(.scconvert_data_error(sprintf(
        "Malformed h5ad: layer '/layers/counts' has shape (%d x %d) but X has shape (%d x %d)",
        nrow(counts_matrix), ncol(counts_matrix),
        nrow(expr_matrix), ncol(expr_matrix)),
        layer = "counts"))
    }
    rownames(counts_matrix) <- object_features
    colnames(counts_matrix) <- object_cells

    # X first moves into data (it currently sits in the counts layer), then
    # the real counts replace it.
    seurat_obj[[assay.name]] <- SetAssayData(
      object = seurat_obj[[assay.name]],
      layer = "data",
      new.data = x_relabeled
    )
    seurat_obj[[assay.name]] <- SetAssayData(
      object = seurat_obj[[assay.name]],
      layer = "counts",
      new.data = counts_matrix
    )
    x_mapped_to <- "data"
    if (verbose) message("  Set layers/counts as counts, X as data (normalized)")
  }

  # Safety net for every source: guarantee a "data" layer exists. Seurat 5's
  # CreateSeuratObject(counts=) populates only the counts layer; FeaturePlot /
  # FetchData require "data" and fail with "layer 'data' is not found" without
  # it. In the counts_source == "X" case this is the documented copy (the user
  # is expected to call NormalizeData() when they need normalized values).
  tryCatch({
    # suppressWarnings: probing a missing data layer makes Seurat 5 warn
    # "Layer 'data' is empty" -- emptiness is exactly what is being tested.
    data_layer <- tryCatch(
      suppressWarnings(
        Seurat::GetAssayData(seurat_obj, assay = assay.name, layer = "data")
      ),
      error = function(e) NULL
    )
    if (is.null(data_layer) || length(data_layer) == 0L) {
      counts_layer <- Seurat::GetAssayData(seurat_obj, assay = assay.name,
                                            layer = "counts")
      seurat_obj[[assay.name]] <- SetAssayData(
        object = seurat_obj[[assay.name]],
        layer = "data",
        new.data = counts_layer
      )
    }
  }, error = function(e) {
    if (verbose) message("  Could not copy counts -> data layer: ", e$message)
  })

  list(object = seurat_obj, x_mapped_to = x_mapped_to)
}

#' Warn (classed) when the counts layer holds non-integer values
#'
#' The silent failure mode behind counts-layer ambiguity: a normalized /X
#' read into the Seurat counts slot breaks every downstream step that assumes
#' integer counts, without any signal. This check makes it loud. The warning
#' carries class \code{scConvert_counts_warning} so pipelines can
#' \code{withCallingHandlers()} on it specifically.
#'
#' @param object Seurat object after counts-convention application
#' @param assay.name Assay to check
#' @param counts_source Where the counts layer came from
#' @param max_check Cap on the number of values inspected
#' @return invisible(TRUE) if a warning was raised, invisible(FALSE) otherwise
#' @keywords internal
#' @noRd
.h5ad_warn_noninteger_counts <- function(object, assay.name, counts_source,
                                         max_check = 1e6L) {
  vals <- tryCatch({
    m <- Seurat::GetAssayData(object, assay = assay.name, layer = "counts")
    if (inherits(m, "dgCMatrix")) m@x else as.numeric(m)
  }, error = function(e) NULL)
  if (is.null(vals) || length(vals) == 0L) return(invisible(FALSE))
  if (length(vals) > max_check) vals <- vals[seq_len(max_check)]
  fractional <- vals != floor(vals)
  if (!isTRUE(any(fractional, na.rm = TRUE))) return(invisible(FALSE))

  msg <- if (identical(counts_source, "X")) {
    paste0(
      "readH5AD: the counts layer holds non-integer values taken from /X. ",
      "The file has no dedicated counts slot (/layers/counts or /raw/X), so ",
      "/X was used as counts even though it appears to contain normalized ",
      "data. Downstream steps that assume integer counts should not trust ",
      "this layer.")
  } else {
    sprintf(paste0(
      "readH5AD: the counts layer (read from /%s) contains non-integer ",
      "values; that slot may hold normalized data rather than raw counts."),
      counts_source)
  }
  warning(warningCondition(msg, class = "scConvert_counts_warning"))
  invisible(TRUE)
}

#' Emit a classed warning for name deduplication
#' @keywords internal
#' @noRd
.h5ad_warn_duplicate_names <- function(what, n_dup) {
  extra <- if (identical(what, "feature")) {
    " The original index is preserved in the 'orig_var_index' feature metadata column."
  } else {
    ""
  }
  warning(warningCondition(sprintf(
    "readH5AD: %d duplicate %s name(s) in the file were made unique with make.unique().%s",
    n_dup, what, extra), class = "scConvert_names_warning"))
}

#' Plan which obsm keys become which Seurat reductions
#'
#' Replaces the old inline \code{gsub("^X_", "", key)} mapping, which silently
#' overwrote reductions when two obsm keys cleaned to the same name (e.g.
#' \code{X_pca} and \code{pca}). Also implements the user-facing remapping
#' interface: \code{reductions = c(scvi = "X_scVI")} loads only \code{X_scVI}
#' and names it \code{scvi} in the Seurat object.
#'
#' @param obsm_keys Character vector of keys present under /obsm
#' @param reductions NULL (load all), or a character vector of obsm keys to
#'   load; names, when present, set the Seurat reduction names
#' @param exclude_spatial Drop spatial keys from the default plan (they are
#'   handled by the spatial/image pipeline instead)
#' @return Named character vector: values = obsm keys, names = reduction names
#' @keywords internal
#' @noRd
.h5ad_plan_reductions <- function(obsm_keys, reductions = NULL,
                                  exclude_spatial = TRUE) {
  if (is.null(reductions)) {
    keys <- obsm_keys
    if (exclude_spatial) {
      keys <- keys[gsub("^X_", "", keys) != "spatial"]
    }
    targets <- gsub("^X_", "", keys)
  } else {
    if (!is.character(reductions) || length(reductions) == 0L) {
      stop("`reductions` must be NULL or a character vector of obsm keys ",
           "(optionally named with the Seurat reduction names to use)",
           call. = FALSE)
    }
    missing_keys <- setdiff(reductions, obsm_keys)
    if (length(missing_keys) > 0L) {
      warning(warningCondition(sprintf(
        "readH5AD: requested obsm key(s) not present in file: %s. Available: %s",
        paste(missing_keys, collapse = ", "),
        paste(obsm_keys, collapse = ", ")),
        class = "scConvert_reduction_warning"))
    }
    keep <- reductions %in% obsm_keys
    keys <- unname(reductions[keep])
    targets <- if (is.null(names(reductions))) {
      rep("", length(keys))
    } else {
      names(reductions)[keep]
    }
    targets[!nzchar(targets)] <- gsub("^X_", "", keys[!nzchar(targets)])
  }
  if (length(keys) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }

  dup <- duplicated(targets)
  if (any(dup)) {
    collided <- unique(targets[dup])
    # First claimant keeps the clean name; later ones fall back to their raw
    # obsm key, then make.unique as a last resort.
    targets[dup] <- keys[dup]
    if (anyDuplicated(targets)) targets <- make.unique(targets)
    warning(warningCondition(sprintf(
      "readH5AD: multiple obsm keys map to the same reduction name (%s); kept all as: %s",
      paste(collided, collapse = ", "),
      paste(sprintf("%s -> %s", keys, targets), collapse = ", ")),
      class = "scConvert_reduction_warning"))
  }
  stats::setNames(keys, targets)
}

#' Restore the `ordered` flag on categorical metadata (C-reader path)
#'
#' The compiled reader builds plain factors (category order preserved, ordered
#' flag dropped). Re-read just the cheap `ordered` attributes via hdf5r and
#' upgrade the matching columns.
#'
#' @param seurat_obj Seurat object with metadata already attached
#' @param h5ad Open hdf5r H5File handle
#' @param assay.name Assay carrying feature metadata
#' @return The Seurat object with ordered factors restored
#' @keywords internal
#' @noRd
.h5ad_restore_ordered_factors <- function(seurat_obj, h5ad, assay.name) {
  ordered_cols <- function(grp_name) {
    out <- character(0)
    if (!isTRUE(tryCatch(h5ad$exists(grp_name), error = function(e) FALSE))) {
      return(out)
    }
    grp <- h5ad[[grp_name]]
    if (!inherits(grp, "H5Group")) return(out)
    for (nm in names(grp)) {
      child <- tryCatch(grp[[nm]], error = function(e) NULL)
      if (is.null(child) || !inherits(child, "H5Group")) next
      enc <- tryCatch(h5attr(child, "encoding-type"), error = function(e) "")
      if (!identical(enc, "categorical")) next
      ord <- tryCatch(isTRUE(as.logical(h5attr(child, "ordered"))[1]),
                      error = function(e) FALSE)
      if (ord) out <- c(out, nm)
    }
    out
  }

  for (col in ordered_cols("obs")) {
    v <- seurat_obj@meta.data[[col]]
    if (is.factor(v) && !is.ordered(v)) {
      seurat_obj@meta.data[[col]] <- factor(v, levels = levels(v), ordered = TRUE)
    }
  }
  for (col in ordered_cols("var")) {
    tryCatch({
      fmeta <- seurat_obj[[assay.name]][[]]
      v <- fmeta[[col]]
      if (is.factor(v) && !is.ordered(v)) {
        v <- factor(v, levels = levels(v), ordered = TRUE)
        names(v) <- rownames(seurat_obj[[assay.name]])
        seurat_obj[[assay.name]][[col]] <- v
      }
    }, error = function(e) NULL)
  }
  seurat_obj
}

#' Preserve the file's original var index when feature names had to change
#'
#' Feature names can drift from the file's var index through
#' \code{make.unique()} deduplication or Seurat's documented underscore-to-dash
#' replacement. When that happens the original identifiers (often the only
#' stable gene identity, e.g. Ensembl IDs) are kept as feature-level metadata
#' column \code{orig_var_index}, so gene identity is never silently lost and
#' round-trips can re-join on it.
#'
#' @param seurat_obj Seurat object with final feature names
#' @param assay.name Assay to annotate
#' @param var_index_original var index exactly as stored in the file
#' @return The Seurat object, possibly with an orig_var_index column added
#' @keywords internal
#' @noRd
.h5ad_preserve_var_identity <- function(seurat_obj, assay.name,
                                        var_index_original) {
  final_names <- rownames(seurat_obj)
  if (length(var_index_original) != length(final_names)) return(seurat_obj)
  if (identical(final_names, var_index_original)) return(seurat_obj)
  existing <- tryCatch(colnames(seurat_obj[[assay.name]][[]]),
                       error = function(e) character(0))
  if ("orig_var_index" %in% existing) return(seurat_obj)
  tryCatch({
    vals <- var_index_original
    names(vals) <- final_names
    seurat_obj[[assay.name]][["orig_var_index"]] <- vals
  }, error = function(e) NULL)
  seurat_obj
}

#' Post-read structural verification
#'
#' Read-back discipline for the reverse conversion: instead of trusting that
#' the load succeeded, assert what can be asserted. Failures raise
#' \code{scConvert_data_error}.
#' \itemize{
#'   \item shape/orientation: object dims must equal the file's (n_var, n_obs)
#'   \item cell identity/order: colnames must be exactly the obs index
#'   \item feature identity/order: rownames must be the var index modulo
#'     Seurat's documented underscore-to-dash replacement (and the
#'     make.unique that can follow it)
#'   \item uniqueness: no duplicate barcodes or genes may survive
#' }
#'
#' @param object Loaded Seurat object
#' @param cell.names Expected cell names (file obs index, post-dedup)
#' @param feature.names Expected feature names (file var index, post-dedup)
#' @param file Source path, for error messages
#' @return invisible(TRUE); raises scConvert_data_error otherwise
#' @keywords internal
#' @noRd
.h5ad_verify_read <- function(object, cell.names, feature.names, file) {
  n_obs <- length(cell.names)
  n_var <- length(feature.names)
  if (ncol(object) != n_obs || nrow(object) != n_var) {
    stop(.scconvert_data_error(sprintf(
      paste0("Post-read verification failed for '%s': loaded object is ",
             "%d features x %d cells but the file declares %d var x %d obs. ",
             "This indicates a transpose/orientation fault or a malformed file."),
      file, nrow(object), ncol(object), n_var, n_obs)))
  }
  if (!identical(colnames(object), cell.names)) {
    stop(.scconvert_data_error(sprintf(
      paste0("Post-read verification failed for '%s': cell names/order in the ",
             "loaded object do not match the file's obs index."),
      file)))
  }
  actual_features <- rownames(object)
  mangled <- gsub("_", "-", feature.names, fixed = TRUE)
  features_ok <- identical(actual_features, feature.names) ||
    identical(actual_features, mangled) ||
    identical(actual_features, make.unique(mangled))
  if (!features_ok) {
    stop(.scconvert_data_error(sprintf(
      paste0("Post-read verification failed for '%s': feature names/order in ",
             "the loaded object do not match the file's var index."),
      file)))
  }
  if (anyDuplicated(colnames(object)) > 0L) {
    stop(.scconvert_data_error(sprintf(
      "Post-read verification failed for '%s': duplicate cell barcodes in loaded object.",
      file)))
  }
  if (anyDuplicated(actual_features) > 0L) {
    stop(.scconvert_data_error(sprintf(
      "Post-read verification failed for '%s': duplicate feature names in loaded object.",
      file)))
  }
  invisible(TRUE)
}

#' Record what the reader actually did in misc$scConvert_read
#'
#' So callers (and any read-back verification downstream) can branch on facts
#' rather than trusting \code{readH5AD()}'s success return: which slot became
#' the counts layer, where X went, whether names had to be de-duplicated, and
#' which scConvert version performed the read.
#'
#' @keywords internal
#' @noRd
.h5ad_record_provenance <- function(object, file, counts_source, x_mapped_to,
                                    dedup_cells = FALSE, dedup_features = FALSE) {
  object@misc[["scConvert_read"]] <- list(
    package_version = as.character(utils::packageVersion("scConvert")),
    source_file = tryCatch(normalizePath(file), error = function(e) file),
    source_format = "h5ad",
    counts_source = counts_source,
    x_mapped_to = x_mapped_to,
    n_obs = ncol(object),
    n_var = nrow(object),
    duplicate_cell_names_renamed = isTRUE(dedup_cells),
    duplicate_feature_names_renamed = isTRUE(dedup_features),
    verified = TRUE
  )
  object
}

#' Stamp /uns/scConvert provenance into a freshly-written h5ad file
#'
#' Records the writing package version and -- by inspecting the file itself,
#' not the writer's intent -- where the raw counts landed
#' (\code{layers/counts}, \code{raw/X}, \code{X}, or \code{none}). Readers use
#' the stamp as the authoritative, version-aware branch for counts resolution,
#' eliminating the layout guessing that made older files ambiguous.
#'
#' Best-effort: returns invisible(FALSE) on any failure rather than failing
#' the write.
#'
#' @param target An open hdf5r H5File in write mode, or a file path
#' @return invisible(TRUE) on success
#' @keywords internal
#' @noRd
.h5ad_stamp_provenance <- function(target) {
  own_handle <- is.character(target)
  h5 <- if (own_handle) {
    tryCatch(hdf5r::H5File$new(target, mode = "r+"), error = function(e) NULL)
  } else {
    target
  }
  if (is.null(h5)) return(invisible(FALSE))

  ok <- tryCatch({
    counts_location <- if (.h5ad_slot_exists(h5, "layers/counts")) {
      "layers/counts"
    } else if (.h5ad_slot_exists(h5, "raw/X")) {
      "raw/X"
    } else if (isTRUE(h5$exists("X"))) {
      "X"
    } else {
      "none"
    }

    if (!h5$exists("uns")) {
      uns <- h5$create_group("uns")
      AddAnndataEncoding(uns, encoding_type = "dict", encoding_version = "0.1.0")
    }
    uns <- h5[["uns"]]
    if (uns$exists("scConvert")) uns$link_delete("scConvert")
    sg <- uns$create_group("scConvert")
    AddAnndataEncoding(sg, encoding_type = "dict", encoding_version = "0.1.0")
    # Scalar strings are written the same way WriteUnsItem writes them (no
    # element encoding attr): anndata's legacy fallback reads them cleanly,
    # which the existing python-validation CI already exercises.
    version_str <- as.character(utils::packageVersion("scConvert"))
    sg$create_dataset("version", robj = version_str,
                      dtype = CachedUtf8Type(), chunk_dims = 1L)
    sg$create_dataset("counts_location", robj = counts_location,
                      dtype = CachedUtf8Type(), chunk_dims = 1L)
    TRUE
  }, error = function(e) FALSE)

  if (own_handle) tryCatch(h5$close_all(), error = function(e) NULL)
  invisible(isTRUE(ok))
}
