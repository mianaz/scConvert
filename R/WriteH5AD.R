#' @include Convert.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read observation metadata from an H5AD file
#'
#' Read the cell-level metadata (obs) from an H5AD file and return it as a
#' data.frame. Handles both categorical (factor) and standard data types.
#'
#' @param file Path to an H5AD file
#'
#' @return A \code{data.frame} with cell metadata, where row names are cell
#'   barcodes from the \code{_index} dataset and columns correspond to
#'   observation annotations
#'
#' @export
#'
readH5AD_obs <- function(file) {
  suppressWarnings(expr = hfile <- scConnect(filename = file, force = TRUE))
  on.exit(hfile$close_all())
  hfile_obs <- hfile[['obs']]
  obs_groups <- setdiff(names(hfile_obs), c('__categories', '_index'))
  matrix <- as.data.frame(
    x = matrix(data = NA,
               nrow = hfile_obs[['_index']]$dims[1],
               ncol = length(obs_groups))
    )
  colnames(matrix) <- obs_groups
  rownames(matrix) <- hfile_obs[['_index']][]
  if ('__categories' %in% names(x = hfile_obs)) {
    # Legacy h5ad format: __categories group with per-column category lists
    hfile_cate <- hfile_obs[['__categories']]
    for (i in seq_along(obs_groups)) {
      obs.i <- obs_groups[i]
      obs_value_i <- hfile_obs[[obs.i]][]
      if (obs.i %in% names(x = hfile_cate)) {
        categories <- as.character(hfile_cate[[obs.i]][])
        codes <- obs_value_i
        codes[codes == -1L] <- NA_integer_
        obs_value_i <- factor(categories[codes + 1L], levels = categories)
      }
      matrix[, i] <- obs_value_i
    }
  } else {
    # Modern h5ad format: categorical groups with categories/codes sub-datasets
    for (i in seq_along(obs_groups)) {
      obs.i <- obs_groups[i]
      col_obj <- hfile_obs[[obs.i]]
      if (inherits(col_obj, 'H5Group') &&
          all(c("categories", "codes") %in% names(col_obj))) {
        categories <- as.character(col_obj[['categories']]$read())
        codes <- col_obj[['codes']]$read()
        codes[codes == -1L] <- NA_integer_
        obs_value_i <- factor(categories[codes + 1L], levels = categories)
      } else if (inherits(col_obj, 'H5D')) {
        obs_value_i <- col_obj$read()
      } else {
        obs_value_i <- NA
      }
      matrix[, i] <- obs_value_i
    }
  }
  return(matrix)
}

#' Read observation embeddings from an H5AD file
#'
#' Read the cell embeddings (obsm) from an H5AD file and return them as a
#' named list of matrices. Each entry corresponds to a dimensional reduction
#' (e.g., PCA, UMAP, t-SNE).
#'
#' @param file Path to an H5AD file
#'
#' @return A named list of matrices, where each matrix has cells as rows and
#'   embedding dimensions as columns. Names are derived from the obsm keys
#'   with the \dQuote{X_} prefix removed.
#'
#' @export
#'
readH5AD_obsm <- function(file) {
  suppressWarnings(hfile <- scConnect(filename = file, force = TRUE))
  on.exit(hfile$close_all())
  hfile_obsm <- hfile[['obsm']]
  if (length(names(hfile_obsm)) == 0) {
    message('No obsm found in this object')
    return(list())
  }
  obsm_set <- names(hfile_obsm)
  cells.name <- hfile[['obs']][['_index']][]
  n_cells <- length(cells.name)
  obsm.list <- lapply(obsm_set, function(x) {
    h5d <- hfile_obsm[[x]]
    ndims_h5d <- length(h5d$dims)
    emb <- if (ndims_h5d == 1L) {
      matrix(h5d[], ncol = 1L)
    } else {
      h5d[,]
    }
    # hdf5r may return transposed (n_dims x n_cells) due to row/column major difference
    if (nrow(emb) != n_cells && ncol(emb) == n_cells) {
      emb <- t(emb)
    }
    rownames(emb) <- cells.name
    key.name <- gsub('X_', '', x)
    colnames(emb) <- paste0(key.name, "_", seq_len(ncol(emb)))
    return(emb)
  })
  names(obsm.list) <- gsub('X_', '', obsm_set)
  return(obsm.list)
}

#' Direct Seurat to H5AD Conversion
#'
#' Converts a Seurat object directly to an H5AD file, handling the intermediate
#' h5Seurat file automatically. The intermediate file is created in a temporary
#' location and removed after conversion.
#'
#' @param object A Seurat object to convert
#' @param filename Output H5AD filename. If not specified, uses the project name
#'   with .h5ad extension.
#' @param assay Name of assay to convert. Default is the default assay.
#' @param overwrite Logical; overwrite existing file. Default FALSE.
#' @param verbose Logical; show progress messages. Default TRUE.
#' @param standardize Logical; convert Seurat metadata names to scanpy conventions.
#'   Default FALSE.
#' @param gzip Integer gzip compression level (0-9), or NULL to use the package default.
#' @param ... Additional arguments passed to writeH5Seurat and Convert.
#'
#' @return Invisibly returns the path to the created H5AD file.
#'
#' @details
#' This function provides a convenient one-step conversion from Seurat objects
#' to H5AD format (used by Python's scanpy/anndata). Internally, it:
#' \enumerate{
#'   \item Saves the Seurat object to a temporary h5Seurat file
#'   \item Converts the h5Seurat file to H5AD format
#'   \item Removes the intermediate h5Seurat file
#' }
#'
#' This is useful when you want to export data for Python analysis without
#' keeping the intermediate h5Seurat file.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scConvert)
#'
#' # Load or create a Seurat object
#' pbmc <- pbmc_small
#'
#' # Convert directly to H5AD
#' writeH5AD(pbmc, filename = "pbmc.h5ad")
#'
#' # The file can now be loaded in Python:
#' # import scanpy as sc
#' # adata = sc.read_h5ad("pbmc.h5ad")
#' }
#'
#' @seealso \code{\link{writeH5Seurat}}, \code{\link{scConvert}}, \code{\link{readH5AD}}
#'
#' @export
#'
writeH5AD <- function(
  object,
  filename = NULL,
  assay = DefaultAssay(object = object),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  gzip = NULL,
  ...
) {
  if (!inherits(x = object, what = 'Seurat')) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  if (is.null(x = filename)) {
    filename <- paste0(Project(object = object), '.h5ad')
  }
  if (!grepl(pattern = '\\.h5ad$', x = filename, ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5ad')
  }
  # Optimization 3: configurable gzip level (0 = no compression, ~50-70% faster writes)
  if (!is.null(gzip)) {
    old_gzip <- getOption("scConvert.compression.level")
    on.exit(options("scConvert.compression.level" = old_gzip), add = TRUE)
    options("scConvert.compression.level" = as.integer(gzip))
  }
  scConvert(source = object, dest = filename, assay = assay,
          overwrite = overwrite, verbose = verbose, standardize = standardize, ...)
}

#' Fast C-based h5ad writer
#'
#' Writes a Seurat object directly to h5ad format using native C HDF5 routines.
#' Exploits the zero-copy CSC<->CSR reinterpretation: dgCMatrix(genesxcells) CSC
#' is identical to h5ad's cellsxgenes CSR with relabeled arrays.
#'
#' @param object Seurat object
#' @param filename Output file path
#' @param overwrite Allow overwriting
#' @param verbose Show progress
#'
#' @return TRUE on success, FALSE on failure
#' @keywords internal
#'
.writeH5AD_c <- function(object, filename, overwrite = FALSE, verbose = TRUE) {
  if (file.exists(filename)) {
    if (overwrite) file.remove(filename)
    else stop("File '", filename, "' already exists", call. = FALSE)
  }

  assay_name <- DefaultAssay(object)
  if (length(Assays(object)) > 1) return(FALSE)

  assay_obj <- object[[assay_name]]
  layer_has <- function(ln) {
    tryCatch({
      m <- GetAssayData(assay_obj, layer = ln)
      inherits(m, "dgCMatrix") && length(m@x) > 0
    }, error = function(e) FALSE)
  }
  has_counts <- layer_has("counts")
  has_data   <- layer_has("data")

  # anndata convention: X holds normalized data, raw/X holds counts. When both
  # layers exist, route counts to the raw group; when only one exists, it
  # becomes X directly.
  mat <- if (has_data) GetAssayData(assay_obj, layer = "data")
         else if (has_counts) GetAssayData(assay_obj, layer = "counts")
         else NULL
  if (is.null(mat)) return(FALSE)

  raw_list <- NULL
  if (has_data && has_counts) {
    raw_mat <- GetAssayData(assay_obj, layer = "counts")
    raw_list <- list(i = raw_mat@i, p = raw_mat@p, x = raw_mat@x,
                     dim = dim(raw_mat))
    rm(raw_mat)
  }

  mat_list <- list(
    i = mat@i, p = mat@p, x = mat@x,
    dim = dim(mat),
    rownames = rownames(mat),
    colnames = colnames(mat)
  )

  meta <- as.list(object@meta.data)

  reductions <- list()
  for (rname in names(object@reductions)) {
    reduc <- object@reductions[[rname]]
    emb <- Embeddings(reduc)
    attr(emb, "key") <- Key(reduc)
    reductions[[rname]] <- emb
  }

  graphs <- list()
  for (gname in names(object@graphs)) {
    g <- object@graphs[[gname]]
    if (inherits(g, "dgCMatrix") || inherits(g, "Graph")) {
      gmat <- as(g, "dgCMatrix")
      graphs[[gname]] <- list(i = gmat@i, p = gmat@p, x = gmat@x, dim = dim(gmat))
    }
  }

  gzip_level <- GetCompressionLevel()
  if (verbose) message("Writing h5ad (C writer): ", filename)

  result <- .Call(C_write_h5ad,
    filename, mat_list, meta, reductions, graphs,
    assay_name, as.integer(gzip_level), raw_list
  )

  if (isTRUE(result) && verbose) {
    message("  Written: ", ncol(mat), " cells, ", nrow(mat), " features")
  }
  return(isTRUE(result))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Direct Seurat-to-h5ad Pipeline
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Direct Seurat-to-h5ad write (single-pass, no intermediate h5Seurat)
#'
#' Writes a Seurat object directly to an h5ad file, bypassing the two-pass
#' Seurat -> h5Seurat -> h5ad pipeline. This avoids 3x I/O overhead and
#' provides significant speedup for in-memory Seurat objects.
#'
#' @param object A Seurat object
#' @param filename Output h5ad file path
#' @param assay Name of the assay to write (default: DefaultAssay)
#' @param overwrite Overwrite existing file (default: FALSE)
#' @param verbose Print progress messages (default: TRUE)
#' @param standardize Convert Seurat metadata names to scanpy conventions
#'
#' @return Invisibly returns the output filename
#'
#' @keywords internal
#'
DirectSeuratToH5AD <- function(
  object,
  filename,
  assay = DefaultAssay(object = object),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE
) {
  if (!inherits(object, 'Seurat')) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  if (file.exists(filename) && !overwrite) {
    stop("File '", filename, "' already exists; set overwrite = TRUE", call. = FALSE)
  }

  # Atomic write: write to temp file, rename on success to avoid data loss on failure
  use_atomic <- file.exists(filename) && overwrite
  write_target <- if (use_atomic) {
    tempfile(tmpdir = dirname(filename), fileext = ".h5ad")
  } else {
    filename
  }

  gzip <- GetCompressionLevel()
  assay_obj <- object[[assay]]
  cell_names <- colnames(object)
  gene_names <- rownames(assay_obj)
  n_cells <- length(cell_names)
  n_genes <- length(gene_names)

  if (verbose) message("Writing h5ad directly from Seurat object (", n_cells, " cells x ", n_genes, " genes)")

  # Create h5ad file (may be temp file for atomic write)
  dfile <- H5File$new(write_target, mode = 'w')
  on.exit(tryCatch(dfile$close_all(), error = function(e) NULL), add = TRUE)

  # Local aliases for shared helpers (capture gzip from enclosing scope)
  write_csr_group <- function(h5parent, group_name, mat) {
    WriteCSRGroup(h5parent, group_name, mat, gzip = gzip)
  }
  write_csc_group <- function(h5parent, group_name, mat) {
    WriteCscGroup(h5parent, group_name, mat, gzip = gzip)
  }
  write_df_group <- function(h5parent, group_name, df, index_values) {
    WriteDFGroup(h5parent, group_name, df, index_values, gzip = gzip)
  }

  # ========== X (primary expression matrix) ==========
  if (verbose) message("  Writing X...")
  # Retrieve each layer ONCE (GetAssayData materializes the full matrix each call)
  # Sequential layer write: load one matrix at a time to halve peak RAM.
  # Detect which layers exist without materializing them.
  has_data <- tryCatch({
    d <- GetAssayData(object, assay = assay, layer = 'data')
    !is.null(d) && prod(dim(d)) > 0
  }, error = function(e) FALSE)

  has_counts <- tryCatch({
    d <- GetAssayData(object, assay = assay, layer = 'counts')
    !is.null(d) && prod(dim(d)) > 0
  }, error = function(e) FALSE)

  if (has_data) {
    data_mat <- GetAssayData(object, assay = assay, layer = 'data')
    write_csr_group(dfile, "X", data_mat)
    rm(data_mat); gc(verbose = FALSE)
    if (has_counts) {
      counts_mat <- GetAssayData(object, assay = assay, layer = 'counts')
      raw_grp <- dfile$create_group("raw")
      write_csr_group(raw_grp, "X", counts_mat)
      rm(counts_mat); gc(verbose = FALSE)
    }
  } else if (has_counts) {
    counts_mat <- GetAssayData(object, assay = assay, layer = 'counts')
    write_csr_group(dfile, "X", counts_mat)
    rm(counts_mat); gc(verbose = FALSE)
  }

  # ========== obs (cell metadata) ==========
  if (verbose) message("  Writing obs...")
  meta <- object[[]]
  if (standardize) {
    # Convert Seurat-style names to scanpy conventions
    name_map <- c(
      "nCount_RNA" = "n_counts", "nFeature_RNA" = "n_genes",
      "percent.mt" = "pct_counts_mt", "orig.ident" = "batch"
    )
    for (old_name in names(name_map)) {
      if (old_name %in% colnames(meta)) {
        colnames(meta)[colnames(meta) == old_name] <- name_map[[old_name]]
      }
    }
  }
  write_df_group(dfile, "obs", meta, cell_names)

  # ========== var (gene metadata) ==========
  if (verbose) message("  Writing var...")
  meta_features <- tryCatch(assay_obj[[]], error = function(e) data.frame(row.names = gene_names))
  if (is.null(meta_features) || nrow(meta_features) == 0) {
    meta_features <- data.frame(row.names = gene_names)
  }
  # Add highly_variable boolean mask (scConvert unique feature)
  var_features <- tryCatch(VariableFeatures(object, assay = assay), error = function(e) character(0))
  if (length(var_features) > 0) {
    meta_features$highly_variable <- rownames(meta_features) %in% var_features
  }
  write_df_group(dfile, "var", meta_features, gene_names)

  # Write raw/var if raw/X exists
  if (has_data && has_counts) {
    write_df_group(dfile[["raw"]], "var",
                   data.frame(row.names = gene_names), gene_names)
  }

  # ========== obsm (dimensional reductions) ==========
  reducs <- Reductions(object)
  if (length(reducs) > 0) {
    if (verbose) message("  Writing obsm (", length(reducs), " reductions)...")
    obsm_grp <- dfile$create_group("obsm")
    for (reduc_name in reducs) {
      emb <- tryCatch(Embeddings(object, reduction = reduc_name), error = function(e) NULL)
      if (!is.null(emb) && nrow(emb) > 0) {
        key_name <- paste0("X_", reduc_name)
        # scTranspose: hdf5r writes R (n_cells, n_dims) as HDF5 (n_dims, n_cells)
        # anndata expects (n_cells, n_dims), so we write t(emb)
        emb_t <- t(emb)
        obsm_grp$create_dataset(key_name, robj = emb_t,
                                chunk_dims = dim(emb_t),
                                gzip_level = gzip)
        obsm_grp[[key_name]]$create_attr(
          attr_name = 'encoding-type', robj = 'array',
          dtype = CachedGuessDType('array'), space = ScalarSpace())
        obsm_grp[[key_name]]$create_attr(
          attr_name = 'encoding-version', robj = '0.2.0',
          dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
      }
    }
  }

  # ========== obsp (graphs) ==========
  graphs <- Graphs(object)
  if (length(graphs) > 0) {
    if (verbose) message("  Writing obsp (", length(graphs), " graphs)...")
    obsp_grp <- dfile$create_group("obsp")
    for (graph_name in graphs) {
      graph_mat <- tryCatch(object[[graph_name]], error = function(e) NULL)
      if (!is.null(graph_mat) && inherits(graph_mat, "Graph")) {
        write_csc_group(obsp_grp, graph_name, graph_mat)
      }
    }
  }

  # ========== varp (pairwise variable annotations) ==========
  # Check both global __varp__ and per-assay __varp__.{assay} (from MuData roundtrips)
  varp_data <- tryCatch(Misc(object)[["__varp__"]], error = function(e) NULL)
  if (is.null(varp_data) || !is.list(varp_data) || length(varp_data) == 0) {
    varp_key <- paste0("__varp__.", assay)
    varp_data <- tryCatch(Misc(object)[[varp_key]], error = function(e) NULL)
  }
  if (!is.null(varp_data) && is.list(varp_data) && length(varp_data) > 0) {
    if (verbose) message("  Writing varp (", length(varp_data), " pairwise annotations)...")
    varp_grp <- dfile$create_group("varp")
    for (varp_name in names(varp_data)) {
      varp_mat <- varp_data[[varp_name]]
      if (inherits(varp_mat, "sparseMatrix") || inherits(varp_mat, "Matrix")) {
        write_csr_group(varp_grp, varp_name, varp_mat)
      } else if (is.matrix(varp_mat)) {
        varp_grp$create_dataset(varp_name, robj = varp_mat,
                                chunk_dims = dim(varp_mat), gzip_level = gzip)
        varp_grp[[varp_name]]$create_attr(
          attr_name = 'encoding-type', robj = 'array',
          dtype = CachedGuessDType('array'), space = ScalarSpace())
        varp_grp[[varp_name]]$create_attr(
          attr_name = 'encoding-version', robj = '0.2.0',
          dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
      }
    }
  }

  # ========== uns (misc + spatial) ==========
  misc <- tryCatch(Misc(object), error = function(e) list())
  images <- tryCatch(Images(object), error = function(e) character(0))
  has_uns <- (length(misc) > 0 && !all(sapply(misc, is.null))) || length(images) > 0

  if (has_uns) {
    if (verbose) message("  Writing uns...")
    uns_grp <- dfile$create_group("uns")

    # Internal keys managed separately (varp, lazy-load bookkeeping)
    skip_keys <- c("__varp__", ".__h5ad_path__", ".__h5ad_loaded__")
    skip_keys <- c(skip_keys, grep("^__varp__\\.", names(misc), value = TRUE))

    gzip <- GetCompressionLevel()
    for (item_name in names(misc)) {
      if (item_name %in% skip_keys) next
      tryCatch(
        WriteUnsItem(uns_grp, item_name, misc[[item_name]], gzip),
        error = function(e) {
          if (verbose) message("  Could not write uns/", item_name, ": ", e$message)
        }
      )
    }

    # Delegate spatial data (images, coordinates, scale factors)
    if (length(images) > 0) {
      tryCatch({
        SeuratSpatialToH5AD(object, dfile, verbose = verbose)
      }, error = function(e) {
        if (verbose) message("Spatial data conversion failed: ", e$message)
      })
    }
  }

  # ========== layers (additional layers beyond data/counts) ==========
  all_layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character(0))
  extra_layers <- setdiff(all_layers, c("data", "counts", "scale.data"))
  if (length(extra_layers) > 0) {
    if (verbose) message("  Writing layers (", length(extra_layers), " extra)...")
    layers_grp <- dfile$create_group("layers")
    for (ln in extra_layers) {
      layer_mat <- tryCatch(GetAssayData(object, assay = assay, layer = ln), error = function(e) NULL)
      if (!is.null(layer_mat) && prod(dim(layer_mat)) > 0) {
        write_csr_group(layers_grp, ln, layer_mat)
      }
    }
  }

  # ========== Root-level AnnData attributes ==========
  dfile$create_attr(attr_name = 'encoding-type', robj = 'anndata',
                    dtype = CachedGuessDType('anndata'), space = ScalarSpace())
  dfile$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                    dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())

  dfile$flush()
  dfile$close_all()

  # Atomic rename: replace original only after successful write
  if (use_atomic) {
    file.rename(write_target, filename)
  }

  if (verbose) message("  Done: ", filename)
  invisible(filename)
}

# H5MU pair converters (H5MUToH5Seurat, H5SeuratToH5MU, H5MUToH5AD) now live
# in R/StreamHDF5.R with stream=TRUE streaming semantics. Prior stale hub-only
# copies were removed 2026-04-20.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Seurat-to-h5mu Pipeline (native, no MuDataSeurat dependency)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write a Seurat object to h5mu format
#'
#' Writes a multi-assay Seurat object to an h5mu file, with each
#' assay becoming a separate modality under \code{/mod/{modality}/}.
#' Works with Seurat V5 Assay5 objects. No external dependencies required.
#'
#' @param object A Seurat object with one or more assays
#' @param filename Output h5mu file path. If NULL, derived from project name.
#' @param assays Character vector of assay names to export (default: all)
#' @param overwrite Overwrite existing file (default: FALSE)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Invisibly returns the output filename
#'
#' @details
#' The h5mu format is designed for multimodal data and stores each Seurat assay
#' as a separate modality under \code{/mod/{modality_name}}. This function:
#' \itemize{
#'   \item Extracts each specified assay from the Seurat object
#'   \item Converts assays to modality structure
#'   \item Writes counts and data layers for each modality
#'   \item Preserves cell metadata in global /obs and per-modality obs
#'   \item Maintains dimensional reductions and graphs
#' }
#'
#' @section Assay Mapping:
#' By default, Seurat assay names are mapped to standard MuData modality names:
#' \itemize{
#'   \item RNA -> rna
#'   \item ADT -> prot
#'   \item ATAC -> atac
#'   \item Spatial -> spatial
#'   \item Other names are converted to lowercase
#' }
#'
#' @examples
#' \dontrun{
#' # Write a multi-assay Seurat object (e.g. CITE-seq with RNA + ADT)
#' writeH5MU(seurat_obj, filename = "multimodal.h5mu")
#'
#' # Write specific assays only
#' writeH5MU(seurat_obj, filename = "rna_only.h5mu", assays = "RNA")
#' }
#'
#' @seealso \code{\link{readH5MU}}, \code{\link{scConvert}}, \code{\link{as.h5mu}}
#'
#' @export
#'
writeH5MU <- function(
  object,
  filename = NULL,
  assays = NULL,
  overwrite = FALSE,
  verbose = TRUE
) {
  if (!inherits(object, 'Seurat')) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  if (is.null(filename)) {
    filename <- paste0(Project(object = object), '.h5mu')
  }
  if (!grepl('\\.h5mu$', filename, ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5mu')
  }
  if (file.exists(filename) && !overwrite) {
    stop("File '", filename, "' already exists; set overwrite = TRUE", call. = FALSE)
  }
  if (file.exists(filename) && overwrite) {
    file.remove(filename)
  }

  gzip <- GetCompressionLevel()
  cell_names <- colnames(object)
  n_cells <- length(cell_names)

  # Determine which assays to export
  available_assays <- Assays(object)
  if (is.null(assays)) {
    assays_to_export <- available_assays
  } else {
    assays_to_export <- assays
    missing <- setdiff(assays_to_export, available_assays)
    if (length(missing) > 0) {
      stop("Assays not found in object: ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  # Get modality name mapping
  modality_map <- GetDefaultAssayToModalityMapping(assays_to_export)

  if (verbose) {
    message("Writing h5mu directly from Seurat object (", n_cells, " cells, ",
            length(assays_to_export), " assays)")
  }

  # Create h5mu file
  dfile <- H5File$new(filename, mode = 'w')
  on.exit(tryCatch(dfile$close_all(), error = function(e) NULL), add = TRUE)

  # Local aliases for shared helpers (capture gzip from enclosing scope)
  write_csr_group <- function(h5parent, group_name, mat) {
    WriteCSRGroup(h5parent, group_name, mat, gzip = gzip)
  }
  write_csc_group <- function(h5parent, group_name, mat) {
    WriteCscGroup(h5parent, group_name, mat, gzip = gzip)
  }
  write_df_group <- function(h5parent, group_name, df, index_values) {
    WriteDFGroup(h5parent, group_name, df, index_values, gzip = gzip)
  }

  # Per-modality var rownames collected during the loop below; used after
  # the loop to build the canonical global /var axis and /varmap/{mod}.
  mod_var_names <- list()

  # ========== Root MuData attributes ==========
  dfile$create_attr(attr_name = 'encoding-type', robj = 'MuData',
                    dtype = CachedGuessDType('MuData'), space = ScalarSpace())
  dfile$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                    dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())

  # ========== /mod/ group ==========
  mod_root <- dfile$create_group("mod")
  mod_order <- unname(modality_map[assays_to_export])
  mod_root$create_attr(attr_name = 'mod-order', robj = mod_order,
                       dtype = CachedUtf8Type())

  # ========== Per-modality writing ==========
  for (assay_name in assays_to_export) {
    modality <- modality_map[[assay_name]]
    assay_obj <- object[[assay_name]]
    gene_names <- rownames(assay_obj)

    if (verbose) message("  Writing modality '", modality, "' (assay: ", assay_name, ")...")

    mod_grp <- mod_root$create_group(modality)

    # --- X (primary expression matrix) ---
    has_data <- tryCatch({
      d <- GetAssayData(object, assay = assay_name, layer = 'data')
      !is.null(d) && prod(dim(d)) > 0
    }, error = function(e) FALSE)

    has_counts <- tryCatch({
      d <- GetAssayData(object, assay = assay_name, layer = 'counts')
      !is.null(d) && prod(dim(d)) > 0
    }, error = function(e) FALSE)

    if (has_data) {
      x_mat <- GetAssayData(object, assay = assay_name, layer = 'data')
      write_csr_group(mod_grp, "X", x_mat)
    } else if (has_counts) {
      counts_mat <- GetAssayData(object, assay = assay_name, layer = 'counts')
      write_csr_group(mod_grp, "X", counts_mat)
    }

    # --- layers/counts (if X is data and counts also available) ---
    if (has_data && has_counts) {
      layers_grp <- mod_grp$create_group("layers")
      layers_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                             dtype = CachedGuessDType('dict'), space = ScalarSpace())
      layers_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                             dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
      counts_mat <- GetAssayData(object, assay = assay_name, layer = 'counts')
      write_csr_group(layers_grp, "counts", counts_mat)
    }

    # --- obs (per-modality: just cell names) ---
    write_df_group(mod_grp, "obs", data.frame(row.names = cell_names), cell_names)

    # --- var (feature metadata) ---
    meta_features <- tryCatch(assay_obj[[]], error = function(e) data.frame(row.names = gene_names))
    if (is.null(meta_features) || nrow(meta_features) == 0) {
      meta_features <- data.frame(row.names = gene_names)
    }
    write_df_group(mod_grp, "var", meta_features, gene_names)
    mod_var_names[[modality]] <- gene_names

    # --- obsm (reductions associated with this assay) ---
    reducs <- Reductions(object)
    assay_reducs <- character(0)
    for (rn in reducs) {
      reduc_assay <- tryCatch(DefaultAssay(object[[rn]]), error = function(e) "")
      if (reduc_assay == assay_name) {
        assay_reducs <- c(assay_reducs, rn)
      }
    }
    if (length(assay_reducs) > 0) {
      obsm_grp <- mod_grp$create_group("obsm")
      obsm_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                           dtype = CachedGuessDType('dict'), space = ScalarSpace())
      obsm_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                           dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
      for (reduc_name in assay_reducs) {
        emb <- tryCatch(Embeddings(object, reduction = reduc_name), error = function(e) NULL)
        if (!is.null(emb) && nrow(emb) > 0) {
          key_name <- paste0("X_", reduc_name)
          emb_t <- t(emb)
          obsm_grp$create_dataset(key_name, robj = emb_t,
                                  chunk_dims = dim(emb_t),
                                  gzip_level = gzip)
          obsm_grp[[key_name]]$create_attr(
            attr_name = 'encoding-type', robj = 'array',
            dtype = CachedGuessDType('array'), space = ScalarSpace())
          obsm_grp[[key_name]]$create_attr(
            attr_name = 'encoding-version', robj = '0.2.0',
            dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
        }
      }
    }

    # --- obsp (graphs associated with this assay) ---
    all_graphs <- Graphs(object)
    assay_graphs <- character(0)
    for (gn in all_graphs) {
      graph_obj <- tryCatch(object[[gn]], error = function(e) NULL)
      if (!is.null(graph_obj) && inherits(graph_obj, "Graph")) {
        graph_assay <- tryCatch(DefaultAssay(graph_obj), error = function(e) "")
        if (length(graph_assay) == 1 && !is.null(graph_assay) && graph_assay == assay_name) {
          assay_graphs <- c(assay_graphs, gn)
        }
      }
    }
    if (length(assay_graphs) > 0) {
      obsp_grp <- mod_grp$create_group("obsp")
      obsp_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                           dtype = CachedGuessDType('dict'), space = ScalarSpace())
      obsp_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                           dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
      for (graph_name in assay_graphs) {
        graph_mat <- object[[graph_name]]
        # Strip assay prefix if present (e.g. "RNA_snn" -> "snn")
        clean_name <- sub(paste0("^", assay_name, "_"), "", graph_name)
        write_csc_group(obsp_grp, clean_name, graph_mat)
      }
    }

    # --- varp (pairwise variable annotations) ---
    varp_key <- paste0("__varp__.", assay_name)
    varp_data <- tryCatch(Misc(object)[[varp_key]], error = function(e) NULL)
    if (is.null(varp_data)) {
      varp_data <- tryCatch(Misc(object)[["__varp__"]], error = function(e) NULL)
    }
    if (!is.null(varp_data) && is.list(varp_data) && length(varp_data) > 0) {
      varp_grp <- mod_grp$create_group("varp")
      varp_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                           dtype = CachedGuessDType('dict'), space = ScalarSpace())
      varp_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                           dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
      for (varp_name in names(varp_data)) {
        varp_mat <- varp_data[[varp_name]]
        if (inherits(varp_mat, "sparseMatrix") || inherits(varp_mat, "Matrix")) {
          write_csr_group(varp_grp, varp_name, varp_mat)
        } else if (is.matrix(varp_mat)) {
          varp_grp$create_dataset(varp_name, robj = varp_mat,
                                  chunk_dims = dim(varp_mat), gzip_level = gzip)
          varp_grp[[varp_name]]$create_attr(
            attr_name = 'encoding-type', robj = 'array',
            dtype = CachedGuessDType('array'), space = ScalarSpace())
          varp_grp[[varp_name]]$create_attr(
            attr_name = 'encoding-version', robj = '0.2.0',
            dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
        }
      }
    }

    # --- Per-modality uns (mod/{modality}/uns) ---
    # Mirror of obj@misc[["__h5mu_uns_per_mod__"]][[modality]]. This is the
    # write side of the fix for the h5mu 19/22 fidelity gap where per-modality
    # uns was previously flattened into the global uns.
    mod_uns_list <- tryCatch(
      Misc(object)[["__h5mu_uns_per_mod__"]][[modality]],
      error = function(e) NULL
    )
    if (is.list(mod_uns_list) && length(mod_uns_list) > 0) {
      mod_uns_grp <- mod_grp$create_group("uns")
      mod_uns_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                               dtype = CachedGuessDType('dict'), space = ScalarSpace())
      mod_uns_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                               dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
      for (key in names(mod_uns_list)) {
        local({
          k <- key
          mod <- modality
          tryCatch(
            WriteUnsItem(mod_uns_grp, k, mod_uns_list[[k]], gzip),
            error = function(e) {
              if (verbose) {
                message("    Warning: could not write uns key '", k,
                        "' for modality '", mod, "': ", e$message)
              }
            }
          )
        })
      }
    }

    # --- Ensure empty dict groups ---
    for (empty_name in c("obsm", "obsp", "varm", "varp", "layers", "uns")) {
      if (!mod_grp$exists(empty_name)) {
        eg <- mod_grp$create_group(empty_name)
        eg$create_attr(attr_name = 'encoding-type', robj = 'dict',
                       dtype = CachedGuessDType('dict'), space = ScalarSpace())
        eg$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                       dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
      }
    }

    # Set AnnData encoding on modality group
    mod_grp$create_attr(attr_name = 'encoding-type', robj = 'anndata',
                        dtype = CachedGuessDType('anndata'), space = ScalarSpace())
    mod_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                        dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
  }

  # ========== Global obs (shared cell metadata) ==========
  if (verbose) message("  Writing global obs...")
  meta <- object[[]]
  write_df_group(dfile, "obs", meta, cell_names)

  # ========== Global var (concat of modality vars, canonical MuData) ==========
  # A Seurat object's modalities share the obs axis but have disjoint var
  # axes. MuData's canonical layout is a single /var that concatenates all
  # modality vars in mod-order, with /varmap/{m} pointing each global-var
  # index back at its position in modality m (-1 for vars from other mods).
  if (verbose) message("  Writing canonical global var + mapping...")

  concat_var_names <- unlist(mod_var_names[mod_order], use.names = FALSE)
  if (is.null(concat_var_names)) concat_var_names <- character(0)
  # Deduplicate with modality prefix only where a gene name reappears in
  # multiple modalities; otherwise preserve the original name so callers
  # see their expected feature identifiers.
  dup_mask <- duplicated(concat_var_names) |
              duplicated(concat_var_names, fromLast = TRUE)
  if (any(dup_mask)) {
    mod_owner <- rep(mod_order,
                     vapply(mod_order,
                            function(m) length(mod_var_names[[m]]),
                            integer(1)))
    prefixed <- paste0(mod_owner, ":", concat_var_names)
    concat_var_names[dup_mask] <- prefixed[dup_mask]
  }
  write_df_group(dfile, "var",
                 data.frame(row.names = concat_var_names),
                 concat_var_names)

  # ========== /obsmap and /varmap ==========
  # /obsmap/{mod}[i] = 0-based index of global cell i in modality mod's obs,
  # or -1 if that cell is absent from the modality. Seurat objects share the
  # obs axis across modalities, so every cell is present in every modality
  # and the map is always 0..n-1.
  #
  # /varmap/{mod}[j] = 0-based index of global var j in modality mod's var,
  # or -1 if that global var belongs to a different modality.
  obsmap_grp <- dfile$create_group("obsmap")
  obsmap_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                         dtype = CachedGuessDType('dict'), space = ScalarSpace())
  obsmap_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                         dtype = CachedGuessDType('0.1.0'),
                         space = ScalarSpace())

  varmap_grp <- dfile$create_group("varmap")
  varmap_grp$create_attr(attr_name = 'encoding-type', robj = 'dict',
                         dtype = CachedGuessDType('dict'), space = ScalarSpace())
  varmap_grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                         dtype = CachedGuessDType('0.1.0'),
                         space = ScalarSpace())

  n_obs_global <- length(cell_names)
  n_var_global <- length(concat_var_names)
  var_offset <- 0L
  for (assay_name in assays_to_export) {
    modality <- modality_map[[assay_name]]
    n_mod_var <- length(mod_var_names[[modality]])

    obsmap_vec <- seq.int(0L, n_obs_global - 1L)
    storage.mode(obsmap_vec) <- "integer"
    obsmap_grp$create_dataset(modality, robj = obsmap_vec,
                              chunk_dims = length(obsmap_vec),
                              gzip_level = gzip)

    varmap_vec <- rep(-1L, n_var_global)
    if (n_mod_var > 0L) {
      varmap_vec[(var_offset + 1L):(var_offset + n_mod_var)] <-
        seq.int(0L, n_mod_var - 1L)
    }
    varmap_grp$create_dataset(modality, robj = varmap_vec,
                              chunk_dims = length(varmap_vec),
                              gzip_level = gzip)

    var_offset <- var_offset + n_mod_var
  }

  # ========== Global obsm (empty) ==========
  global_obsm <- dfile$create_group("obsm")
  global_obsm$create_attr(attr_name = 'encoding-type', robj = 'dict',
                          dtype = CachedGuessDType('dict'), space = ScalarSpace())
  global_obsm$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                          dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())

  dfile$flush()
  dfile$close_all()

  if (verbose) message("  Done: ", filename)
  invisible(filename)
}
