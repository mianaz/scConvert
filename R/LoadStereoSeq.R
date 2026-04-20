#' @include zzz.R
#' @importFrom hdf5r H5File h5attr
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stereo-seq GEF reader (pure R, no reticulate)
#
# Implements native reading of BGI's GEF/cellbin.gef HDF5 schema documented at
# https://www.stomics.tech/service/saw_8_1/docs/gao-ji-she-zhi/expression-matrix-format.html
#
# Square-bin .gef:
#   /geneExp/bin{N}/gene        compound {geneID(S64), geneName(S64),
#                                         offset(u32), count(u32)}
#   /geneExp/bin{N}/expression  compound {x(i32), y(i32), count(uN)}
# Cell-bin .cellbin.gef:
#   /cellBin/cell    compound {id, x, y, offset, geneCount, expCount, ...}
#   /cellBin/gene    compound {geneID, geneName, offset, cellcount, expCount, ...}
#   /cellBin/cellExp compound {geneID, count}
#
# TODO(stereoseq, 2026-04-20): the reader is unit-tested against synthetic
# .gef fixtures but has not yet been validated end-to-end against a real
# MOSTA section because the CNGB STOmicsDB portal gates .gef downloads and
# the Zenodo MOSTA mirror (record 10698963) hosts truncated h5ad files.
# Remaining work before claiming Stereo-seq Task-5 fidelity:
#   1. Obtain a conforming public .gef exemplar (CNGB account, BGI rep, or
#      collaborator handoff) and stage it under data/raw/stereoseq/.
#   2. Extend benchmark/run_real_modalities.R with a Stereo-seq block that
#      ingests via LoadStereoSeqGef(), round-trips through writeH5AD /
#      readH5AD, and writes a Task-5 row.
#   3. Restore the stereoseq_mosta manifest entry in
#      data/dataset_manifest.yaml (see git log commit 8926b72 for the
#      previous entry, which noted the CNGB access path and stereopy
#      fallback).
#   4. Add a short paragraph to scConvert.tex §Spatial and §Limitations
#      reintroducing Stereo-seq as a validated vendor-native path.
# Coordinates / downstream spatial plotting against squidpy and stereopy
# AnnData outputs have not been cross-checked yet — add a fidelity check
# for spatial_x / spatial_y consistency in helpers/fidelity_schema.R when
# this work is resumed.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Convert an hdf5r compound-dataset string field to a character vector.
# hdf5r may return fixed-length strings as raw bytes (one row per cell) rather
# than R characters; this helper normalises either representation.
.read_compound_string <- function(field) {
  if (is.character(field)) return(field)
  if (is.factor(field))    return(as.character(field))
  if (is.list(field)) {
    return(vapply(field, function(el) {
      if (is.raw(el)) {
        keep <- el[el != as.raw(0)]
        if (length(keep) == 0) "" else rawToChar(keep)
      } else {
        as.character(el)
      }
    }, character(1)))
  }
  if (is.raw(field)) {
    keep <- field[field != as.raw(0)]
    if (length(keep) == 0) "" else rawToChar(keep)
  } else {
    as.character(field)
  }
}

# Abort early if the gene/expression tables are inconsistent. The GEF spec
# states sum(gene$count) == nrow(expression); a mismatch means the file is
# corrupt or truncated.
.validate_gef_consistency <- function(genes, expr, where = "geneExp") {
  expected <- sum(as.numeric(genes$count))
  observed <- nrow(expr)
  if (!isTRUE(expected == observed)) {
    stop(sprintf(
      "Corrupt Stereo-seq GEF: %s gene$count sums to %.0f but expression has %.0f rows",
      where, expected, observed), call. = FALSE)
  }
  invisible(TRUE)
}

# Read the root `bin_type` attribute to decide between square-bin and cell-bin.
# Missing attribute defaults to square-bin (older GEFs).
.detect_gef_bin_type <- function(h5) {
  if (!("bin_type" %in% hdf5r::list.attributes(h5))) return("bin")
  val <- tryCatch(h5attr(h5, "bin_type"), error = function(e) "bin")
  if (is.raw(val))       val <- rawToChar(val[val != as.raw(0)])
  if (is.list(val))      val <- val[[1]]
  if (is.null(val) || !nzchar(val)) return("bin")
  # BGI writes either "bin"/"b'bin'" for square and "CellBin"/"cellbin" for cell.
  if (grepl("cellbin", tolower(val), fixed = TRUE)) "cellbin" else "bin"
}

# Pack (x, y) integer coordinate pairs into a single 64-bit key without R's
# slow paste() path. Avoids the 10M+ spot RAM blowup described in the plan.
.pack_spot_key <- function(x, y) {
  xi <- as.integer(x)
  yi <- as.integer(y)
  xmin <- min(xi); ymin <- min(yi)
  xi <- xi - xmin
  yi <- yi - ymin
  ymax <- max(yi)
  shift <- max(1L, as.integer(ceiling(log2(ymax + 2))))
  keys <- bitwShiftL(xi, shift) + yi
  list(keys = keys, shift = shift, xmin = xmin, ymin = ymin)
}

# Recover (x, y) coordinates from the bit-packed keys.
.unpack_spot_keys <- function(keys, shift, xmin, ymin) {
  y_mask <- bitwShiftL(1L, shift) - 1L
  y <- bitwAnd(keys, y_mask) + ymin
  x <- bitwShiftR(keys, shift) + xmin
  cbind(x = x, y = y)
}

# Read a square-bin GEF and build a gene x spot sparse matrix wrapped in a
# Seurat object. Internal; called by LoadStereoSeqGef().
.load_stereoseq_squarebin <- function(h5, bin_size, assay, verbose) {
  bin_path <- sprintf("geneExp/bin%d", bin_size)
  if (!h5$exists(bin_path)) {
    available <- if (h5$exists("geneExp")) names(h5[["geneExp"]]) else character(0)
    stop(sprintf(
      "Stereo-seq GEF does not contain bin%d. Available bins: %s",
      bin_size,
      if (length(available) == 0) "<none>" else paste(available, collapse = ", ")),
      call. = FALSE)
  }

  if (verbose) message("  Reading gene table at ", bin_path, "/gene")
  genes_ds <- h5[[paste0(bin_path, "/gene")]]
  genes <- genes_ds$read()
  genes_ds$close()

  if (verbose) message("  Reading expression table at ", bin_path, "/expression")
  expr_ds <- h5[[paste0(bin_path, "/expression")]]
  expr <- expr_ds$read()
  expr_ds$close()

  # Compound string fields may arrive as raw bytes; normalise.
  gene_names <- .read_compound_string(genes$geneName)
  gene_ids   <- .read_compound_string(genes$geneID)

  .validate_gef_consistency(genes, expr, where = bin_path)

  n_genes <- nrow(genes)
  n_expr  <- nrow(expr)

  # Expand the per-gene (offset, count) run into a per-record gene index.
  gene_idx <- rep.int(seq_len(n_genes), times = as.integer(genes$count))

  # Bit-pack spot coordinates to avoid building a 10M-element character vector.
  packed <- .pack_spot_key(expr$x, expr$y)
  spot_fac <- factor(packed$keys)
  spot_idx <- as.integer(spot_fac)
  n_spots  <- nlevels(spot_fac)

  spot_coords <- .unpack_spot_keys(
    as.integer(levels(spot_fac)),
    shift = packed$shift, xmin = packed$xmin, ymin = packed$ymin
  )
  spot_names <- sprintf("spot_%d_%d", spot_coords[, "x"], spot_coords[, "y"])

  # Dedup gene names (Stereo-seq can have multi-mapped probes sharing a name).
  if (anyDuplicated(gene_names)) {
    gene_names <- make.unique(gene_names, sep = ".")
  }

  mat <- sparseMatrix(
    i = gene_idx,
    j = spot_idx,
    x = as.numeric(expr$count),
    dims = c(n_genes, n_spots),
    dimnames = list(gene_names, spot_names)
  )

  seurat_obj <- suppressWarnings(
    CreateSeuratObject(counts = mat, assay = assay, project = "stereoseq")
  )

  seurat_obj@meta.data$spatial_x <- spot_coords[, "x"]
  seurat_obj@meta.data$spatial_y <- spot_coords[, "y"]
  seurat_obj@misc$spatial_technology <- "StereoSeq"

  # Root-level provenance metadata for downstream consumers.
  root_attrs <- list()
  for (name in c("sn", "resolution", "offsetX", "offsetY", "version",
                 "gef_area", "bin_type", "omics")) {
    if (name %in% hdf5r::list.attributes(h5)) {
      val <- tryCatch(h5attr(h5, name), error = function(e) NULL)
      if (!is.null(val)) root_attrs[[name]] <- val
    }
  }
  seurat_obj@misc$stereo_seq <- list(
    bin_size    = bin_size,
    n_genes     = n_genes,
    n_spots     = n_spots,
    gene_ids    = gene_ids,
    coords_xmin = packed$xmin,
    coords_ymin = packed$ymin,
    root_attrs  = root_attrs
  )

  if (verbose) {
    message(sprintf("  Built %d genes x %d spots sparse matrix (%d nonzero)",
                    n_genes, n_spots, length(mat@x)))
  }
  seurat_obj
}

# Read a .cellbin.gef: per-cell metadata + flat cellExp stream.
# Cell exp records are indexed by cell$offset/cell$geneCount.
.load_stereoseq_cellbin <- function(h5, assay, verbose) {
  if (!h5$exists("cellBin/cell") || !h5$exists("cellBin/gene") ||
      !h5$exists("cellBin/cellExp")) {
    stop("Stereo-seq cellbin.gef missing required groups under /cellBin/",
         call. = FALSE)
  }

  if (verbose) message("  Reading cellBin/cell")
  cells_ds <- h5[["cellBin/cell"]]
  cells <- cells_ds$read()
  cells_ds$close()

  if (verbose) message("  Reading cellBin/gene")
  genes_ds <- h5[["cellBin/gene"]]
  genes <- genes_ds$read()
  genes_ds$close()

  if (verbose) message("  Reading cellBin/cellExp")
  expr_ds <- h5[["cellBin/cellExp"]]
  expr <- expr_ds$read()
  expr_ds$close()

  gene_names <- .read_compound_string(genes$geneName)
  gene_ids   <- .read_compound_string(genes$geneID)

  n_genes <- nrow(genes)
  n_cells <- nrow(cells)
  n_expr  <- nrow(expr)

  expected <- sum(as.numeric(cells$geneCount))
  if (!isTRUE(expected == n_expr)) {
    stop(sprintf(
      "Corrupt .cellbin.gef: cell$geneCount sums to %.0f but cellExp has %.0f rows",
      expected, n_expr), call. = FALSE)
  }

  # Expand per-cell (offset, geneCount) run into a record->cell index.
  cell_idx <- rep.int(seq_len(n_cells), times = as.integer(cells$geneCount))
  gene_idx <- as.integer(expr$geneID) + 1L

  if (anyDuplicated(gene_names)) {
    gene_names <- make.unique(gene_names, sep = ".")
  }

  cell_names <- sprintf("cell_%d", as.integer(cells$id))

  mat <- sparseMatrix(
    i = gene_idx,
    j = cell_idx,
    x = as.numeric(expr$count),
    dims = c(n_genes, n_cells),
    dimnames = list(gene_names, cell_names)
  )

  seurat_obj <- suppressWarnings(
    CreateSeuratObject(counts = mat, assay = assay, project = "stereoseq_cellbin")
  )

  seurat_obj@meta.data$spatial_x <- as.numeric(cells$x)
  seurat_obj@meta.data$spatial_y <- as.numeric(cells$y)
  if (!is.null(cells$area)) seurat_obj@meta.data$cell_area <- as.numeric(cells$area)
  if (!is.null(cells$geneCount)) seurat_obj@meta.data$n_genes_gef <- as.integer(cells$geneCount)
  if (!is.null(cells$expCount)) seurat_obj@meta.data$n_umi_gef <- as.integer(cells$expCount)

  seurat_obj@misc$spatial_technology <- "StereoSeqCellBin"
  seurat_obj@misc$stereo_seq <- list(
    n_genes  = n_genes,
    n_cells  = n_cells,
    gene_ids = gene_ids
  )

  if (verbose) {
    message(sprintf("  Built %d genes x %d cells sparse matrix (%d nonzero)",
                    n_genes, n_cells, length(mat@x)))
  }
  seurat_obj
}

#' Load a Stereo-seq GEF file into a Seurat object
#'
#' Pure-R native reader for BGI Stereo-seq \code{.gef} and \code{.cellbin.gef}
#' HDF5 files. No Python or stereopy dependency.
#'
#' The function auto-detects square-bin versus cell-bin layouts from the root
#' \code{bin_type} attribute. Square-bin files are read at the requested bin
#' size (typically 50 or 100 for downstream analysis; bin1 is raw DNB-level
#' data and can be very large). Cell-bin files are read into a gene x cell
#' sparse matrix using the cell segmentation table.
#'
#' @param file Path to a \code{.gef} or \code{.cellbin.gef} file
#' @param bin_size Bin size to read for square-bin files (default 50, ignored
#'   for cellbin)
#' @param assay Seurat assay name (default "Spatial")
#' @param verbose Print progress messages
#'
#' @return A Seurat object with spatial coordinates in
#'   \code{meta.data$spatial_x/y} and provenance in \code{@misc$stereo_seq}
#'
#' @examples
#' \dontrun{
#' seurat_obj <- LoadStereoSeqGef("mouse_embryo.gef", bin_size = 50)
#' scConvert(seurat_obj, "mouse_embryo.h5ad")
#' }
#'
#' @export
LoadStereoSeqGef <- function(file, bin_size = 50, assay = "Spatial",
                             verbose = TRUE) {
  if (!file.exists(file)) {
    stop("Stereo-seq file not found: ", file, call. = FALSE)
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("LoadStereoSeqGef requires the 'hdf5r' package", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("LoadStereoSeqGef requires the 'Matrix' package", call. = FALSE)
  }

  h5 <- H5File$new(file, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  bin_type <- .detect_gef_bin_type(h5)

  if (verbose) {
    message("Loading Stereo-seq GEF: ", basename(file),
            " (bin_type=", bin_type, ")")
  }

  result <- if (bin_type == "cellbin") {
    .load_stereoseq_cellbin(h5, assay = assay, verbose = verbose)
  } else {
    .load_stereoseq_squarebin(h5, bin_size = bin_size,
                              assay = assay, verbose = verbose)
  }

  result@misc$stereo_seq$source_file <- normalizePath(file, mustWork = FALSE)
  result
}
