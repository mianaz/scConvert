#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 10x Xenium in-situ bundle reader
#
# Wraps Seurat::LoadXenium() which already parses the canonical Xenium output
# bundle (cells.csv.gz, cell_feature_matrix.h5, transcripts.parquet,
# nucleus_boundaries.parquet, cell_boundaries.parquet). Adds bundle validation,
# spatial_technology tagging, and scConvert format-registry glue so the result
# flows cleanly into WriteH5Seurat / WriteH5AD / WriteZarr.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.xenium_canonical_files <- c(
  "cell_feature_matrix",
  "cells.csv",
  "cells.parquet",
  "transcripts.parquet",
  "cell_boundaries.parquet",
  "nucleus_boundaries.parquet"
)

#' @keywords internal
.is_xenium_dir <- function(path) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) return(FALSE)
  if (!dir.exists(path)) return(FALSE)
  files <- list.files(path, recursive = FALSE, full.names = FALSE)
  if (length(files) == 0) return(FALSE)
  hits <- vapply(.xenium_canonical_files, function(needle) {
    any(grepl(needle, files, fixed = TRUE))
  }, logical(1))
  sum(hits) >= 2L
}

#' Load a 10x Xenium output bundle into a Seurat object
#'
#' Thin R-only wrapper around \code{Seurat::LoadXenium()} that validates the
#' bundle layout, sets the scConvert spatial-technology tag, and returns a
#' Seurat object ready for conversion to h5ad / h5Seurat / Zarr via
#' \code{\link{scConvert}}.
#'
#' The bundle directory must contain the canonical 10x Xenium output files
#' (\code{cell_feature_matrix.h5} or \code{cell_feature_matrix/}, \code{cells.csv.gz}
#' or \code{cells.parquet}, and optional \code{transcripts.parquet} /
#' \code{cell_boundaries.parquet} / \code{nucleus_boundaries.parquet}).
#'
#' @param data.dir Path to the Xenium output directory
#' @param fov FOV name attached to the resulting Seurat object (default "fov")
#' @param assay Seurat assay name (default "Xenium")
#' @param verbose Print progress messages
#'
#' @return A Seurat object with an FOV slot containing cell centroids,
#'   segmentation boundaries, and transcript molecules
#'
#' @examples
#' \dontrun{
#' seurat_obj <- LoadXenium("/path/to/xenium_lung/")
#' scConvert(seurat_obj, "xenium_lung.h5ad")
#' }
#'
#' @export
LoadXenium <- function(data.dir, fov = "fov", assay = "Xenium",
                       verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("LoadXenium requires the 'Seurat' package (>= 5.0)", call. = FALSE)
  }
  if (!exists("LoadXenium", envir = asNamespace("Seurat"), inherits = FALSE)) {
    stop("Seurat::LoadXenium is unavailable; upgrade to Seurat >= 5.0",
         call. = FALSE)
  }
  if (!dir.exists(data.dir)) {
    stop("Xenium bundle directory not found: ", data.dir, call. = FALSE)
  }
  if (!.is_xenium_dir(data.dir)) {
    stop("Directory does not look like a Xenium bundle: ", data.dir,
         "\n  Expected at least two of: ",
         paste(.xenium_canonical_files, collapse = ", "),
         call. = FALSE)
  }

  if (verbose) {
    message("Loading Xenium bundle via Seurat::LoadXenium: ", data.dir)
  }

  seurat_obj <- Seurat::LoadXenium(data.dir = data.dir, fov = fov,
                                   assay = assay)

  seurat_obj@misc$spatial_technology <- "Xenium"
  seurat_obj@misc$xenium <- list(
    bundle_dir = normalizePath(data.dir, mustWork = FALSE),
    fov_name   = fov,
    assay_name = assay
  )

  if (verbose) {
    message(sprintf("  Xenium: %d features x %d cells, assay '%s'",
                    nrow(seurat_obj), ncol(seurat_obj), assay))
  }
  seurat_obj
}
