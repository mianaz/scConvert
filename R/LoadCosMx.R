#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CosMx SMI (NanoString) flat-file bundle reader
#
# Wraps Seurat's LoadNanostring() which already parses the canonical bundle
# (exprMat / metadata / fov_positions / tx_file / polygons CSVs). Adds bundle
# validation, spatial_technology tagging, and scConvert format registry glue.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Canonical substrings found in NanoString CosMx CSV filenames.
.cosmx_canonical_files <- c(
  "exprMat_file",
  "metadata_file",
  "fov_positions_file",
  "tx_file"
)

# Predicate: does `path` look like a CosMx bundle directory? True if it is a
# directory containing at least two canonical CosMx CSV stems. Matches the
# detection logic in src/main.c::is_cosmx_dir() so CLI and R agree.
#
#' @keywords internal
.is_cosmx_dir <- function(path) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) return(FALSE)
  if (!dir.exists(path)) return(FALSE)
  files <- list.files(path, pattern = "\\.(csv|csv\\.gz)$",
                      full.names = FALSE, ignore.case = TRUE)
  if (length(files) == 0) return(FALSE)
  hits <- vapply(.cosmx_canonical_files, function(needle) {
    any(grepl(needle, files, fixed = TRUE))
  }, logical(1))
  sum(hits) >= 2L
}

#' Load a NanoString CosMx SMI bundle into a Seurat object
#'
#' Thin R-only wrapper around \code{Seurat::LoadNanostring()} that validates
#' the bundle layout, sets the scConvert spatial-technology tag, and returns
#' a Seurat object ready for conversion to h5ad / h5Seurat / Zarr via
#' \code{\link{scConvert}}.
#'
#' The bundle directory must contain the canonical NanoString CosMx CSV
#' files (\code{*exprMat_file*.csv}, \code{*metadata_file*.csv},
#' \code{*fov_positions_file*.csv}) and optionally \code{*tx_file*.csv} and
#' polygon files. These are the files produced by NanoString's AtoMx export.
#'
#' @param data.dir Path to the CosMx bundle directory
#' @param fov FOV (field-of-view) name attached to the resulting Seurat object
#'   (default "cosmx")
#' @param assay Seurat assay name (default "Nanostring")
#' @param verbose Print progress messages
#'
#' @return A Seurat object with an FOV slot containing cell centroids,
#'   segmentation boundaries, and transcript molecules
#'
#' @examples
#' \dontrun{
#' seurat_obj <- LoadCosMx("/path/to/cosmx_nsclc/")
#' scConvert(seurat_obj, "cosmx_nsclc.h5ad")
#' }
#'
#' @export
LoadCosMx <- function(data.dir, fov = "cosmx", assay = "Nanostring",
                      verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("LoadCosMx requires the 'Seurat' package (>= 5.0)", call. = FALSE)
  }
  if (!dir.exists(data.dir)) {
    stop("CosMx bundle directory not found: ", data.dir, call. = FALSE)
  }
  if (!.is_cosmx_dir(data.dir)) {
    stop("Directory does not look like a CosMx bundle: ", data.dir,
         "\n  Expected CSV files matching at least two of: ",
         paste(.cosmx_canonical_files, collapse = ", "),
         call. = FALSE)
  }

  if (verbose) {
    message("Loading CosMx bundle via Seurat::LoadNanostring: ", data.dir)
  }

  # Seurat::LoadNanostring hardcodes type = c("centroids", "segmentations")
  # and errors on bundles without a *-polygons.csv file. We inline its body
  # so we can request centroids-only when no polygons are staged.
  has_polygons <- length(list.files(data.dir, pattern = "-polygons\\.csv$",
                                    recursive = FALSE)) > 0
  req_type <- if (has_polygons) c("centroids", "segmentations") else "centroids"
  fov_type <- if (has_polygons) c("segmentation", "centroids") else "centroids"

  data <- Seurat::ReadNanostring(data.dir = data.dir, type = req_type)
  cents <- SeuratObject::CreateCentroids(data$centroids)
  if (has_polygons) {
    segs <- SeuratObject::CreateSegmentation(data$segmentations)
    seg_data <- list(centroids = cents, segmentation = segs)
  } else {
    seg_data <- list(centroids = cents)
  }
  coords <- SeuratObject::CreateFOV(coords = seg_data, type = fov_type,
                                    molecules = data$pixels, assay = assay)
  seurat_obj <- Seurat::CreateSeuratObject(counts = data$matrix, assay = assay)
  cells <- Reduce(intersect, lapply(fov_type, function(bt) {
    SeuratObject::Cells(x = coords, boundary = bt)
  }))
  cells <- intersect(SeuratObject::Cells(seurat_obj), cells)
  coords <- subset(x = coords, cells = cells)
  seurat_obj[[fov]] <- coords

  seurat_obj@misc$spatial_technology <- "CosMx"
  seurat_obj@misc$cosmx <- list(
    bundle_dir  = normalizePath(data.dir, mustWork = FALSE),
    fov_name    = fov,
    assay_name  = assay
  )

  if (verbose) {
    message(sprintf("  CosMx: %d features x %d cells, assay '%s'",
                    nrow(seurat_obj), ncol(seurat_obj), assay))
  }
  seurat_obj
}
