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
#' The bundle directory must contain the canonical NanoString CosMx CSV files
#' (\code{*exprMat_file*.csv}, \code{*metadata_file*.csv},
#' \code{*fov_positions_file*.csv}) and optionally \code{*tx_file*.csv} and
#' polygon files. These are the files produced by NanoString's AtoMx export.
#'
#' The \code{*metadata_file*.csv} is read and attached to the Seurat object
#' after construction. This populates cell morphology columns (\code{Area},
#' \code{AspectRatio}, \code{Width}, \code{Height},
#' \code{CenterX_global_px}, \code{CenterY_global_px}) and co-registered
#' immunofluorescence intensity columns (\code{Mean.PanCK}, \code{Mean.CD45},
#' \code{Mean.CD3}, \code{Mean.MembraneStain}, etc.) that are not returned by
#' \code{Seurat::ReadNanostring()} alone.
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

  # Attach per-cell metadata from metadata_file.csv (Area, AspectRatio,
  # CenterX/Y_global_px, Width, Height, immunofluorescence Mean.*/Max.*
  # columns). ReadNanostring only returns the count matrix and centroids;
  # the metadata file is the only source of these columns.
  meta_files <- list.files(data.dir, pattern = "metadata_file",
                           full.names = TRUE, ignore.case = TRUE)
  meta_files <- meta_files[grepl("\\.csv(\\.gz)?$", meta_files,
                                 ignore.case = TRUE)]
  if (length(meta_files) > 0L) {
    meta_df <- tryCatch(
      read.csv(meta_files[[1L]], check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(meta_df) &&
        all(c("cell_ID", "fov") %in% colnames(meta_df))) {
      rownames(meta_df) <- paste0(meta_df[["cell_ID"]], "_",
                                  meta_df[["fov"]])
      keep <- intersect(rownames(meta_df), SeuratObject::Cells(seurat_obj))
      meta_df <- meta_df[keep, , drop = FALSE]
      if (nrow(meta_df) > 0L) {
        seurat_obj <- SeuratObject::AddMetaData(seurat_obj, meta_df)
        if (verbose) {
          message(sprintf("  Attached %d metadata columns from metadata_file",
                          ncol(meta_df)))
        }
      }
    }
  }

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
