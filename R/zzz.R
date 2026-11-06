#' @importFrom rlang %||%
#' @importFrom hdf5r H5T_COMPOUND
#' @importFrom methods setOldClass is
#' @importFrom stats na.omit
#' @importFrom utils data read.csv
#' @useDynLib scConvert, .registration = TRUE
#'
NULL

#' @docType package
#' @name scConvert-package
#' @rdname scConvert-package
#'
#' @description
#' Universal single-cell format conversion for R. Converts between h5ad
#' (AnnData), h5Seurat, h5mu (MuData), Loom, Zarr, RDS, and
#' SingleCellExperiment formats through a hub architecture using Seurat as
#' an intermediate. No Python dependency required.
#'
#' @section Package options:
#' Set with \code{\link[base]{options}}:
#' \describe{
#'  \item{\code{scConvert.dtypes.logical_to_int}}{Store logical vectors as
#'   integers (0/1, NA = 2) so other languages can read them (default
#'   \code{TRUE})}
#'  \item{\code{scConvert.dtypes.dataframe_as_group}}{Store every data frame
#'   as a group of column datasets, also when it has no factor columns
#'   (default \code{TRUE})}
#'  \item{\code{scConvert.chunking.MARGIN}}{Axis used when streaming a dataset
#'   in blocks: \code{"largest"} (default), \code{"smallest"},
#'   \code{"first"} or \code{"last"}}
#'  \item{\code{scConvert.dimreducs.allglobal}}{Load every dimensional
#'   reduction as global (default \code{FALSE})}
#'  \item{\code{scConvert.compression.level}}{gzip level (0-9) for new HDF5
#'   datasets; 0 disables compression (default 4)}
#'  \item{\code{scConvert.use_cli}}{Use the optional compiled \code{scconvert}
#'   binary for file-to-file conversions (default \code{FALSE}). It is fast
#'   but drops feature metadata and spatial images; for large files prefer
#'   the streaming R path (\code{scConvert(..., stream = TRUE)})}
#'  \item{\code{scConvert.stream}}{Default for the \code{stream} argument of
#'   \code{\link{scConvert}}: \code{"auto"} (stream files larger than
#'   \code{scConvert.stream_threshold} bytes, 2 GB by default), \code{TRUE} or
#'   \code{FALSE}}
#'  \item{\code{scConvert.stream_chunk_bytes}}{Memory budget per block when
#'   streaming matrices (default 64 MB)}
#' }
#'
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default.options <- list(
  "scConvert.dtypes.logical_to_int" = TRUE,
  "scConvert.dtypes.dataframe_as_group" = TRUE,
  "scConvert.chunking.MARGIN" = c("largest", "smallest", "first", "last"),
  "scConvert.dimreducs.allglobal" = FALSE,
  "scConvert.compression.level" = 4L,
  "scConvert.use_cli" = FALSE
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# hdf5r open modes that create a file / open an existing one
modes <- list(new = c('w', 'w-', 'x'), existing = c('r', 'r+'))

# x.y.z, optionally with a .9xxx development suffix
version.regex <- '^\\d+(\\.\\d+){2}(\\.9\\d{3})?$'

# Registry of scdisk format classes (see RegisterSCDisk)
scdisk.types <- new.env()

# h5Seurat layout versions: spatial images, and Seurat v5 assays
spatial.version <- '3.1.5.9900'
v5.version <- '5.0.0'

# Standard Seurat assay layer names
.STANDARD_LAYERS <- c("counts", "data", "scale.data")

# Threshold for determining if matrix should be stored as sparse (90% zeros)
.SPARSITY_THRESHOLD <- 0.9

# AnnData encoding defaults
.ANNDATA_ENCODING_VERSION <- "0.2.0"

# Cached HDF5 objects (populated in .onLoad)
.h5_scalar_space <- NULL
.h5_utf8_dtype <- NULL
.h5_guessdtype_cache <- new.env(parent = emptyenv())

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal utility functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert BPCells matrix objects to dgCMatrix
#'
#' @param mat Matrix object (potentially BPCells IterableMatrix or RenameDims)
#' @param verbose Show conversion message
#' @return dgCMatrix or original matrix if not BPCells
#' @keywords internal
#'
ConvertBPCellsMatrix <- function(mat, verbose = FALSE) {
  if (is.null(mat)) return(NULL)

  # Check if it's a BPCells object
  is_bpcells <- inherits(mat, c("IterableMatrix", "RenameDims")) ||
    (is.object(mat) && !inherits(mat, c("dgCMatrix", "dgTMatrix", "matrix")))

  if (!is_bpcells) return(mat)

  if (verbose) message("Converting BPCells object to dgCMatrix")

  tryCatch(
    as(mat, "dgCMatrix"),
    error = function(e) {
      warning("BPCells conversion failed: ", e$message, call. = FALSE)
      mat
    }
  )
}


#' File type of a path or URI
#'
#' The lower-cased extension of \code{file} (or the whole name when it has
#' none, so a bare type such as \code{"h5ad"} maps to itself). URI schemes
#' (\code{soma://...}), \code{.spatialdata.zarr}, Stereo-seq GEF files and
#' CosMx export directories are recognised specially.
#'
#' @param file Path, URI or format name
#'
#' @return A lower-case type string such as \code{"h5ad"}
#' @keywords internal
#' @importFrom tools file_ext
#'
#' @examples
#' \donttest{
#' scConvert:::FileType("pbmc3k.h5Seurat")  # "h5seurat"
#' scConvert:::FileType("h5ad")             # "h5ad"
#' }
#'
FileType <- function(file) {
  # URI schemes (e.g. soma://collection/measurement)
  if (grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", file)) {
    return(tolower(sub("://.*", "", file)))
  }
  # Multi-part extensions (e.g. .spatialdata.zarr)
  bn <- basename(file)
  if (grepl("\\.spatialdata\\.zarr$", bn, ignore.case = TRUE)) {
    return("spatialdata.zarr")
  }
  # Stereo-seq cellbin GEF must be checked before generic .gef.
  if (grepl("\\.cellbin\\.gef$", bn, ignore.case = TRUE)) {
    return("cellbin.gef")
  }
  if (grepl("\\.gef$", bn, ignore.case = TRUE)) {
    return("gef")
  }
  # CosMx SMI bundles are distributed as a directory of CSV files, not a
  # single-extension container. Delegate detection to .is_cosmx_dir().
  if (dir.exists(file) && exists(".is_cosmx_dir", mode = "function") &&
      isTRUE(.is_cosmx_dir(file))) {
    return("cosmx")
  }
  ext <- file_ext(x = file)
  ext <- ifelse(test = nchar(x = ext), yes = ext, no = bn)
  return(tolower(x = ext))
}


#' Add AnnData encoding attributes to an HDF5 dataset or group
#'
#' Helper function to add standard AnnData encoding-type and encoding-version
#' attributes to an HDF5 object. This is used for string arrays, categorical
#' variables, and dataframes to ensure compatibility with AnnData/scanpy.
#'
#' @param h5obj An hdf5r H5D or H5Group object
#' @param encoding_type The encoding type (e.g., 'string-array', 'categorical', 'dataframe')
#' @param encoding_version The encoding version (default: '0.2.0')
#'
#' @return NULL (modifies h5obj in place)
#'
#' @keywords internal
#'
AddAnndataEncoding <- function(h5obj, encoding_type = 'string-array', encoding_version = '0.2.0') {
  h5obj$create_attr(
    attr_name = 'encoding-type',
    robj = encoding_type,
    dtype = CachedGuessDType(x = encoding_type),
    space = ScalarSpace()
  )
  h5obj$create_attr(
    attr_name = 'encoding-version',
    robj = encoding_version,
    dtype = CachedGuessDType(x = encoding_version),
    space = ScalarSpace()
  )
  invisible(NULL)
}


#' Get compression level for HDF5 dataset creation
#'
#' @return Integer gzip compression level (0-9)
#'
#' @keywords internal
#'
GetCompressionLevel <- function() {
  level <- getOption("scConvert.compression.level", default = 4L)
  level <- as.integer(level)
  if (is.na(level) || level < 0L || level > 9L) {
    warning("Invalid compression level, using default (4)", immediate. = TRUE)
    level <- 4L
  }
  return(level)
}

#' Get cached scalar HDF5 dataspace
#'
#' Returns a singleton H5S scalar dataspace, avoiding repeated allocation.
#'
#' @return An \code{\link[hdf5r]{H5S}} object of type "scalar"
#'
#' @keywords internal
#'
ScalarSpace <- function() {
  if (is.null(.h5_scalar_space)) {
    .h5_scalar_space <<- H5S$new(type = "scalar")
  }
  .h5_scalar_space
}

#' Get cached UTF-8 string dtype
#'
#' Returns a singleton UTF-8 string type, avoiding repeated allocation.
#'
#' @return An \code{\link[hdf5r]{H5T_STRING}} object
#'
#' @keywords internal
#'
CachedUtf8Type <- function() {
  if (is.null(.h5_utf8_dtype)) {
    .h5_utf8_dtype <<- StringType('utf8')
  }
  .h5_utf8_dtype
}

#' Cached version of GuessDType for known string constants
#'
#' Checks an internal cache before calling \code{GuessDType()}.
#' This avoids repeated H5T allocation for the same AnnData encoding strings.
#'
#' @param x An R object to guess the HDF5 dtype for
#' @param stype String encoding type (default 'utf8')
#'
#' @return An object of class \code{\link[hdf5r]{H5T}}
#'
#' @keywords internal
#'
CachedGuessDType <- function(x, stype = 'utf8') {
  if (is.character(x) && length(x) == 1L) {
    cached <- .h5_guessdtype_cache[[x]]
    if (!is.null(cached)) return(cached)
    dtype <- GuessDType(x = x, stype = stype)
    .h5_guessdtype_cache[[x]] <- dtype
    return(dtype)
  }
  GuessDType(x = x, stype = stype)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# V5 Safe Helper Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Helper: create CLI-first direct path with R streaming fallback
# Defined outside .onLoad to avoid R CMD check NOTE about message() in startup

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Format converter functions (defined outside .onLoad to avoid R CMD check NOTE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.h5seurat_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  readH5Seurat(file = file, verbose = verbose, ...)
}

.h5seurat_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  writeH5Seurat(object = object, filename = filename,
               overwrite = overwrite, verbose = verbose, ...)
  invisible(filename)
}

.h5ad_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  readH5AD(file = file, assay.name = assay, verbose = verbose, ...)
}

.h5ad_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE,
                        standardize = FALSE, ...) {
  assay_name <- DefaultAssay(object = object)
  if (isTRUE(standardize)) {
    # the C writer has no renaming step, so rename up front
    object@meta.data <- .ScanpyMetaNames(object@meta.data)
    standardize <- FALSE
  }

  # Detect BPCells on-disk layers to avoid materializing them in RAM
  bpcells_layers <- list()
  if (requireNamespace("BPCells", quietly = TRUE)) {
    for (ln in SeuratObject::Layers(object[[assay_name]])) {
      dat <- tryCatch(
        SeuratObject::GetAssayData(object, assay = assay_name, layer = ln),
        error = function(e) NULL
      )
      if (!is.null(dat) && inherits(dat, c("IterableMatrix", "RenameDims"))) {
        bpcells_layers[[ln]] <- dat
      }
    }
  }

  if (length(bpcells_layers) > 0) {
    if (verbose) message("BPCells on-disk matrices detected: streaming to h5ad (no full materialization)")

    # Replace BPCells layers with zero-entry placeholders (correct dims, ~zero RAM)
    obj_write <- object
    for (ln in names(bpcells_layers)) {
      bp_mat <- bpcells_layers[[ln]]
      nr <- nrow(bp_mat); nc <- ncol(bp_mat)
      placeholder <- new("dgCMatrix",
                         i = integer(0), p = rep(0L, nc + 1L), x = numeric(0),
                         Dim = c(as.integer(nr), as.integer(nc)),
                         Dimnames = list(rownames(bp_mat), colnames(bp_mat)))
      obj_write <- SeuratObject::SetAssayData(
        obj_write, assay = assay_name, layer = ln, new.data = placeholder
      )
    }

    # Write everything except the on-disk matrices (placeholders), then
    # stream those matrices into the file below
    DirectSeuratToH5AD(object = obj_write, filename = filename,
                       assay = assay_name, overwrite = overwrite,
                       verbose = verbose, standardize = standardize)
    rm(obj_write)

    gzip <- GetCompressionLevel()
    has_data_layer <- "data" %in% names(bpcells_layers)
    for (ln in names(bpcells_layers)) {
      if (ln == "data") {
        group_name <- "X"
      } else if (ln == "counts" && has_data_layer) {
        group_name <- "layers/counts"
      } else if (ln == "counts") {
        group_name <- "X"
      } else {
        group_name <- paste0("layers/", ln)
      }

      h5 <- hdf5r::H5File$new(filename, mode = "r+")
      if (h5$exists(group_name)) h5$link_delete(group_name)
      h5$close_all()

      get("write_matrix_anndata_hdf5", envir = asNamespace("BPCells"))(
        mat = bpcells_layers[[ln]],
        path = filename,
        group = group_name,
        gzip_level = gzip
      )

      if (verbose) message("  Streamed '", ln, "' -> ", group_name, " via BPCells (no materialization)")
    }
  } else {
    if (file.exists(filename) && !overwrite) {
      stop("File '", filename, "' already exists; set overwrite = TRUE", call. = FALSE)
    }
    # Fast path: Direct C h5ad writer (no h5seurat intermediate, no transpose)
    # Enabled by default when compiled code is available. Disable with
    # options(scConvert.use_c_writer = FALSE)
    c_h5ad <- !isFALSE(getOption("scConvert.use_c_writer")) &&
              is.loaded("C_write_h5ad", PACKAGE = "scConvert") &&
              length(Assays(object)) == 1
    if (c_h5ad) {
      c_ok <- tryCatch(
        .writeH5AD_c(object, filename = filename, overwrite = overwrite,
                      verbose = verbose),
        error = function(e) FALSE
      )
      if (isTRUE(c_ok)) {
        # C writer handles matrix/meta/reductions/graphs but not misc/uns.
        # Append misc and spatial data via R if present.
        misc <- tryCatch(Misc(object), error = function(e) list())
        images <- tryCatch(Images(object), error = function(e) character(0))
        skip_keys <- c("__varp__", ".__h5ad_path__", ".__h5ad_loaded__",
                       "scConvert_read")
        skip_keys <- c(skip_keys, grep("^__varp__\\.", names(misc), value = TRUE))
        misc_to_write <- misc[!names(misc) %in% skip_keys]
        misc_to_write <- misc_to_write[!vapply(misc_to_write, is.null, logical(1))]
        has_varp <- !is.null(tryCatch(Misc(object)[["__varp__"]], error = function(e) NULL))
        if (!(length(misc_to_write) > 0 || length(images) > 0 || has_varp)) {
          # Nothing to append; stamp provenance by path.
          .h5ad_stamp_provenance(filename)
        } else {
          h5 <- hdf5r::H5File$new(filename, mode = "r+")
          on.exit(h5$close_all(), add = TRUE)
          if (!h5$exists("uns")) h5$create_group("uns")
          uns_grp <- h5[["uns"]]
          gzip <- GetCompressionLevel()
          for (item_name in names(misc_to_write)) {
            tryCatch(
              WriteUnsItem(uns_grp, item_name, misc_to_write[[item_name]], gzip),
              error = function(e) {
                if (verbose) message("  Could not write uns/", item_name, ": ", e$message)
              }
            )
          }
          if (length(images) > 0) {
            tryCatch(
              SeuratSpatialToH5AD(object, h5, verbose = verbose),
              error = function(e) {
                if (verbose) message("Spatial data conversion failed: ", e$message)
              }
            )
          }
          # Write varp if present
          varp_data <- tryCatch(Misc(object)[["__varp__"]], error = function(e) NULL)
          if (!is.null(varp_data) && is.list(varp_data) && length(varp_data) > 0) {
            if (!h5$exists("varp")) h5$create_group("varp")
            varp_grp <- h5[["varp"]]
            for (vn in names(varp_data)) {
              tryCatch({
                v <- varp_data[[vn]]
                if (inherits(v, c("dgCMatrix", "dgRMatrix", "dgTMatrix"))) {
                  WriteCSRGroup(varp_grp, vn, v, gzip = gzip)
                } else if (is.matrix(v)) {
                  varp_grp$create_dataset(vn, robj = v,
                                          chunk_dims = dim(v), gzip_level = gzip)
                }
              }, error = function(e) NULL)
            }
          }
          # Stamp last so /uns/scConvert reflects the final file (a stale
          # misc$scConvert echoed back through WriteUnsItem is replaced).
          .h5ad_stamp_provenance(h5)
        }
        return(invisible(filename))
      }
      if (verbose) message("C h5ad writer failed, falling back to R writer")
    }
    DirectSeuratToH5AD(object = object, filename = filename,
                       assay = assay_name, overwrite = overwrite,
                       verbose = verbose, standardize = standardize)
  }
  invisible(filename)
}

.loom_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  readLoom(file = file, verbose = verbose, ...)
}

.loom_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  writeLoom(object = object, filename = filename,
           overwrite = overwrite, verbose = verbose, ...)
  invisible(filename)
}

.h5mu_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  readH5MU(file = file, verbose = verbose)
}

.h5mu_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  if (file.exists(filename) && !overwrite) {
    stop("File '", filename, "' already exists; set overwrite = TRUE", call. = FALSE)
  }
  writeH5MU(object = object, filename = filename,
                   overwrite = overwrite, verbose = verbose)
  invisible(filename)
}

.rds_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  obj <- readRDS(file = file)
  if (!inherits(x = obj, what = 'Seurat')) {
    obj <- tryCatch(
      expr = as.Seurat(x = obj, verbose = verbose, ...),
      error = function(e) {
        stop("RDS file does not contain a Seurat-coercible object: ",
             conditionMessage(e), call. = FALSE)
      }
    )
  }
  obj
}

.rds_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  if (file.exists(filename) && !overwrite) {
    stop("Destination RDS file exists", call. = FALSE)
  }
  saveRDS(object = object, file = filename)
  if (verbose) message("Saved Seurat object to ", filename)
  invisible(filename)
}

.zarr_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  readZarr(file = file, assay.name = assay, verbose = verbose, ...)
}

.zarr_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  writeZarr(object = object, filename = filename,
           overwrite = overwrite, verbose = verbose, ...)
  invisible(filename)
}


# File-to-file conversion between two HDF5 formats. `source` may be a path
# or an open connection. Arguments a path does not use (e.g. `standardize`
# for non-h5ad targets) are ignored.
.make_direct <- function(stype, dtype) {
  force(stype)
  force(dtype)
  function(source, dest, assay = 'RNA', overwrite = FALSE, verbose = TRUE,
           standardize = FALSE, stream = NULL, ...) {
    src_path <- if (inherits(source, c('H5File', 'h5Seurat', 'loom'))) source$filename else source
    if (isTRUE(getOption("scConvert.use_cli")) && !isTRUE(standardize)) {
      ok <- tryCatch(scConvert_cli(input = src_path, output = dest, assay = assay %||% 'RNA',
                                   overwrite = overwrite, verbose = verbose),
                     error = function(e) FALSE)
      if (isTRUE(ok)) return(invisible(dest))
      if (verbose) message("CLI conversion failed, using the R conversion")
    }
    if (.use_streaming(stream, src_path)) {
      return(.stream_convert(src_path, dest, stype, dtype, assay = assay, overwrite = overwrite,
                             verbose = verbose, standardize = standardize))
    }
    extra <- if (dtype == 'h5ad') list(standardize = standardize) else list()
    do.call(.hub_from_connection, c(list(source, dest, stype = stype, dtype = dtype,
                                         assay = assay %||% 'RNA', overwrite = overwrite,
                                         verbose = verbose), extra))
  }
}

# Without the CLI, file-to-file conversion goes through an in-memory Seurat
# object. `source` may be a path or an open H5File / h5Seurat connection.
.hub_from_connection <- function(source, dest, stype, dtype, assay = "RNA",
                                 overwrite = FALSE, verbose = TRUE, ...) {
  if (file.exists(dest) && !overwrite) {
    stop("Destination file '", dest, "' already exists; set overwrite = TRUE",
         call. = FALSE)
  }
  obj <- if (inherits(source, "h5Seurat")) {
    readH5Seurat(file = source, verbose = verbose)
  } else {
    path <- if (inherits(source, "H5File")) source$filename else source
    GetLoader(ext = stype)(file = path, assay = assay, verbose = verbose)
  }
  if (stype == "h5seurat" && !is.null(assay) && assay %in% Assays(obj)) {
    DefaultAssay(obj) <- assay
  }
  GetSaver(ext = dtype)(object = obj, filename = dest, overwrite = overwrite,
                        verbose = verbose, ...)
  invisible(dest)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Loading handler
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  # let S4 methods dispatch on the R6 connection classes
  setOldClass(Classes = c('scdisk', 'h5Seurat'))
  for (cls in list(h5Seurat, loom)) RegisterSCDisk(r6class = cls)

  # package defaults for options the user has not set
  unset <- setdiff(names(default.options), names(options()))
  if (length(unset)) options(default.options[unset])

  # Pre-populate HDF5 caches for common AnnData encoding strings
  .h5_scalar_space <<- H5S$new(type = "scalar")
  .h5_utf8_dtype <<- StringType('utf8')
  for (s in c('string-array', '0.2.0', '0.1.0', 'categorical', 'array',
              'csr_matrix', 'csc_matrix', 'dataframe', 'numeric-scalar',
              '_index')) {
    .h5_guessdtype_cache[[s]] <<- GuessDType(x = s)
  }

  # The compiled `scconvert` binary is opt-in (options(scConvert.use_cli =
  # TRUE)): it is fast but does not carry everything the R paths do
  # (feature metadata, spatial images); see ?scConvert_cli.

  # Register file format loaders and savers
  RegisterFormat(ext = 'h5seurat', loader = .h5seurat_loader, saver = .h5seurat_saver)
  RegisterFormat(ext = 'h5ad', loader = .h5ad_loader, saver = .h5ad_saver)
  RegisterFormat(ext = 'loom', loader = .loom_loader, saver = .loom_saver)
  RegisterFormat(ext = 'h5mu', loader = .h5mu_loader, saver = .h5mu_saver)
  RegisterFormat(ext = 'rds', loader = .rds_loader, saver = .rds_saver)
  RegisterFormat(ext = 'zarr', loader = .zarr_loader, saver = .zarr_saver)

  # Register direct HDF5-level conversion paths (h5ad <-> h5seurat)
  # File-to-file paths between the HDF5 formats: the opt-in CLI, streaming
  # (stream = TRUE / "auto" for large files) or an in-memory conversion
  hdf5_formats <- c('h5ad', 'h5seurat', 'h5mu', 'loom')
  for (s in hdf5_formats) for (d in setdiff(hdf5_formats, s)) {
    RegisterDirectPath(s, d, .make_direct(s, d))
  }

  # Register SOMA format
  RegisterFormat(ext = 'soma', loader = function(source, ...) {
    readSOMA(source, ...)
  }, saver = function(object, filename, ...) {
    writeSOMA(object, filename, ...)
  })

  # Register SpatialData zarr format (uses .spatialdata.zarr extension)
  RegisterFormat(ext = 'spatialdata.zarr', loader = function(source, ...) {
    readSpatialData(source, ...)
  }, saver = function(object, filename, ...) {
    writeSpatialData(object, filename, ...)
  })

  # Register vendor raw spatial formats (read-only, no corresponding saver).
  # Stereo-seq square-bin and cellbin GEFs and NanoString CosMx flat-file
  # bundles are upstream of the scverse container ecosystem; scConvert loads
  # them into Seurat and then routes through the standard write paths.
  RegisterFormat(ext = 'gef',
    loader = function(source, assay = 'Spatial', verbose = TRUE, ...) {
      bin_size <- list(...)$bin_size
      if (is.null(bin_size)) bin_size <- 50L
      LoadStereoSeqGef(file = source, bin_size = bin_size,
                       assay = assay, verbose = verbose)
    },
    saver = NULL)

  RegisterFormat(ext = 'cellbin.gef',
    loader = function(source, assay = 'Spatial', verbose = TRUE, ...) {
      LoadStereoSeqGef(file = source, assay = assay, verbose = verbose)
    },
    saver = NULL)

  RegisterFormat(ext = 'cosmx',
    loader = function(source, assay = 'Nanostring', verbose = TRUE, ...) {
      fov <- list(...)$fov
      if (is.null(fov)) fov <- 'cosmx'
      LoadCosMx(data.dir = source, fov = fov, assay = assay, verbose = verbose)
    },
    saver = NULL)

  invisible(x = NULL)
}
