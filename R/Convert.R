#' @include zzz.R
#' @include FormatRegistry.R
#' @include Bindings.R
#' @importFrom utils setTxtProgressBar
#' @importFrom hdf5r H5File h5attr H5S
#' @importFrom tools file_path_sans_ext
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert single-cell datasets between formats
#'
#' Universal converter between single-cell file formats and object types.
#' Supports arbitrary source/destination pairs by routing through Seurat as
#' a hub format. Direct HDF5-level paths (h5ad <-> h5seurat) are used when
#' available for memory efficiency.
#'
#' @param source Source dataset: a Seurat object, SingleCellExperiment, loom
#'   connection, filename path, or H5File connection
#' @param dest Name/path of destination file or format. Supported formats:
#'   h5seurat, h5ad, h5mu, loom, rds. Also accepts \code{"sce"} to return
#'   a SingleCellExperiment object (in-memory, no file created).
#' @param assay For h5Seurat -> other formats: name of assay to convert.
#'   For other formats -> h5Seurat: name to assign to the assay.
#'   Default is "RNA".
#' @param overwrite Logical; if \code{TRUE}, overwrite an existing destination file.
#'   Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} (default), show progress updates
#' @param standardize Logical; if \code{TRUE}, convert Seurat-style metadata column names
#'   to scanpy/AnnData conventions when converting to h5ad format. For example,
#'   \code{nCount_RNA} becomes \code{n_counts}, \code{nFeature_RNA} becomes \code{n_genes}.
#'   Only applicable for conversions to h5ad format. Default is \code{FALSE}.
#' @param ... Arguments passed to specific conversion methods. File-to-file
#'   conversion between h5ad, h5Seurat, h5mu and loom files also takes
#'   \code{stream}: \code{TRUE} copies
#'   the expression matrices between the files in bounded-memory blocks instead
#'   of loading the object (see \emph{Large datasets}), \code{FALSE} loads it,
#'   and \code{"auto"} (the default, \code{options(scConvert.stream)}) streams
#'   files larger than \code{options(scConvert.stream_threshold)} (2 GB).
#'
#' @section Large datasets:
#' With \code{stream = TRUE}, only the small parts of the file (cell and
#' feature metadata, embeddings, graphs, images, unstructured data) are read
#' into memory and converted by the regular reader and writer; the
#' expression matrices (h5ad X, layers and raw, h5Seurat layers, h5mu
#' modalities, loom matrix and layers) are copied between the files by HDF5
#' itself, or in blocks of \code{options(scConvert.stream_chunk_bytes)} (64 MB)
#' when they have to be transposed or densified, or features subset or
#' reordered. The result is the same
#' as an in-memory conversion; peak memory depends on the number of cells,
#' not on the size of the matrices. Use \code{\link{scConvertBatch}} to convert
#' many files in parallel.
#'
#' @return For file destinations, invisibly returns the destination filename.
#'   For \code{dest = "sce"}, returns a SingleCellExperiment object.
#'
#'
#' @details
#' \strong{Supported Formats:}
#' \itemize{
#'   \item \strong{R objects}: Seurat, SingleCellExperiment (requires
#'     \pkg{SingleCellExperiment}), loom (R6 connection)
#'   \item \strong{File formats}: h5seurat, h5ad, h5mu, loom, rds
#' }
#'
#' Any source format can be converted to any destination format. Conversions
#' without a direct path go through Seurat as a universal hub:
#' \code{Source -> Seurat -> Destination}.
#'
#' \strong{Direct Paths} (memory-efficient, no full dataset loading):
#' \itemize{
#'   \item \code{h5ad <-> h5seurat}: Direct HDF5-level copy
#' }
#'
#' \strong{Key Features:}
#' \itemize{
#'   \item Preserves expression matrices, metadata, and dimensional reductions
#'   \item For Visium/spatial data: reconstructs images with scale factors
#'   \item Handles multiple data layers (V5 compatibility)
#' }
#'
#' @seealso
#' \code{\link{writeH5AD}} for direct Seurat to h5ad convenience function
#' \code{\link{writeH5Seurat}} to save Seurat objects
#' \code{\link{readH5Seurat}} to load h5Seurat files
#' \code{\link{readH5AD}} to directly load h5ad files
#' \code{\link{readH5MU}} to load h5mu files
#' \code{\link{scConnect}} to establish file connections
#'
#' @examples
#' \dontrun{
#' library(scConvert)
#' library(Seurat)
#'
#' # --- Any format to any format ---
#' scConvert("data.h5ad", dest = "data.h5seurat")     # h5ad -> h5seurat
#' scConvert("data.h5ad", dest = "data.rds")           # h5ad -> RDS
#' scConvert("data.h5ad", dest = "data.loom")           # h5ad -> loom
#' scConvert("data.h5mu", dest = "data.h5ad")           # h5mu -> h5ad
#' scConvert("data.loom", dest = "data.h5seurat")       # loom -> h5seurat
#' scConvert("data.rds",  dest = "data.h5ad")           # RDS -> h5ad
#'
#' # --- From R objects ---
#' scConvert(seurat_obj, dest = "output.h5ad")          # Seurat -> h5ad
#' scConvert(seurat_obj, dest = "output.loom")           # Seurat -> loom
#' scConvert(seurat_obj, dest = "output.rds")            # Seurat -> RDS
#' scConvert(sce_obj, dest = "output.h5ad")              # SCE -> h5ad
#' sce <- scConvert(seurat_obj, dest = "sce")            # Seurat -> SCE
#' }
#'
#' @name scConvert
#' @rdname scConvert
#'
#' @export
#'
scConvert <- function(source, dest, assay, overwrite = FALSE, verbose = TRUE, standardize = FALSE, ...) {
  if (!missing(x = dest) && !is.character(x = dest)) {
    stop("'dest' must be a filename or type", call. = FALSE)
  }
  UseMethod(generic = 'scConvert', object = source)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom hdf5r H5File
#'
#' @rdname scConvert
#' @method scConvert character
#' @export
#'
scConvert.character <- function(
  source,
  dest,
  assay,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  is_uri <- grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", source)
  if (!is_uri && !file.exists(source) && !dir.exists(source)) {
    stop("Source file not found: ", source, call. = FALSE)
  }
  stype <- FileType(file = source)

  # RDS: read the object and dispatch on its class
  if (stype == 'rds') {
    obj <- readRDS(file = source)
    if (!inherits(x = obj, what = 'Seurat')) {
      obj <- tryCatch(
        expr = as.Seurat(x = obj, verbose = verbose),
        error = function(e) {
          stop("RDS file does not contain a Seurat-coercible object: ",
               conditionMessage(e), call. = FALSE)
        }
      )
    }
    if (missing(x = assay)) assay <- DefaultAssay(object = obj)
    return(scConvert(
      source = obj, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, standardize = standardize, ...
    ))
  }

  # Vendor raw spatial formats (.gef, .cellbin.gef, CosMx bundles): invoke the
  # registered loader, then recurse on the resulting Seurat object so we reuse
  # the hub write path. These formats have no corresponding saver because they
  # sit upstream of the scverse container ecosystem.
  if (stype %in% c('gef', 'cellbin.gef', 'cosmx')) {
    loader <- GetLoader(ext = stype)
    if (is.null(loader)) {
      stop("No loader registered for vendor format: ", stype, call. = FALSE)
    }
    if (missing(x = assay)) {
      assay <- if (stype == 'cosmx') 'Nanostring' else 'Spatial'
    }
    obj <- loader(source = source, assay = assay, verbose = verbose, ...)
    return(scConvert(
      source = obj, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, standardize = standardize, ...
    ))
  }

  # Loom without a file-to-file path (e.g. to rds or zarr): open as a loom
  # connection and convert through scConvert.loom
  if (stype == 'loom' && is.null(GetDirectPath(stype = 'loom', dtype = FileType(file = dest)))) {
    lfile <- scConnect(filename = source, type = 'loom', mode = 'r')
    on.exit(expr = lfile$close_all(), add = TRUE)
    if (missing(x = assay)) {
      assay <- tryCatch(
        expr = DefaultAssay(object = lfile),
        error = function(...) 'RNA'
      )
    }
    return(scConvert(
      source = lfile, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, standardize = standardize, ...
    ))
  }

  # HDF5-based formats: check for direct path first, otherwise use hub
  dtype <- FileType(file = dest)

  # h5ad inputs written by anndata >= 0.11 may use encodings the HDF5-level
  # converters do not understand (nullable-string-array obs/var index and
  # columns, `null` entries). Materialise a converter-compatible temporary
  # copy in that case; see .h5ad_normalize_for_convert().
  if (stype == 'h5ad') {
    normalized <- .h5ad_normalize_for_convert(source, verbose = verbose)
    if (isTRUE(normalized$temporary)) {
      source <- normalized$path
      on.exit(expr = unlink(normalized$path), add = TRUE)
    }
  }
  if (missing(x = assay)) {
    # Try to read default assay from the source file
    hfile_tmp <- tryCatch(
      expr = scConnect(filename = source, force = TRUE),
      error = function(...) NULL
    )
    if (!is.null(x = hfile_tmp)) {
      assay <- tryCatch(
        expr = DefaultAssay(object = hfile_tmp),
        error = function(...) {
          message("'assay' not set, setting to 'RNA'")
          "RNA"
        }
      )
      tryCatch(hfile_tmp$close_all(), error = function(e) NULL)
    } else {
      assay <- 'RNA'
    }
  }

  # Check for a registered direct HDF5 path (memory-efficient)
  direct_fn <- GetDirectPath(stype = stype, dtype = dtype)
  if (!is.null(x = direct_fn)) {
    if (tolower(x = dest) == dtype) {
      dest <- paste(file_path_sans_ext(x = source), dtype, sep = '.')
    }
    hfile <- scConnect(filename = source, force = TRUE)
    on.exit(expr = tryCatch(hfile$close_all(), error = function(e) NULL), add = TRUE)
    dfile <- direct_fn(
      source = hfile, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose,
      standardize = standardize, ...
    )
    if (is.character(x = dfile)) {
      return(invisible(x = dfile))
    }
    on.exit(expr = tryCatch(dfile$close_all(), error = function(e) NULL), add = TRUE)
    return(invisible(x = dfile$filename))
  }

  # No direct path: use hub conversion (loader handles its own file I/O)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source), dtype, sep = '.')
  }
  return(HubConvert(
    source_file = source,
    dest_file = dest,
    stype = stype,
    dtype = dtype,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname scConvert
#' @method scConvert H5File
#' @export
#'
scConvert.H5File <- function(source, dest = 'h5seurat', assay = 'RNA',
                             overwrite = FALSE, verbose = TRUE, ...) {
  .convert_open_file(source, dest = dest, stype = FileType(file = source$filename),
                     assay = assay, overwrite = overwrite, verbose = verbose, ...)
}

#' @rdname scConvert
#' @method scConvert h5Seurat
#' @export
#'
scConvert.h5Seurat <- function(source, dest = 'h5ad', assay = DefaultAssay(object = source),
                               overwrite = FALSE, verbose = TRUE, standardize = FALSE, ...) {
  .convert_open_file(source, dest = dest, stype = 'h5seurat', assay = assay,
                     overwrite = overwrite, verbose = verbose, standardize = standardize, ...)
}

# Convert an open file connection: use a registered direct (file-to-file)
# path when there is one, otherwise load into Seurat and save. A bare
# format name as `dest` writes next to the source with that extension.
.convert_open_file <- function(source, dest, stype, assay, overwrite, verbose, ...) {
  dtype <- FileType(file = dest)
  if (identical(tolower(dest), dtype)) {
    dest <- paste0(file_path_sans_ext(x = source$filename), '.', dtype)
  }
  direct <- GetDirectPath(stype = stype, dtype = dtype)
  if (!is.null(direct)) {
    if (verbose && stype != 'h5seurat') {
      message("Using direct ", toupper(stype), " -> ", toupper(dtype), " conversion")
    }
    return(direct(source = source, dest = dest, assay = assay, overwrite = overwrite,
                  verbose = verbose, ...))
  }
  extra <- list(...)
  extra$standardize <- NULL  # only the h5ad writers understand it
  do.call(HubConvert, c(list(source_file = source$filename, dest_file = dest, stype = stype,
                             dtype = dtype, assay = assay, overwrite = overwrite,
                             verbose = verbose), extra))
}

#' @rdname scConvert
#' @method scConvert Seurat
#' @export
scConvert.Seurat <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (missing(x = dest)) {
    stop("'dest' must be provided for Seurat object conversion", call. = FALSE)
  }
  type <- FileType(file = dest)

  # Handle in-memory SCE destination (no file created)
  if (type == 'sce') {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' is required to convert to SCE", call. = FALSE)
    }
    if (verbose) message("Converting Seurat to SingleCellExperiment...")
    return(Seurat::as.SingleCellExperiment(x = source, assay = assay))
  }

  if (tolower(x = dest) == type) {
    dest <- paste(Project(object = source), type, sep = '.')
  }

  # Use registered saver if available
  saver <- GetSaver(ext = type)
  if (!is.null(x = saver)) {
    if (type == 'h5ad') {
      saver(object = source, filename = dest, overwrite = overwrite,
            verbose = verbose, standardize = standardize, ...)
    } else {
      saver(object = source, filename = dest, overwrite = overwrite, verbose = verbose, ...)
    }
    if (verbose) message("Successfully created: ", dest)
    return(invisible(x = dest))
  }

  stop("Cannot convert Seurat objects to '", type, "' format", call. = FALSE)
}

#' @rdname scConvert
#' @method scConvert loom
#' @export
#'
scConvert.loom <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (missing(x = dest)) {
    stop("'dest' must be provided for loom conversion", call. = FALSE)
  }
  if (verbose) message("Loading loom file as Seurat object...")
  seurat_obj <- as.Seurat(x = source, verbose = verbose)
  return(scConvert(
    source = seurat_obj, dest = dest, assay = assay,
    overwrite = overwrite, verbose = verbose, standardize = standardize, ...
  ))
}

#' @rdname scConvert
#' @method scConvert SingleCellExperiment
#' @export
#'
scConvert.SingleCellExperiment <- function(
  source,
  dest,
  assay = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required", call. = FALSE)
  }
  if (missing(x = dest)) {
    stop("'dest' must be provided for SingleCellExperiment conversion", call. = FALSE)
  }
  if (verbose) message("Converting SingleCellExperiment to Seurat...")
  # Detect available SCE assays and set conversion params accordingly.
  # as.Seurat defaults to data="logcounts" which may not exist.
  sce_assays <- SummarizedExperiment::assayNames(source)
  sce_args <- list(x = source)
  if ("counts" %in% sce_assays) sce_args$counts <- "counts"
  if ("logcounts" %in% sce_assays) {
    sce_args$data <- "logcounts"
  } else {
    # Use alist to preserve NULL (list$key <- NULL removes the key in R)
    sce_args <- c(sce_args, alist(data = NULL))
  }
  extra_args <- list(...)
  sce_args <- c(sce_args, extra_args[!names(extra_args) %in% names(sce_args)])
  seurat_obj <- do.call(Seurat::as.Seurat, sce_args)
  effective_assay <- assay %||% DefaultAssay(object = seurat_obj)
  return(scConvert(
    source = seurat_obj, dest = dest, assay = effective_assay,
    overwrite = overwrite, verbose = verbose, standardize = standardize
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hub Conversion Helper
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Internal hub converter: load any format -> Seurat -> save to any format.
# Used when no direct HDF5-level path exists between (stype, dtype).
#
# @param source_file Path to source file
# @param dest_file Path to destination file
# @param stype Source format extension (lowercase)
# @param dtype Destination format extension (lowercase)
# @param assay Assay name
# @param overwrite Whether to overwrite existing dest
# @param verbose Show progress messages
# @param ... Additional arguments passed to saver
#
# @return Invisible destination filepath
#
# @keywords internal
#
HubConvert <- function(source_file, dest_file, stype, dtype,
                       assay = 'RNA', overwrite = FALSE, verbose = TRUE, ...) {
  loader <- GetLoader(ext = stype)
  saver <- GetSaver(ext = dtype)

  if (is.null(x = loader)) {
    stop("No loader registered for source format '", stype,
         "'. Supported: ", paste(ListFormats()$loaders, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(x = saver)) {
    stop("No saver registered for destination format '", dtype,
         "'. Supported: ", paste(ListFormats()$savers, collapse = ", "),
         call. = FALSE)
  }

  if (verbose) {
    message("Converting ", toupper(x = stype), " -> Seurat -> ", toupper(x = dtype))
  }

  seurat_obj <- loader(file = source_file, assay = assay, verbose = verbose)
  saver(object = seurat_obj, filename = dest_file,
        overwrite = overwrite, verbose = verbose, ...)

  return(invisible(x = dest_file))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Shared HDF5 Write Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Convert 'shape' attribute to 'dims' and add assay.used on an HDF5 graph group
#
# @param graph_group An H5Group representing the graph
# @param assay Character assay name to tag
#
# @keywords internal
#

# Write a sparse matrix as CSR group in HDF5
#
# @param h5parent H5Group parent
# @param group_name Name for the new group
# @param mat Sparse matrix (dgCMatrix or coercible)
# @param gzip Gzip compression level
#
# @keywords internal
#
WriteCSRGroup <- function(h5parent, group_name, mat, gzip = GetCompressionLevel()) {
  mat <- ConvertBPCellsMatrix(mat)
  mat <- as(mat, "dgCMatrix")
  csr <- DeconstructSparseCSR(mat, coerce = FALSE)  # already dgCMatrix, skip re-coercion
  grp <- h5parent$create_group(group_name)
  # Use fixed chunk size (64K elements) for streaming compression instead of

  # single-chunk which forces entire array into memory at once
  chunk_size <- 65536L
  grp$create_dataset("data", robj = csr$data,
                     chunk_dims = max(1L, min(length(csr$data), chunk_size)), gzip_level = gzip)
  grp$create_dataset("indices", robj = csr$indices,
                     chunk_dims = max(1L, min(length(csr$indices), chunk_size)), gzip_level = gzip)
  grp$create_dataset("indptr", robj = csr$indptr,
                     chunk_dims = max(1L, min(length(csr$indptr), chunk_size)), gzip_level = gzip)
  grp$create_attr(attr_name = 'encoding-type', robj = 'csr_matrix',
                  dtype = CachedGuessDType('csr_matrix'), space = ScalarSpace())
  grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                  dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
  grp$create_attr(attr_name = 'shape', robj = csr$shape,
                  dtype = GuessDType(csr$shape))
  invisible(grp)
}

# Write a sparse matrix as CSC group in HDF5 (for obsp graphs)
#
# Writes dgCMatrix directly as CSC without shape reversal.
# Use this for obsp matrices to avoid transposing non-symmetric graphs.
#
# @param h5parent H5Group parent
# @param group_name Name for the new group
# @param mat Sparse matrix (dgCMatrix or coercible)
# @param gzip Gzip compression level
#
# @keywords internal
#
WriteCscGroup <- function(h5parent, group_name, mat, gzip = GetCompressionLevel()) {
  mat <- ConvertBPCellsMatrix(mat)
  mat <- as(mat, "dgCMatrix")
  grp <- h5parent$create_group(group_name)
  chunk_size <- 65536L
  grp$create_dataset("data", robj = mat@x,
                     chunk_dims = max(1L, min(length(mat@x), chunk_size)), gzip_level = gzip)
  grp$create_dataset("indices", robj = mat@i,
                     chunk_dims = max(1L, min(length(mat@i), chunk_size)), gzip_level = gzip)
  grp$create_dataset("indptr", robj = mat@p,
                     chunk_dims = max(1L, min(length(mat@p), chunk_size)), gzip_level = gzip)
  grp$create_attr(attr_name = 'encoding-type', robj = 'csc_matrix',
                  dtype = CachedGuessDType('csc_matrix'), space = ScalarSpace())
  grp$create_attr(attr_name = 'encoding-version', robj = '0.1.0',
                  dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
  grp$create_attr(attr_name = 'shape', robj = as.integer(dim(mat)),
                  dtype = GuessDType(as.integer(dim(mat))))
  invisible(grp)
}

# Write a data frame as an AnnData-encoded HDF5 group (obs or var)
#
# @param h5parent H5Group parent
# @param group_name Name for the new group
# @param df Data frame to write
# Recursive writer: R object -> anndata-compatible HDF5 dataset/group in uns.
# Handles scalars, vectors, matrices, sparse matrices, data frames, factors,
# and named lists (recursed as HDF5 groups). Called from both the C writer
# post-processing path (zzz.R) and the R writer (DirectSeuratToH5AD).
#
# @param parent H5Group to write into
# @param name Dataset/group name
# @param value R object to write
# @param gzip Gzip compression level
#
# @keywords internal
#
WriteUnsItem <- function(parent, name, value, gzip = GetCompressionLevel()) {
  if (is.null(value)) return(invisible(NULL))

  if (is.data.frame(value)) {
    WriteDFGroup(parent, name, value,
                 rownames(value) %||% as.character(seq_len(nrow(value))),
                 gzip = gzip)
  } else if (is.list(value) && !is.null(names(value))) {
    grp <- parent$create_group(name)
    for (sub_name in names(value)) {
      tryCatch(
        WriteUnsItem(grp, sub_name, value[[sub_name]], gzip),
        error = function(e) NULL
      )
    }
  } else if (inherits(value, c("dgCMatrix", "dgRMatrix", "dgTMatrix"))) {
    WriteCSRGroup(parent, name, value, gzip = gzip)
  } else if (is.matrix(value)) {
    parent$create_dataset(name, robj = value,
                          chunk_dims = dim(value), gzip_level = gzip)
  } else if (is.character(value)) {
    parent$create_dataset(name, robj = value, dtype = CachedUtf8Type(),
                          chunk_dims = max(1L, length(value)), gzip_level = gzip)
    if (length(value) > 1) AddAnndataEncoding(parent[[name]], "string-array")
  } else if (is.logical(value)) {
    parent$create_dataset(name, robj = BoolToInt(value),
                          chunk_dims = max(1L, length(value)), gzip_level = gzip)
    if (length(value) > 1) AddAnndataEncoding(parent[[name]], "array")
  } else if (is.numeric(value) || is.integer(value)) {
    if (length(value) == 1) {
      parent$create_dataset(name, robj = value)
    } else {
      parent$create_dataset(name, robj = value,
                            chunk_dims = max(1L, length(value)), gzip_level = gzip)
      AddAnndataEncoding(parent[[name]], "array")
    }
  } else if (is.factor(value)) {
    enc <- EncodeCategorical(value)
    cat_grp <- parent$create_group(name)
    cat_grp$create_dataset("codes", robj = enc$codes,
                           chunk_dims = max(1L, length(enc$codes)), gzip_level = gzip)
    cat_grp$create_dataset("categories", robj = enc$categories,
                           dtype = CachedUtf8Type(),
                           chunk_dims = max(1L, length(enc$categories)),
                           gzip_level = gzip)
    cat_grp$create_attr("encoding-type", robj = "categorical",
                        dtype = CachedGuessDType("categorical"), space = ScalarSpace())
    cat_grp$create_attr("encoding-version", robj = "0.2.0",
                        dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
  }
  invisible(NULL)
}

# @param index_values Character vector of row names / index
# @param gzip Gzip compression level
#
# @keywords internal
#
WriteDFGroup <- function(h5parent, group_name, df, index_values,
                         gzip = GetCompressionLevel()) {
  grp <- h5parent$create_group(group_name)
  grp$create_dataset("_index", robj = index_values,
                     dtype = CachedUtf8Type(),
                     chunk_dims = length(index_values), gzip_level = gzip)
  grp$create_attr(attr_name = '_index', robj = '_index',
                  dtype = CachedGuessDType('_index'), space = ScalarSpace())
  col_order <- colnames(df)
  for (col_name in col_order) {
    col_data <- df[[col_name]]
    if (is.factor(col_data)) {
      enc <- EncodeCategorical(col_data)
      cat_grp <- grp$create_group(col_name)
      cat_grp$create_dataset("codes", robj = enc$codes,
                             chunk_dims = length(enc$codes), gzip_level = gzip)
      cat_grp$create_dataset("categories", robj = enc$categories,
                             dtype = CachedUtf8Type(),
                             chunk_dims = length(enc$categories), gzip_level = gzip)
      cat_grp$create_attr(attr_name = 'encoding-type', robj = 'categorical',
                          dtype = CachedGuessDType('categorical'), space = ScalarSpace())
      cat_grp$create_attr(attr_name = 'encoding-version', robj = '0.2.0',
                          dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
      # numpy-bool enum, so anndata/pandas see a real bool flag
      .h5w_attr(cat_grp, 'ordered', is.ordered(col_data))
    } else if (is.character(col_data) && anyNA(col_data)) {
      # string-array cannot hold missing values (NA would become the
      # literal "NA"); a categorical with -1 codes is readable by every
      # anndata version and decodes back to NA
      .h5ad_write_categorical_encoded(grp, col_name, factor(col_data))
    } else if (is.character(col_data)) {
      grp$create_dataset(col_name, robj = col_data,
                         dtype = CachedUtf8Type(),
                         chunk_dims = length(col_data), gzip_level = gzip)
      AddAnndataEncoding(grp[[col_name]], encoding_type = 'string-array')
    } else if (is.logical(col_data) && anyNA(col_data)) {
      .h5ad_write_nullable_encoded(grp, col_name, col_data, 'nullable-boolean')
    } else if (is.logical(col_data)) {
      # numpy-bool enum (was int with NA -> 2, read back as integers)
      d <- .h5w_dataset(grp, col_name, col_data)
      .h5w_encoding(d, 'array', '0.2.0')
    } else if (is.integer(col_data) && anyNA(col_data)) {
      # NA_integer_ would otherwise be written as INT_MIN
      .h5ad_write_nullable_encoded(grp, col_name, col_data, 'nullable-integer')
    } else {
      grp$create_dataset(col_name, robj = col_data,
                         chunk_dims = length(col_data), gzip_level = gzip)
      AddAnndataEncoding(grp[[col_name]], encoding_type = 'array')
    }
  }
  # anndata iterates `column-order`; hdf5r needs an explicit (possibly
  # empty) 1-D space when the data frame has no columns
  grp$create_attr(attr_name = 'column-order', robj = as.character(col_order),
                  dtype = CachedUtf8Type(),
                  space = H5S$new(type = 'simple', dims = length(col_order),
                                  maxdims = length(col_order)))
  grp$create_attr(attr_name = 'encoding-type', robj = 'dataframe',
                  dtype = CachedGuessDType('dataframe'), space = ScalarSpace())
  grp$create_attr(attr_name = 'encoding-version', robj = '0.2.0',
                  dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
  invisible(grp)
}

# Write one categorical column into an h5ad dataframe group.
#
# Creates /<parent_grp>/<name>/ with:
#   codes       integer array, 0-based with -1 for NA, dtype sized to n_levels
#   categories  UTF-8 string array of level names
#   attrs: encoding-type="categorical", encoding-version="0.2.0", ordered=FALSE
#
# Caller is responsible for name mapping (e.g. standardize) and for
# slicing `values` with any subset index before the call.
#
# @param parent_grp hdf5r H5Group to create the child group under (e.g. /obs)
# @param name       character(1), child group name
# @param levels     character vector, factor level labels
# @param values     integer vector, 1-based factor codes (NA becomes -1)
#
# @return Invisibly, the created categorical group
# @keywords internal
#

# Seurat image-slot scale-factor names vs. AnnData/Scanpy attribute names.
# The only source of truth for this 4-entry table. Previously inlined in
# three separate sites (H5ADToH5Seurat, H5SeuratToH5AD, SeuratSpatialToH5AD);
# the slot-to-attribute and attribute-to-slot inverses are both derived from
# this vector so they cannot drift.
SCALEFACTOR_NAME_SEURAT_TO_H5AD <- c(
  spot     = "spot_diameter_fullres",
  fiducial = "fiducial_diameter_fullres",
  hires    = "tissue_hires_scalef",
  lowres   = "tissue_lowres_scalef"
)
SCALEFACTOR_NAME_H5AD_TO_SEURAT <- setNames(
  names(SCALEFACTOR_NAME_SEURAT_TO_H5AD),
  unname(SCALEFACTOR_NAME_SEURAT_TO_H5AD)
)

# Resolve the h5Seurat assay-layout paths for data / counts / scale.data.
# Seurat v5 stores expression matrices under `layers/` (layers/data,
# layers/counts); v4 stores them at the assay root (data, counts). `scale.data`
# sits at the assay root in both versions. The same inline detection had to be
# patched in every site when Seurat v5 shipped; routing every caller through
# this helper prevents future drift.
#
# @param assay_group hdf5r H5Group for an assay (e.g. `assays/RNA`)
#
# @return List with components:
#   - is_v5: logical, TRUE if a layers/ H5Group is present
#   - data_path: "layers/data" / "data" / NULL (relative to assay_group)
#   - counts_path: "layers/counts" / "counts" / NULL
#   - scale_data_path: "scale.data" / NULL
#
# @keywords internal
DetectSeuratLayerPaths <- function(assay_group) {
  is_v5 <- assay_group$exists(name = 'layers') &&
    inherits(assay_group[['layers']], 'H5Group')
  if (is_v5) {
    layers <- assay_group[['layers']]
    data_path <- if (layers$exists(name = 'data')) 'layers/data' else NULL
    counts_path <- if (layers$exists(name = 'counts')) 'layers/counts' else NULL
  } else {
    data_path <- if (assay_group$exists(name = 'data')) 'data' else NULL
    counts_path <- if (assay_group$exists(name = 'counts')) 'counts' else NULL
  }
  scale_data_path <- if (assay_group$exists(name = 'scale.data')) 'scale.data' else NULL
  list(
    is_v5 = is_v5,
    data_path = data_path,
    counts_path = counts_path,
    scale_data_path = scale_data_path
  )
}

# Resolve rownames in an h5ad dataframe group. Older h5ad files encode the
# index as an `_index` attribute; newer ones write `_index` or `index` as a
# dataset. Returns either the attribute value or the dataset name.
#
# @keywords internal

# Decode categorical feature_name group (categories + codes) into a character
# vector. Memoises on the caller-provided environment so the same group is
# not re-read for each assay.
#
# @keywords internal

# Read a compound HDF5 dataset (old-style obs/var stored as a flat table) into
# a data.frame. Decodes `__categories` sibling-group encodings to characters.
#
# @keywords internal

# Decompress an LZF-compressed sparse matrix into a temp h5 file so downstream
# H5Ocopy can handle it without the LZF filter.
#
# @keywords internal

# Convert old-style h5ad `__categories` sibling-group encoding into h5Seurat
# factor groups (levels + 1-based values). Booleans are specialised to a
# 1-D dataset since two-level factors produce noise.
#
# @keywords internal

# Write placeholder cell.names into `dfile` when the h5ad has no obs index.
# Uses the X matrix shape to determine cell count.
#
# @keywords internal

# Sanitize column names in an HDF5 dataframe group and flatten nullable
# (mask + values) sub-groups into plain datasets.
#
# @keywords internal

# Convert h5ad-style categorical columns (categories + codes, 0-based) in a
# copied meta.data group to h5Seurat factor columns (levels + values, 1-based
# with NA for -1 codes).
#
# @keywords internal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Implementations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AnnData Preprocessing Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Sanitize a single column name to valid R/Seurat format
# Replaces /, spaces, commas, semicolons, colons, backslashes with underscore
# Removes consecutive and trailing underscores
# Ensures valid R name (can't start with number)
#
# @param x Character string to sanitize
#
# @return Sanitized character string
#
# @keywords internal
#
SanitizeColumnName <- function(x) {
  # Replace problematic characters with underscore
  x <- gsub(pattern = '[/\\s,;:\\\\]', replacement = '_', x = x, perl = TRUE)
  # Remove consecutive underscores
  x <- gsub(pattern = '_+', replacement = '_', x = x)
  # Remove leading/trailing underscores
  x <- gsub(pattern = '^_|_$', replacement = '', x = x)
  # Ensure valid R name (can't start with number)
  if (grepl('^[0-9]', x)) {
    x <- paste0('X', x)
  }
  return(x)
}

# Sanitize all column names, handling duplicates
# Returns a named vector mapping original -> sanitized names
#
# @param names Character vector of column names to sanitize
#
# @return Named character vector where names are original and values are sanitized
#
# @keywords internal
#
SanitizeColumnNames <- function(names) {
  sanitized <- vapply(names, SanitizeColumnName, character(1), USE.NAMES = FALSE)
  # Handle duplicates by appending __dupN
  if (anyDuplicated(sanitized)) {
    counts <- table(sanitized)
    dups <- names(counts[counts > 1])
    for (dup in dups) {
      idx <- which(sanitized == dup)
      sanitized[idx[-1]] <- paste0(dup, '__dup', seq_along(idx[-1]))
    }
  }
  names(sanitized) <- names  # Map original -> sanitized
  return(sanitized)
}

# Flatten nullable dtype (mask+values group) to simple vector
# h5ad stores nullable dtypes as groups with 'mask' (boolean) and 'values' arrays
# In h5ad nullable dtypes, mask=TRUE means the value is MISSING (NA)
#
# @param col_group H5Group object that may contain mask+values structure
#
# @return Flattened vector with NA values where mask is TRUE, or NULL if not a mask+values structure
#
# @keywords internal
#
FlattenNullable <- function(col_group) {
  if (!inherits(col_group, 'H5Group')) {
    return(NULL)  # Not a group, skip
  }

  # Check if this is a mask+values structure (not categorical which has categories/codes)
  if (col_group$exists('mask') && col_group$exists('values')) {
    # Categoricals (codes + categories) are handled elsewhere; every
    # nullable-* encoding (nullable-integer, nullable-boolean and, since
    # anndata 0.11, nullable-string-array) is a values + mask pair.
    if (col_group$exists('categories') || col_group$exists('codes')) {
      return(NULL)
    }
    values <- .h5ad_read_nullable(col_group)
    if (is.character(values)) {
      # h5Seurat string datasets cannot hold NA; a factor (NA-capable
      # levels/values group) keeps the missing values intact.
      values <- factor(values)
    }
    return(values)
  }

  return(NULL)  # Not a mask+values structure
}


