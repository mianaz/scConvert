#' @include zzz.R
#' @include FormatRegistry.R
#' @include scConnect.R
#' @include TestObject.R
#' @include scTranspose.R
#' @include PadMatrix.R
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
#' @param ... Arguments passed to specific conversion methods
#'
#' @return For file destinations, invisibly returns the destination filename.
#'   For \code{dest = "sce"}, returns a SingleCellExperiment object.
#'
#' @inheritSection H5ADToH5Seurat AnnData/H5AD to h5Seurat
#' @inheritSection H5SeuratToH5AD h5Seurat to AnnData/H5AD
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

  # Loom: open as loom R6 object and dispatch to Convert.loom
  if (stype == 'loom') {
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
      hfile_tmp$close_all()
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
    on.exit(expr = hfile$close_all(), add = TRUE)
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
scConvert.H5File <- function(
  source,
  dest = 'h5seurat',
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  stype <- FileType(file = source$filename)
  dtype <- FileType(file = dest)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source$filename), dtype, sep = '.')
  }

  # 1. Try registered direct HDF5-level path (memory-efficient)
  direct_fn <- GetDirectPath(stype = stype, dtype = dtype)
  if (!is.null(x = direct_fn)) {
    if (verbose) {
      message("Using direct ", toupper(x = stype), " -> ", toupper(x = dtype), " conversion")
    }
    dfile <- direct_fn(
      source = source, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose, ...
    )
    return(dfile)
  }

  # 2. Hub conversion: load source -> Seurat -> save dest
  return(HubConvert(
    source_file = source$filename,
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
#' @method scConvert h5Seurat
#' @export
#'
scConvert.h5Seurat <- function(
  source,
  dest = 'h5ad',
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
) {
  dtype <- FileType(file = dest)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source$filename), dtype, sep = '.')
  }

  # Try direct path first (h5seurat -> h5ad is registered)
  direct_fn <- GetDirectPath(stype = 'h5seurat', dtype = dtype)
  if (!is.null(x = direct_fn)) {
    return(direct_fn(
      source = source, dest = dest, assay = assay,
      overwrite = overwrite, verbose = verbose,
      standardize = standardize, ...
    ))
  }

  # Hub: load h5seurat -> Seurat -> save dest
  return(HubConvert(
    source_file = source$filename,
    dest_file = dest,
    stype = 'h5seurat',
    dtype = dtype,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
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
    saver(object = source, filename = dest, overwrite = overwrite, verbose = verbose, ...)
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
FixGraphAttrs <- function(graph_group, assay) {
  if (isTRUE(x = AttrExists(x = graph_group, name = 'shape'))) {
    if (!isTRUE(x = AttrExists(x = graph_group, name = 'dims'))) {
      shape_val <- h5attr(x = graph_group, which = 'shape')
      graph_group$create_attr(
        attr_name = 'dims',
        robj = shape_val,
        dtype = GuessDType(x = shape_val)
      )
    }
    graph_group$attr_delete(attr_name = 'shape')
  }
  if (!isTRUE(x = AttrExists(x = graph_group, name = 'assay.used'))) {
    graph_group$create_attr(
      attr_name = 'assay.used',
      robj = assay,
      dtype = GuessDType(x = assay)
    )
  }
}

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
      cat_grp$create_attr(attr_name = 'ordered', robj = is.ordered(col_data),
                          dtype = GuessDType(TRUE), space = ScalarSpace())
    } else if (is.character(col_data)) {
      grp$create_dataset(col_name, robj = col_data,
                         dtype = CachedUtf8Type(),
                         chunk_dims = length(col_data), gzip_level = gzip)
      AddAnndataEncoding(grp[[col_name]], encoding_type = 'string-array')
    } else if (is.logical(col_data)) {
      grp$create_dataset(col_name, robj = BoolToInt(col_data),
                         chunk_dims = length(col_data), gzip_level = gzip)
    } else {
      grp$create_dataset(col_name, robj = col_data,
                         chunk_dims = length(col_data), gzip_level = gzip)
    }
  }
  grp$create_attr(attr_name = 'column-order', robj = col_order,
                  dtype = CachedUtf8Type())
  grp$create_attr(attr_name = 'encoding-type', robj = 'dataframe',
                  dtype = CachedGuessDType('dataframe'), space = ScalarSpace())
  grp$create_attr(attr_name = 'encoding-version', robj = '0.2.0',
                  dtype = CachedGuessDType('0.2.0'), space = ScalarSpace())
  invisible(grp)
}

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
    # Skip if this looks like a categorical (has categories attribute or encoding-type)
    if (col_group$exists('categories') ||
        isTRUE(x = AttrExists(x = col_group, name = 'encoding-type'))) {
      return(NULL)
    }

    values <- col_group[['values']][]
    mask <- col_group[['mask']][]

    # In h5ad nullable dtypes, mask=TRUE means the value is MISSING (NA)
    values[mask] <- NA
    return(values)
  }

  return(NULL)  # Not a mask+values structure
}

#' Convert AnnData/H5AD files to h5Seurat files
#'
#' @inheritParams scConvert
#'
#' @return Returns a handle to \code{dest} as an \code{\link{h5Seurat}} object
#'
#' @importFrom Seurat Project<- DefaultAssay<-
#'
#' @section AnnData/H5AD to h5Seurat:
#' The AnnData/H5AD to h5Seurat conversion will try to automatically fill in
#' datasets based on data presence. It works in the following manner:
#' \subsection{Expression data}{
#'  The expression matrices \code{counts}, \code{data}, and \code{scale.data}
#'  are filled by \code{/X} and \code{/raw/X} in the following manner:
#'  \itemize{
#'   \item \code{counts} will be filled with \code{/raw/X} if present;
#'   otherwise, it will be filled with \code{/X}
#'   \item \code{data} will be filled with \code{/raw/X} if \code{/raw/X} is
#'   present and \code{/X} is dense; otherwise, it will be filled with \code{/X}
#'   \item \code{scale.data} will be filled with \code{/X} if it dense;
#'   otherwise, it will be empty
#'  }
#'  Feature names are taken from the feature-level metadata
#' }
#' \subsection{Feature-level metadata}{
#' Feature-level metadata is added to the \code{meta.features} datasets in each
#' assay. Feature names are taken from the dataset specified by the
#' \dQuote{_index} attribute, the \dQuote{_index} dataset, or the \dQuote{index}
#' dataset, in that order. Metadata is populated with \code{/raw/var} if
#' present, otherwise with \code{/var}; if both \code{/raw/var} and \code{/var}
#' are present, then \code{meta.features} will be populated with \code{/raw/var}
#' first, then \code{/var} will be added to it. For columns present in both
#' \code{/raw/var} and \code{/var}, the values in \code{/var} will be used
#' instead. \strong{Note}: it is possible for \code{/var} to have fewer features
#' than \code{/raw/var}; if this is the case, then only the features present in
#' \code{/var} will be overwritten, with the metadata for features \emph{not}
#' present in \code{/var} remaining as they were in \code{/raw/var} or empty
#' }
#' \subsection{Cell-level metadata}{
#' Cell-level metadata is added to \code{meta.data}; the row names of the
#' metadata (as determined by the value of the \dQuote{_index} attribute, the
#' \dQuote{_index} dataset, or the \dQuote{index} dataset, in that order) are
#' added to the \dQuote{cell.names} dataset instead. If the
#' \dQuote{__categories} dataset is present, each dataset within
#' \dQuote{__categories} will be stored as a factor group. Cell-level metadata
#' will be added as an HDF5 group unless factors are \strong{not} present and
#' the \code{scConvert.dtypes.dataframe_as_group} option is \code{FALSE}
#' }
#' \subsection{Dimensional reduction information:}{
#'  Cell embeddings are taken from \code{/obsm}; dimensional reductions are
#'  named based on their names from \code{obsm} by removing the preceding
#'  \dQuote{X_}.For example, if a dimensional reduction is named \dQuote{X_pca}
#'  in \code{/obsm}, the resulting dimensional reduction information will be
#'  named \dQuote{pca}. The key will be set to one of the following:
#'  \itemize{
#'   \item \dQuote{PC_} if \dQuote{pca} is present in the dimensional reduction
#'   name (\code{grepl("pca", reduction.name, ignore.case = TRUE)})
#'   \item \dQuote{tSNE_} if \dQuote{tsne} is present in the dimensional
#'   reduction name (\code{grepl("tsne", reduction.name, ignore.case = TRUE)})
#'   \item \code{reduction.name_} for all other reductions
#'  }
#'  Remember that the preceding \dQuote{X_} will be removed from the reduction
#'  name before converting to a key. Feature loadings are taken from
#'  \code{/varm} and placed in the associated dimensional reduction. The
#'  dimensional reduction is determine from the loadings name in \code{/varm}:
#'  \itemize{
#'   \item \dQuote{PCs} will be added to a dimensional reduction named
#'   \dQuote{pca}
#'   \item All other loadings in \code{/varm} will be added to a dimensional
#'   reduction named \code{tolower(loading)} (eg. a loading named \dQuote{ICA}
#'   will be added to a dimensional reduction named \dQuote{ica})
#'  }
#'  If a dimensional reduction cannot be found according to the rules above, the
#'  loading will not be taken from the AnnData/H5AD file. Miscellaneous
#'  information will be taken from \code{/uns/reduction} where \code{reduction}
#'  is the name of the reduction in \code{/obsm} without the preceding
#'  \dQuote{X_}; if no dimensional reduction information present, then
#'  miscellaneous information will not be taken from the AnnData/H5AD file.
#'  Standard deviations are taken from a dataset \code{/uns/reduction/variance};
#'  the variances will be converted to standard deviations and added to the
#'  \code{stdev} dataset of a dimensional reduction
#' }
#' \subsection{Nearest-neighbor graph}{
#'  If a nearest neighbor graph is present in \code{/uns/neighbors/distances},
#'  it will be added as a graph dataset in the h5Seurat file and associated with
#'  \code{assay}; if a value is present in \code{/uns/neighbors/params/method},
#'  the name of the graph will be \code{assay_method}, otherwise, it will be
#'  \code{assay_anndata}
#' }
#' \subsection{Layers}{
#'  TODO: add this
#' }
#' \subsection{Miscellaneous information}{
#'  All groups and datasets from \code{/uns} will be copied to \code{misc} in
#'  the h5Seurat file except for the following:
#'  \itemize{
#'   \item Any group or dataset named the same as a dimensional reduction (eg.
#'   \code{/uns/pca})
#'   \item \code{/uns/neighbors}
#'  }
#' }
#'
#' @keywords internal
#'
H5ADToH5Seurat <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE
) {
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination h5Seurat file exists", call. = FALSE)
    }
  }
  dfile <- h5Seurat$new(filename = dest, mode = WriteMode(overwrite = FALSE))
  # Get rownames from an H5AD data frame
  #
  # @param dset Name of data frame
  #
  # @return Returns the name of the dataset that contains the rownames
  #
  GetRownames <- function(dset) {
    if (inherits(x = source[[dset]], what = 'H5Group')) {
      rownames <- if (isTRUE(x = AttrExists(x = source[[dset]], name = '_index'))) {
        h5attr(x = source[[dset]], which = '_index')
      } else if (source[[dset]]$exists(name = '_index')) {
        '_index'
      } else if (source[[dset]]$exists(name = 'index')) {
        'index'
      } else {
        stop("Cannot find rownames in ", dset, call. = FALSE)
      }
    } else {
      stop("Don't know how to handle datasets", call. = FALSE)
    }
    return(rownames)
  }
  # Read categorical feature names from H5AD
  #
  # @param dset Name of var dataframe group (e.g., 'raw/var' or 'var')
  #
  # @return Returns feature names (gene symbols) if categorical encoding exists, NULL otherwise
  #
  # Cache for ReadCategoricalFeatures results to avoid redundant HDF5 reads
  .categorical_cache <- new.env(parent = emptyenv())

  ReadCategoricalFeatures <- function(dset) {
    # Return cached result if available
    cache_key <- dset
    if (exists(cache_key, envir = .categorical_cache)) {
      return(get(cache_key, envir = .categorical_cache))
    }

    # Check if feature_name exists and is a categorical group
    feature_name_path <- paste0(dset, '/feature_name')
    if (!source$exists(name = feature_name_path)) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    feature_name_group <- source[[feature_name_path]]
    if (!inherits(x = feature_name_group, what = 'H5Group')) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    # Check for categorical encoding using simple attr_exists (skip space check)
    if (!feature_name_group$attr_exists(attr_name = 'encoding-type')) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    encoding_type <- h5attr(x = feature_name_group, which = 'encoding-type')
    if (encoding_type != 'categorical') {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    # Read categories and codes
    if (!feature_name_group$exists(name = 'categories') ||
        !feature_name_group$exists(name = 'codes')) {
      assign(cache_key, NULL, envir = .categorical_cache)
      return(NULL)
    }

    categories <- feature_name_group[['categories']]$read()
    codes <- feature_name_group[['codes']]$read()

    # Decode: categories[codes] (R uses 1-based indexing, codes are 0-based)
    feature_names <- as.character(categories[codes + 1])

    # Cache and return
    assign(cache_key, feature_names, envir = .categorical_cache)
    return(feature_names)
  }
  # Read compound HDF5 dataset using hdf5r
  #
  # Some older h5ad files store obs/var as compound HDF5 datasets (a single flat
  # table rather than a group with sub-datasets). hdf5r reads these natively.
  #
  # @param h5file Path to the h5ad file
  # @param dataset_path Path to the dataset within the file (e.g., 'obs', 'var')
  #
  # @return A data.frame with the dataset contents, or NULL on failure
  #
  ReadCompoundDataset <- function(h5file, dataset_path) {
    tryCatch({
      h5 <- H5File$new(h5file, mode = "r")
      on.exit(h5$close_all(), add = TRUE)
      if (!h5$exists(dataset_path)) return(NULL)
      ds <- h5[[dataset_path]]
      if (!inherits(ds, "H5D")) return(NULL)

      data <- ds$read()
      if (is.list(data)) {
        # Compound dataset returns a named list; convert to data.frame
        df <- as.data.frame(data, stringsAsFactors = FALSE)
        # Decode raw/bytes columns to character
        for (col in names(df)) {
          if (is.raw(df[[col]]) || is.list(df[[col]])) next
          if (is.character(df[[col]])) next
          # Leave numeric columns as-is
        }

        # Handle __categories if present (old-style categorical encoding)
        cats_path <- paste0(dataset_path, "/__categories")
        if (h5$exists(cats_path)) {
          cats_grp <- h5[[cats_path]]
          for (col in names(cats_grp)) {
            if (col %in% names(df)) {
              categories <- as.character(cats_grp[[col]]$read())
              codes <- as.integer(df[[col]])
              decoded <- rep(NA_character_, length(codes))
              valid <- !is.na(codes) & codes >= 0L & codes < length(categories)
              decoded[valid] <- categories[codes[valid] + 1L]
              df[[col]] <- decoded
            }
          }
        }
        return(df)
      }
      NULL
    }, error = function(e) {
      if (verbose) {
        message("  Warning: Failed to read compound dataset: ", conditionMessage(e))
      }
      NULL
    })
  }
  # Read LZF-compressed sparse matrix and write decompressed to temp file
  #
  # Some older h5ad files use LZF compression. hdf5r reads these natively
  # when the HDF5 library has LZF filter support (standard on modern installs).
  # The data is read and rewritten uncompressed so H5Ocopy can work downstream.
  #
  # @param src_file Path to source h5ad file
  # @param src_path Path to sparse matrix group in source (e.g., 'X')
  #
  # @return Path to temporary H5 file with decompressed data, or NULL if failed
  #
  CopySparseMatrixDecompressed <- function(src_file, src_path) {
    temp_file <- tempfile(fileext = ".h5")
    tryCatch({
      src_h5 <- H5File$new(src_file, mode = "r")
      on.exit(src_h5$close_all(), add = TRUE)
      dst_h5 <- H5File$new(temp_file, mode = "w")
      on.exit(dst_h5$close_all(), add = TRUE)

      src_grp <- src_h5[[src_path]]
      dst_grp <- dst_h5$create_group("matrix")

      # Read each component (hdf5r handles LZF decompression if HDF5 supports it)
      for (comp in c("data", "indices", "indptr")) {
        if (src_grp$exists(comp)) {
          vals <- src_grp[[comp]]$read()
          dst_grp$create_dataset(comp, robj = vals)
        }
      }

      # Copy attributes
      for (attr_name in h5attr_names(src_grp)) {
        attr_val <- h5attr(src_grp, attr_name)
        dst_grp$create_attr(
          attr_name = attr_name, robj = attr_val,
          dtype = GuessDType(x = if (length(attr_val) > 0) attr_val[1] else attr_val)
        )
      }

      return(temp_file)
    }, error = function(e) {
      if (verbose) {
        message("  Warning: LZF decompression failed: ", conditionMessage(e),
                "\n  Ensure HDF5 was built with LZF filter support")
      }
      unlink(temp_file)
      return(NULL)
    })
  }
  ColToFactor <- function(dfgroup) {
    if (dfgroup$exists(name = '__categories')) {
      for (i in names(x = dfgroup[['__categories']])) {
        tname <- basename(path = tempfile(tmpdir = ''))
        dfgroup$obj_copy_to(dst_loc = dfgroup, dst_name = tname, src_name = i)
        dfgroup$link_delete(name = i)
        # Because AnnData stores logicals as factors, but have too many levels
        # for factors
        bool.check <- dfgroup[['__categories']][[i]]$dims == 2
        if (isTRUE(x = bool.check)) {
          bool.check <- all(sort(x = dfgroup[['__categories']][[i]][]) == c('False', 'True'))
        }
        if (isTRUE(x = bool.check)) {
          dfgroup$create_dataset(
            name = i,
            robj = dfgroup[[tname]][] + 1L,
            dtype = dfgroup[[tname]]$get_type()
          )
        } else {
          dfgroup$create_group(name = i)
          dfgroup[[i]]$create_dataset(
            name = 'values',
            robj = dfgroup[[tname]][] + 1L,
            dtype = dfgroup[[tname]]$get_type()
          )
          if (IsDType(x = dfgroup[['__categories']][[i]], dtype = 'H5T_STRING')) {
            dfgroup$obj_copy_to(
              dst_loc = dfgroup,
              dst_name = paste0(i, '/levels'),
              src_name = paste0('__categories/', i)
            )
          } else {
            dfgroup[[i]]$create_dataset(
              name = 'levels',
              robj = as.character(x = dfgroup[[H5Path('__categories', i)]][]),
              dtype = StringType()
            )
          }
        }
        dfgroup$link_delete(name = tname)
      }
      dfgroup$link_delete(name = '__categories')
    }
    return(invisible(x = NULL))
  }
  # Helper to create fake cell names when metadata is unavailable
  CreateFakeCellNames <- function() {
    ncells <- if (inherits(x = assay.group[['data']], what = 'H5Group')) {
      assay.group[['data/indptr']]$dims - 1
    } else {
      assay.group[['data']]$dims[2]
    }
    dfile$create_group(name = 'meta.data')
    dfile$create_dataset(
      name = 'cell.names',
      robj = paste0('Cell', seq.default(from = 1, to = ncells)),
      dtype = CachedGuessDType(x = 'Cell1')
    )
  }
  # Sanitize column names and handle nullable dtypes in HDF5 dataframes
  SanitizeH5DFColumns <- function(dfgroup) {
    col_names <- names(dfgroup)
    # Skip internal names
    col_names <- col_names[!col_names %in% c('__categories', '_index')]
    col_name_map <- SanitizeColumnNames(col_names)
    # Only process columns that need renaming or are nullable groups
    needs_rename <- vapply(names(col_name_map), function(n) n != col_name_map[[n]], logical(1))
    for (old_name in names(col_name_map)) {
      new_name <- col_name_map[[old_name]]
      col_obj <- dfgroup[[old_name]]

      # Handle nullable dtype (mask+values) groups
      if (inherits(col_obj, 'H5Group')) {
        flattened <- FlattenNullable(col_obj)
        if (!is.null(flattened)) {
          # Delete old group and create flattened dataset
          dfgroup$link_delete(name = old_name)
          dfgroup$create_dataset(
            name = new_name,
            robj = flattened,
            dtype = GuessDType(x = if (length(flattened) > 0) flattened[1] else flattened)
          )
          next
        }
        # Skip groups that don't need renaming (categorical groups handled elsewhere)
        if (old_name == new_name) next
      }

      # Rename if needed
      if (old_name != new_name) {
        dfgroup$obj_copy_from(
          src_loc = dfgroup,
          src_name = old_name,
          dst_name = new_name
        )
        dfgroup$link_delete(name = old_name)
      }
    }
  }
  ds.map <- c(
    scale.data = if (inherits(x = source[['X']], what = 'H5D')) {
      'X'
    } else {
      NULL
    },
    # Always use X for data slot - X contains normalized data when raw/X has counts
    data = 'X',
    counts = if (source$exists(name = 'raw')) {
      'raw/X'
    } else {
      'X'
    }
  )
  # Add assay data
  assay.group <- dfile[['assays']]$create_group(name = assay)
  for (i in seq_along(along.with = ds.map)) {
    if (verbose) {
      message("Adding ", ds.map[[i]], " as ", names(x = ds.map)[i])
    }
    dst <- names(x = ds.map)[i]

    # Check if this is a sparse matrix with LZF compression
    has_lzf <- FALSE
    if (inherits(x = source[[ds.map[[i]]]], what = 'H5Group')) {
      # Check if this sparse matrix group has LZF-compressed components
      if (source[[ds.map[[i]]]]$exists(name = 'data')) {
        tryCatch({
          # Try to read first element - will fail if LZF compressed
          test_read <- source[[ds.map[[i]]]][['data']][1]
        }, error = function(e) {
          if (grepl("can't synchronously read data|Read failed", conditionMessage(e))) {
            has_lzf <<- TRUE
          }
        })
      }
    }

    # Handle LZF compression with Python
    if (has_lzf) {
      if (verbose) {
        message("  Detected LZF compression; using Python to decompress")
      }

      # Get source file path
      src_file <- if (inherits(x = source, what = 'H5File')) {
        source$filename
      } else {
        as.character(source)
      }

      # Use Python to decompress to temporary file
      temp_file <- CopySparseMatrixDecompressed(
        src_file = src_file,
        src_path = ds.map[[i]]
      )

      if (!is.null(temp_file)) {
        # Copy from temporary file
        temp_h5 <- H5File$new(filename = temp_file, mode = 'r')
        assay.group$obj_copy_from(
          src_loc = temp_h5,
          src_name = 'matrix',
          dst_name = dst
        )
        temp_h5$close_all()
        unlink(temp_file)
      } else {
        if (verbose) {
          message("  Warning: LZF decompression failed, attempting normal copy")
        }
        # Fallback to normal copy (will likely fail on load, but conversion completes)
        assay.group$obj_copy_from(
          src_loc = source,
          src_name = ds.map[[i]],
          dst_name = dst
        )
      }
    } else {
      # Normal copy for non-LZF data
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = ds.map[[i]],
        dst_name = dst
      )
    }

    dst_obj <- assay.group[[dst]]
    if (dst_obj$attr_exists(attr_name = 'shape')) {
      dims <- rev(x = h5attr(x = dst_obj, which = 'shape'))
      dst_obj$create_attr(
        attr_name = 'dims',
        robj = dims,
        dtype = GuessDType(x = dims)
      )
      dst_obj$attr_delete(attr_name = 'shape')
    }
  }
  features.source <- ifelse(
    test = source$exists(name = 'raw') && source$exists(name = 'raw/var'),
    yes = 'raw/var',
    no = 'var'
  )
  # Track if we've already handled compound var dataset
  var_compound_handled <- FALSE

  if (inherits(x = source[[features.source]], what = 'H5Group')) {
    # Try to read categorical feature_name (gene symbols) first
    categorical_features <- ReadCategoricalFeatures(dset = features.source)

    if (!is.null(categorical_features)) {
      # Use gene symbols from categorical encoding
      if (verbose) {
        message("Using gene symbols from categorical feature_name encoding")
      }
      assay.group$create_dataset(
        name = 'features',
        robj = categorical_features,
        dtype = GuessDType(x = categorical_features[1])
      )
    } else {
      # Fallback to index dataset (typically Ensembl IDs)
      features.dset <- GetRownames(dset = features.source)
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = paste(features.source, features.dset, sep = '/'),
        dst_name = 'features'
      )
    }
  } else {
    # When var is an H5D dataset (e.g., compound datatype), try to read it with Python
    if (verbose) {
      message("Detected compound var dataset; attempting to read with Python")
    }

    # Get the source file path (need to handle both filename and H5File)
    h5file_path <- if (inherits(x = source, what = 'H5File')) {
      source$filename
    } else {
      as.character(source)
    }

    var_data <- ReadCompoundDataset(h5file_path, features.source)

    if (!is.null(var_data) && nrow(var_data) > 0) {
      # Successfully read compound dataset with Python
      if (verbose) {
        message("  Successfully read compound var dataset with ", nrow(var_data), " features")
      }

      # Extract feature names: prioritize feature_name (gene symbols) over index (Ensembl IDs)
      feature_names <- if ('feature_name' %in% names(var_data)) {
        as.character(var_data$feature_name)
      } else if ('index' %in% names(var_data)) {
        as.character(var_data$index)
      } else {
        as.character(var_data[[1]])
      }

      # Store the var data for later use in meta.features
      # Create meta.features group and add columns
      if (!assay.group$exists(name = 'meta.features')) {
        assay.group$create_group(name = 'meta.features')
      }
      # Sanitize column names for var metadata
      col_name_map <- SanitizeColumnNames(names(var_data))
      for (col_name in names(var_data)) {
        # Skip the column used for feature names (rownames)
        # If using feature_name for rownames, preserve index (Ensembl IDs) as metadata
        if ('feature_name' %in% names(var_data)) {
          if (col_name == 'feature_name') {
            next  # Skip feature_name as it's used for rownames
          }
        } else if (col_name == 'index' || col_name == names(var_data)[1]) {
          next  # Skip the index column as it's used for rownames
        }
        col_data <- var_data[[col_name]]
        sanitized_name <- col_name_map[[col_name]]
        assay.group[['meta.features']]$create_dataset(
          name = sanitized_name,
          robj = col_data,
          dtype = GuessDType(x = if (length(col_data) > 0) col_data[1] else col_data)
        )
      }
      # Add rownames to meta.features
      assay.group[['meta.features']]$create_dataset(
        name = '_index',
        robj = feature_names,
        dtype = GuessDType(x = feature_names[1])
      )
      assay.group[['meta.features']]$create_attr(
        attr_name = '_index',
        robj = '_index',
        dtype = CachedGuessDType(x = '_index')
      )
      assay.group[['meta.features']]$create_attr(
        attr_name = 's4class',
        robj = 'DFrame',
        dtype = StringType()
      )
      # Set a flag to indicate we've handled meta.features from compound var
      var_compound_handled <- TRUE
    } else {
      # Fallback: create generic feature names
      feature_names <- tryCatch(
        expr = {
          # Try to get rownames directly
          rn <- rownames(x = source[[features.source]])
          # Check if rownames actually returned valid data
          if (is.null(rn) || length(rn) == 0) {
            stop("rownames returned NULL or empty")
          }
          rn
        },
        error = function(e1) {
          # If that fails, try to get dimensions and create generic names
          tryCatch(
            expr = {
              # Get number of features from var dataset dimensions
              nfeatures <- source[[features.source]]$dims
              if (is.null(nfeatures) || length(nfeatures) == 0) {
                stop("dims returned NULL or empty")
              }
              if (verbose) {
                message(
                  "  Warning: Could not read var with Python; using generic feature names ",
                  "(Feature1, Feature2, ..., Feature", nfeatures, ")"
                )
              }
              paste0("Feature", seq_len(nfeatures))
            },
            error = function(e2) {
              # Last resort: get from X matrix dimensions
              if (inherits(x = source[['X']], what = 'H5Group')) {
                # Sparse matrix - check shape attribute
                if (source[['X']]$attr_exists(attr_name = 'shape')) {
                  shape <- h5attr(x = source[['X']], which = 'shape')
                  nfeatures <- shape[1]
                } else {
                  stop("Cannot determine number of features from h5ad file", call. = FALSE)
                }
              } else {
                # Dense matrix
                nfeatures <- source[['X']]$dims[1]
              }
              if (verbose) {
                message(
                  "  Warning: Cannot extract feature names from var; using generic names ",
                  "(Feature1, Feature2, ..., Feature", nfeatures, ")"
                )
              }
              paste0("Feature", seq_len(nfeatures))
            }
          )
        }
      )
    }

    assay.group$create_dataset(
      name = 'features',
      robj = feature_names,
      dtype = GuessDType(x = feature_names[1])
    )
  }
  scaled <- !is.null(x = ds.map['scale.data']) && !is.na(x = ds.map['scale.data'])
  if (scaled) {
    if (inherits(x = source[['var']], what = 'H5Group')) {
      # Try to read categorical feature_name for scaled features
      scaled_categorical_features <- ReadCategoricalFeatures(dset = 'var')

      if (!is.null(scaled_categorical_features)) {
        # Use gene symbols from categorical encoding
        if (verbose) {
          message("Using gene symbols from categorical feature_name encoding for scaled features")
        }
        assay.group$create_dataset(
          name = 'scaled.features',
          robj = scaled_categorical_features,
          dtype = GuessDType(x = scaled_categorical_features[1])
        )
      } else {
        # Fallback to index dataset
        scaled.dset <- GetRownames(dset = 'var')
        assay.group$obj_copy_from(
          src_loc = source,
          src_name = paste0('var/', scaled.dset),
          dst_name = 'scaled.features'
        )
      }
    } else {
      # When var is an H5D dataset, use the same feature names as before
      scaled_feature_names <- tryCatch(
        expr = {
          rn <- rownames(x = source[['var']])
          # Check if rownames actually returned valid data
          if (is.null(rn) || length(rn) == 0) {
            stop("rownames returned NULL or empty")
          }
          rn
        },
        error = function(e) {
          # Get number of scaled features from var or X
          nfeatures <- if (!is.null(source[['var']]$dims) && length(source[['var']]$dims) > 0) {
            source[['var']]$dims
          } else if (inherits(x = source[['X']], what = 'H5Group')) {
            if (source[['X']]$attr_exists(attr_name = 'shape')) {
              h5attr(x = source[['X']], which = 'shape')[1]
            } else {
              stop("Cannot determine number of scaled features", call. = FALSE)
            }
          } else {
            source[['X']]$dims[1]
          }
          paste0("Feature", seq_len(nfeatures))
        }
      )

      assay.group$create_dataset(
        name = 'scaled.features',
        robj = scaled_feature_names,
        dtype = GuessDType(x = scaled_feature_names[1])
      )
    }
  }
  assay.group$create_attr(
    attr_name = 'key',
    robj = paste0(tolower(x = assay), '_'),
    dtype = GuessDType(x = assay)
  )
  # Set default assay
  DefaultAssay(object = dfile) <- assay
  # Add feature-level metadata
  if (!getOption(x = "scConvert.dtypes.dataframe_as_group", default = FALSE)) {
    message("Adding feature-level metadata as a compound is not yet supported")
  }
  # TODO: Support compound metafeatures
  if (Exists(x = source, name = 'raw/var')) {
    if (inherits(x = source[['raw/var']], what = 'H5Group')) {
      if (verbose) {
        message("Adding meta.features from raw/var")
      }
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = 'raw/var',
        dst_name = 'meta.features'
      )
      if (scaled) {
        features.use <- SafeH5DRead(assay.group[['features']]) %in% SafeH5DRead(assay.group[['scaled.features']])
        features.use <- which(x = features.use)
        meta.scaled <- names(x = source[['var']])
        meta.scaled <- meta.scaled[!meta.scaled %in% c('__categories', scaled.dset)]
        for (mf in meta.scaled) {
          # Skip if not a dataset (e.g., skip groups like feature_types, genome)
          if (!inherits(x = source[['var']][[mf]], what = 'H5D')) {
            if (verbose) {
              message("Skipping ", mf, " (not a dataset)")
            }
            next
          }
          if (!mf %in% names(x = assay.group[['meta.features']])) {
            if (verbose) {
              message("Adding ", mf, " from scaled feature-level metadata")
            }
            assay.group[['meta.features']]$create_dataset(
              name = mf,
              dtype = source[['var']][[mf]]$get_type(),
              space = H5S$new(dims = assay.group[['features']]$dims)
            )
          } else if (verbose) {
            message("Merging ", mf, " from scaled feature-level metadata")
          }
          assay.group[['meta.features']][[mf]][features.use] <- source[['var']][[mf]]$read()
        }
      }
    } else if (!var_compound_handled) {
      message("Cannot yet add feature-level metadata from compound datasets")
      if (!assay.group$exists(name = 'meta.features')) {
        assay.group$create_group(name = 'meta.features')
      }
    }
  } else {
    if (inherits(x = source[['var']], what = 'H5Group')) {
      if (verbose) {
        message("Adding meta.features from var")
      }
      assay.group$obj_copy_from(
        src_loc = source,
        src_name = 'var',
        dst_name = 'meta.features'
      )
    } else if (!var_compound_handled) {
      # Only warn if we didn't already handle compound var with Python
      message("Cannot yet add feature-level metadata from compound datasets")
      if (!assay.group$exists(name = 'meta.features')) {
        assay.group$create_group(name = 'meta.features')
      }
    }
  }
  ColToFactor(dfgroup = assay.group[['meta.features']])
  # Sanitize column names and handle nullable dtypes in meta.features
  if (assay.group$exists(name = 'meta.features')) {
    SanitizeH5DFColumns(dfgroup = assay.group[['meta.features']])
  }
  mf_group <- assay.group[['meta.features']]
  if (mf_group$attr_exists(attr_name = 'column-order')) {
    colnames <- h5attr(x = mf_group, which = 'column-order')
    # Sanitize column names in the column-order attribute
    colnames <- vapply(colnames, SanitizeColumnName, character(1), USE.NAMES = FALSE)
    mf_group$create_attr(
      attr_name = 'colnames',
      robj = colnames,
      dtype = GuessDType(x = if (length(colnames) > 0) colnames[1] else colnames)
    )
  }
  if (inherits(x = source[['var']], what = 'H5Group')) {
    rownames_to_delete <- GetRownames(dset = 'var')
    # Only delete if the link exists in meta.features
    if (assay.group[['meta.features']]$exists(name = rownames_to_delete)) {
      assay.group[['meta.features']]$link_delete(name = rownames_to_delete)
    } else if (verbose) {
      message("Note: Rownames '", rownames_to_delete, "' not found in meta.features, skipping deletion")
    }
  }

  # Extract variable features from highly_variable column if available
  # Only check if var is an H5Group (compound datasets don't support sub-paths)
  if (inherits(x = source[['var']], what = 'H5Group') &&
      source$exists(name = 'var/highly_variable')) {
    if (verbose) {
      message("Converting highly_variable to variable.features")
    }

    tryCatch({
      # Check if highly_variable is a categorical (H5Group with categories/codes)
      # CellxGene stores as categorical with string "True"/"False" values
      hv_obj <- source[['var/highly_variable']]

      # Check object type
      is_h5group <- inherits(hv_obj, 'H5Group')
      has_categories <- is_h5group && hv_obj$exists('categories')
      has_codes <- is_h5group && hv_obj$exists('codes')

      if (is_h5group && has_categories && has_codes) {
        # Decode categorical: categories[codes + 1]
        categories <- hv_obj[['categories']][]
        codes <- hv_obj[['codes']][]
        decoded_values <- categories[codes + 1L]
        # Convert string "True"/"False" to logical
        highly_variable <- decoded_values == "True"
        if (verbose) {
          message("Decoded categorical highly_variable: ", sum(highly_variable), " of ", length(highly_variable), " are True")
        }
      } else if (inherits(hv_obj, 'H5D')) {
        # Direct read for non-categorical H5D dataset
        highly_variable <- hv_obj[]

        # Convert to logical if needed
        if (is.numeric(highly_variable)) {
          highly_variable <- as.logical(highly_variable)
        }
      } else {
        # Unsupported type
        if (verbose) {
          message("highly_variable is not a recognized type (H5Group categorical or H5D): ", class(hv_obj)[1])
        }
        highly_variable <- logical(0)
      }

      # Get feature names from /var (same dimension as highly_variable)
      # Use the same feature names as the Seurat object: prioritize feature_name (gene symbols)
      categorical_features <- ReadCategoricalFeatures(dset = 'var')
      if (!is.null(categorical_features)) {
        all_features <- categorical_features
      } else {
        # Fallback to index dataset (Ensembl IDs) or compound labels
        rownames_dset <- GetRownames(dset = 'var')
        if (source[['var']]$exists(rownames_dset)) {
          all_features <- source[[paste('var', rownames_dset, sep = '/')]][]
        } else {
          all_features <- source[['var']]$get_cpd_labels()
        }
      }

      # Get variable features
      if (is.logical(highly_variable) && length(highly_variable) == length(all_features)) {
        variable_features <- all_features[highly_variable]

        # Write to assay group
        if (length(variable_features) > 0) {
          assay.group$create_dataset(
            name = 'variable.features',
            robj = variable_features,
            dtype = GuessDType(x = variable_features)
          )
          if (verbose) {
            message("Added ", length(variable_features), " variable features")
          }
        }
      } else if (verbose) {
        message("Warning: highly_variable format incompatible with features (logical: ",
                is.logical(highly_variable), ", length: ", length(highly_variable),
                ", all_features length: ", length(all_features), ")")
      }
    }, error = function(e) {
      if (verbose) {
        message("Could not convert highly_variable to variable.features: ", conditionMessage(e))
      }
    })
  }

  # Add cell-level metadata
  if (source$exists(name = 'obs') && inherits(x = source[['obs']], what = 'H5Group')) {
    if (!source[['obs']]$exists(name = '__categories') && !getOption(x = "scConvert.dtypes.dataframe_as_group", default = TRUE)) {
      message("Conversion from H5AD to h5Seurat allowing compound datasets is not yet implemented")
    }
    dfile$obj_copy_from(
      src_loc = source,
      src_name = 'obs',
      dst_name = 'meta.data'
    )
    # Normalize h5ad categorical format (categories/codes) to h5Seurat format (levels/values)
    # h5ad uses 0-based indexing for codes, while R factors use 1-based indexing
    NormalizeH5ADCategorical <- function(dfgroup) {
      col_names <- names(dfgroup)
      # Skip internal/non-column entries
      col_names <- col_names[!col_names %in% c('__categories', '_index', 'index')]
      for (col_name in col_names) {
        col_obj <- dfgroup[[col_name]]
        # Check if this is an h5ad categorical group (has 'categories' and 'codes')
        if (!inherits(col_obj, 'H5Group')) next
        sub_names <- names(col_obj)
        if (!all(c('categories', 'codes') %in% sub_names)) next
        # Rename 'categories' to 'levels'
        if (!'levels' %in% sub_names) {
          col_obj$obj_copy_from(
            src_loc = col_obj,
            src_name = 'categories',
            dst_name = 'levels'
          )
          col_obj$link_delete(name = 'categories')
        }
        # Convert 'codes' to 'values' with 0-based to 1-based indexing conversion
        if (!'values' %in% sub_names) {
          codes_0based <- col_obj[['codes']]$read()
          codes_dtype <- col_obj[['codes']]$get_type()
          values_1based <- codes_0based + 1L
          values_1based[codes_0based == -1L] <- NA_integer_
          col_obj$create_dataset(
            name = 'values',
            robj = values_1based,
            dtype = codes_dtype
          )
          col_obj$link_delete(name = 'codes')
        }
      }
    }
    NormalizeH5ADCategorical(dfgroup = dfile[['meta.data']])
    ColToFactor(dfgroup = dfile[['meta.data']])
    SanitizeH5DFColumns(dfgroup = dfile[['meta.data']])
    md_group <- dfile[['meta.data']]
    if (md_group$attr_exists(attr_name = 'column-order')) {
      colnames <- h5attr(x = md_group, which = 'column-order')
      # Sanitize column names in the column-order attribute
      colnames <- vapply(colnames, SanitizeColumnName, character(1), USE.NAMES = FALSE)
      md_group$create_attr(
        attr_name = 'colnames',
        robj = colnames,
        dtype = GuessDType(x = if (length(colnames) > 0) colnames[1] else colnames)
      )
    }
    rownames <- GetRownames(dset = 'obs')
    dfile$obj_copy_from(
      src_loc = dfile,
      src_name = paste0('meta.data/', rownames),
      dst_name = 'cell.names'
    )
    dfile[['meta.data']]$link_delete(name = rownames)
  } else if (source$exists(name = 'obs') && inherits(x = source[['obs']], what = 'H5D')) {
    # Handle compound obs dataset using Python
    if (verbose) {
      message("Detected compound obs dataset; attempting to read with Python")
    }

    # Get the source file path
    h5file_path <- if (inherits(x = source, what = 'H5File')) {
      source$filename
    } else {
      as.character(source)
    }

    obs_data <- ReadCompoundDataset(h5file_path, 'obs')

    if (!is.null(obs_data) && nrow(obs_data) > 0) {
      # Successfully read compound dataset with Python
      if (verbose) {
        message("  Successfully read compound obs dataset with ", nrow(obs_data), " cells")
      }

      # Extract cell names from the 'index' column or first column
      cell_names <- if ('index' %in% names(obs_data)) {
        as.character(obs_data$index)
      } else {
        as.character(obs_data[[1]])
      }

      # Create meta.data group and add columns
      dfile$create_group(name = 'meta.data')
      # Sanitize column names for obs metadata
      obs_col_name_map <- SanitizeColumnNames(names(obs_data))
      for (col_name in names(obs_data)) {
        if (col_name == 'index' || col_name == names(obs_data)[1]) {
          next  # Skip the index column as it's used for cell names
        }
        col_data <- obs_data[[col_name]]
        sanitized_name <- obs_col_name_map[[col_name]]
        dfile[['meta.data']]$create_dataset(
          name = sanitized_name,
          robj = col_data,
          dtype = GuessDType(x = if (length(col_data) > 0) col_data[1] else col_data)
        )
      }
      # Add rownames to meta.data
      dfile[['meta.data']]$create_dataset(
        name = '_index',
        robj = cell_names,
        dtype = GuessDType(x = cell_names[1])
      )
      dfile[['meta.data']]$create_attr(
        attr_name = '_index',
        robj = '_index',
        dtype = CachedGuessDType(x = '_index')
      )
      dfile[['meta.data']]$create_attr(
        attr_name = 's4class',
        robj = 'DFrame',
        dtype = StringType()
      )

      # Create cell.names dataset
      dfile$create_dataset(
        name = 'cell.names',
        robj = cell_names,
        dtype = GuessDType(x = cell_names[1])
      )
    } else {
      # Fallback if Python reading fails
      message("Could not read compound obs dataset with Python; creating fake cell names")
      CreateFakeCellNames()
    }
  } else {
    message("No cell-level metadata present, creating fake cell names")
    CreateFakeCellNames()
  }

  # Restore Seurat-specific metadata from uns['seurat'] if available
  if (source$exists(name = 'uns') && source[['uns']]$exists(name = 'seurat')) {
    if (verbose) {
      message("Restoring Seurat-specific metadata from uns['seurat']")
    }

    seurat_group <- source[['uns/seurat']]

    # Check if this file was originally from Seurat
    if (seurat_group$exists(name = 'version')) {
      version_info <- seurat_group[['version']][]
      if (verbose) {
        message("  Found Seurat metadata from ", version_info)
      }
    }

    # Restore original assay name if different from current
    if (seurat_group$exists(name = 'assay_name')) {
      original_assay <- seurat_group[['assay_name']][]
      if (original_assay != assay && verbose) {
        message("  Original assay name was '", original_assay, "', now using '", assay, "'")
      }
    }

    # Note: Graph names are stored for reference, but graphs themselves
    # will be restored from uns/neighbors/distances in the existing code
    if (seurat_group$exists(name = 'graph_names') && verbose) {
      graph_names <- seurat_group[['graph_names']][]
      message("  Original graph names: ", paste(graph_names, collapse = ", "))
    }
  }

  # Add dimensional reduction information
  if (source$exists(name = 'obsm')) {
    # Add cell embeddings
    if (inherits(x = source[['obsm']], what = 'H5Group')) {
      for (reduc in names(x = source[['obsm']])) {
        sreduc <- gsub(pattern = '^X_', replacement = '', x = reduc)

        # Skip spatial - it will be stored in boundaries/centroids for VisiumV2
        if (sreduc == 'spatial') {
          message("Skipping ", reduc, " - spatial data stored in boundaries/centroids")
          next
        }

        reduc.group <- dfile[['reductions']]$create_group(name = sreduc)
        message("Adding ", reduc, " as cell embeddings for ", sreduc)
        # Read embeddings and fix orientation if needed.
        # Python h5ad obsm shape is (n_obs, n_dims) but hdf5r may read it
        # as (n_dims, n_obs) due to C vs Fortran order differences.
        obsm_item <- source[['obsm']][[reduc]]
        if (inherits(obsm_item, "H5D")) {
          emb_data <- obsm_item$read()
        } else {
          emb_data <- as.matrix(ReadSparseMatrix(obsm_item))
        }
        n_cells <- dfile[['cell.names']]$dims
        if (nrow(emb_data) != n_cells && ncol(emb_data) == n_cells) {
          emb_data <- t(emb_data)
        }
        reduc.group$create_dataset(
          name = 'cell.embeddings',
          robj = emb_data,
          dtype = GuessDType(x = emb_data[1, 1])
        )
        reduc.group$create_group(name = 'misc')
        reduc.group$create_attr(
          attr_name = 'active.assay',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
        key <- paste0(
          if (grepl(pattern = 'pca', x = sreduc, ignore.case = TRUE)) {
            'PC'
          } else if (grepl(pattern = 'tsne', x = sreduc, ignore.case = TRUE)) {
            'tSNE'
          } else {
            sreduc
          },
          '_'
        )
        reduc.group$create_attr(
          attr_name = 'key',
          robj = key,
          dtype = GuessDType(x = reduc)
        )
        global <- BoolToInt(x = grepl(
          pattern = 'tsne|umap',
          x = sreduc,
          ignore.case = TRUE
        ))
        reduc.group$create_attr(
          attr_name = 'global',
          robj = global,
          dtype = GuessDType(x = global)
        )
      }
    } else {
      message("Reading compound dimensional reductions not yet supported, please update your H5AD file")
    }
    # Add feature loadings
    if (source$exists(name = 'varm')) {
      if (inherits(x = source[['varm']], what = 'H5Group')) {
        for (reduc in names(x = source[['varm']])) {
          sreduc <- switch(EXPR = reduc, 'PCs' = 'pca', tolower(x = reduc))
          if (!isTRUE(x = sreduc %in% names(x = dfile[['reductions']]))) {
            message("Cannot find a reduction named ", sreduc, " (", reduc, " in varm)")
            next
          }
          if (isTRUE(x = verbose)) {
            message("Adding ", reduc, " as feature loadings fpr ", sreduc)
          }
          scTranspose(
            x = source[['varm']][[reduc]],
            dest = dfile[['reductions']][[sreduc]],
            dname = 'feature.loadings',
            verbose = FALSE
          )
          reduc.features <- dfile[['reductions']][[sreduc]][['feature.loadings']]$dims[1]
          assay.features <- if (prod(assay.group[['features']]$dims) == reduc.features) {
            'features'
          } else if (assay.group$exists(name = 'scaled.features') && prod(assay.group[['scaled.features']]$dims) == reduc.features) {
            'scaled.features'
          } else {
            NULL
          }
          if (is.null(x = assay.features)) {
            message("Cannot find features for feature loadings, will not be able to load")
          } else {
            dfile[['reductions']][[sreduc]]$obj_copy_from(
              src_loc = assay.group,
              src_name = assay.features,
              dst_name = 'features'
            )
          }
        }
      } else {
        message("Reading compound dimensional reductions not yet supported")
      }
    }
    # Add miscellaneous information
    if (source$exists(name = 'uns')) {
      for (reduc in names(x = source[['uns']])) {
        if (!isTRUE(x = reduc %in% names(x = dfile[['reductions']]))) {
          next
        }
        if (verbose) {
          message("Adding miscellaneous information for ", reduc)
        }
        dfile[['reductions']][[reduc]]$link_delete(name = 'misc')
        dfile[['reductions']][[reduc]]$obj_copy_from(
          src_loc = source[['uns']],
          src_name = reduc,
          dst_name = 'misc'
        )
        if ('variance' %in% names(x = dfile[['reductions']][[reduc]][['misc']])) {
          if (verbose) {
            message("Adding standard deviations for ", reduc)
          }
          dfile[['reductions']][[reduc]]$create_dataset(
            name = 'stdev',
            robj = sqrt(x = dfile[['reductions']][[reduc]][['misc']][['variance']][]),
            dtype = GuessDType(x = 1.0)
          )
        }
      }
    }
  }
  # Add project and cell identities
  Project(object = dfile) <- 'AnnData'
  idents <- dfile$create_group(name = 'active.ident')
  idents$create_dataset(
    name = 'values',
    dtype = GuessDType(x = 1L),
    space = H5S$new(dims = dfile[['cell.names']]$dims)
  )
  idents$create_dataset(
    name = 'levels',
    robj = 'AnnData',
    dtype = CachedGuessDType(x = 'AnnData')
  )
  idents[['values']]$write(
    args = list(seq.default(from = 1, to = idents[['values']]$dims)),
    value = 1L
  )
  # Add nearest-neighbor graphs
  # In AnnData: uns/neighbors/distances = k-NN distance matrix
  #             uns/neighbors/connectivities = weighted adjacency/connectivity matrix (like SNN)
  # In Seurat:  assay_nn = nearest neighbor graph (distances)
  #             assay_snn = shared nearest neighbor graph (connectivities)

  if (Exists(x = source, name = 'uns/neighbors')) {
    # Save distances as NN graph
    if (source[['uns/neighbors']]$exists('distances')) {
      nn.graph.name <- paste(assay, 'nn', sep = '_')
      if (verbose) {
        message("Saving nearest-neighbor distances as ", nn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'uns/neighbors/distances',
        dst_name = nn.graph.name
      )
      FixGraphAttrs(dfile[['graphs']][[nn.graph.name]], assay)
    }

    # Save connectivities as SNN graph
    if (source[['uns/neighbors']]$exists('connectivities')) {
      snn.graph.name <- paste(assay, 'snn', sep = '_')
      if (verbose) {
        message("Saving nearest-neighbor connectivities as ", snn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'uns/neighbors/connectivities',
        dst_name = snn.graph.name
      )
      FixGraphAttrs(dfile[['graphs']][[snn.graph.name]], assay)
    }
  }

  # Also check obsp for neighbor graphs (scanpy default location)
  # This handles h5ad files where neighbors are stored directly in obsp
  # rather than in uns/neighbors
  if (source$exists(name = 'obsp')) {
    obsp_names <- names(x = source[['obsp']])
    nn.graph.name <- paste(assay, 'nn', sep = '_')
    snn.graph.name <- paste(assay, 'snn', sep = '_')

    # Transfer distances as NN graph if not already done from uns/neighbors
    if ('distances' %in% obsp_names && !dfile[['graphs']]$exists(nn.graph.name)) {
      if (verbose) {
        message("Saving obsp/distances as ", nn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'obsp/distances',
        dst_name = nn.graph.name
      )
      FixGraphAttrs(dfile[['graphs']][[nn.graph.name]], assay)
    }

    # Transfer connectivities as SNN graph if not already done from uns/neighbors
    if ('connectivities' %in% obsp_names && !dfile[['graphs']]$exists(snn.graph.name)) {
      if (verbose) {
        message("Saving obsp/connectivities as ", snn.graph.name)
      }
      dfile[['graphs']]$obj_copy_from(
        src_loc = source,
        src_name = 'obsp/connectivities',
        dst_name = snn.graph.name
      )
      FixGraphAttrs(dfile[['graphs']][[snn.graph.name]], assay)
    }

    # Transfer any remaining obsp graphs not already handled
    handled <- c('distances', 'connectivities')
    for (oname in setdiff(obsp_names, handled)) {
      if (!dfile[['graphs']]$exists(oname)) {
        if (verbose) message("Saving obsp/", oname, " as graph")
        dfile[['graphs']]$obj_copy_from(
          src_loc = source,
          src_name = paste0('obsp/', oname),
          dst_name = oname
        )
        FixGraphAttrs(dfile[['graphs']][[oname]], assay)
      }
    }
  }

  # Add miscellaneous information
  if (source$exists(name = 'uns')) {
    misc <- setdiff(
      x = names(x = source[['uns']]),
      y = c('neighbors', names(x = dfile[['reductions']]))
    )
    for (i in misc) {
      if (verbose) {
        message("Adding ", i, " to miscellaneous data")
      }
      tryCatch({
        dfile[['misc']]$obj_copy_from(
          src_loc = source[['uns']],
          src_name = i,
          dst_name = i
        )
      }, error = function(e) {
        if (verbose) {
          message("  Warning: Could not copy ", i, " to miscellaneous data: ", conditionMessage(e))
        }
      })
    }
  }
  # Add spatial images from uns/spatial
  if (source$exists(name = 'uns') && source[['uns']]$exists(name = 'spatial')) {
    if (verbose) {
      message("Adding spatial images from uns/spatial")
    }
    spatial_uns <- source[['uns/spatial']]
    library_ids <- names(x = spatial_uns)

    # Detect coordinate format early to determine if image needs transformation
    # Native h5ad (scanpy/CellXGene): shape is (n_cells, 2)
    # h5ad from Seurat conversion: shape is (2, n_cells)
    is_native_h5ad <- FALSE
    spatial_key <- if (source$exists(name = 'obsm/X_spatial')) {
      'obsm/X_spatial'
    } else if (source$exists(name = 'obsm/spatial')) {
      'obsm/spatial'
    } else {
      NULL
    }
    if (!is.null(spatial_key)) {
      coords_dims <- dim(source[[spatial_key]]$read())
      is_native_h5ad <- coords_dims[2] == 2  # (n_cells, 2) format
      if (verbose) {
        message("  Detected ", if(is_native_h5ad) "native h5ad" else "Seurat-converted", " coordinate format")
      }
    }

    for (lib_id in library_ids) {
      lib_group <- spatial_uns[[lib_id]]

      # Skip non-group items (e.g., 'is_single' dataset)
      if (!inherits(lib_group, 'H5Group')) {
        next
      }

      # Check if images exist for this library
      if (lib_group$exists(name = 'images')) {
        images_group <- lib_group[['images']]
        image_types <- names(x = images_group)

        # Create image group in h5seurat (using library_id as image name)
        # Remove illegal characters from library_id for use as image name
        img_name <- gsub(pattern = '[^A-Za-z0-9_]', replacement = '_', x = lib_id)

        if (verbose) {
          message("  Processing library '", lib_id, "' as image '", img_name, "'")
        }

        # Create image group if it doesn't exist
        if (!dfile[['images']]$exists(name = img_name)) {
          img_h5 <- dfile[['images']]$create_group(name = img_name)

          # Add required attributes for image indexing
          # assay attribute
          img_h5$create_attr(
            attr_name = 'assay',
            robj = assay,  # Use the assay parameter
            dtype = GuessDType(x = assay)
          )

          # s4class attribute (VisiumV2 for spatial data)
          img_h5$create_attr(
            attr_name = 's4class',
            robj = 'VisiumV2',
            dtype = CachedGuessDType(x = 'VisiumV2')
          )

          # global attribute (make image globally available)
          img_h5$create_attr(
            attr_name = 'global',
            robj = 1L,  # Use integer 1 for TRUE
            dtype = GuessDType(x = 1L)
          )
        } else {
          img_h5 <- dfile[['images']][[img_name]]
        }

        # Add key dataset (required by Seurat)
        if (!img_h5$exists(name = 'key')) {
          img_h5$create_dataset(
            name = 'key',
            robj = paste0(img_name, '_'),
            dtype = GuessDType(x = paste0(img_name, '_'))
          )
        }

        # Store lowres as the primary 'image' dataset (used for visualization)
        # Note: hires can be added later if needed
        if ('lowres' %in% image_types) {
          if (verbose) {
            message("    Adding lowres as primary image")
          }

          # Read lowres image data from h5ad
          img_data <- images_group[['lowres']]$read()

          # Check if this is a 3D array
          if (length(dim(img_data)) == 3L) {
            # h5ad/anndata stores images as (channels, width, height) in HDF5
            # Seurat expects (height, width, channels)
            # aperm(c(3,2,1)) converts (C,W,H) -> (H,W,C)
            img_arr <- aperm(img_data, c(3L, 2L, 1L))

            # Write to h5seurat as 'image' dataset
            if (img_h5$exists(name = 'image')) {
              img_h5$link_delete(name = 'image')
            }
            img_h5$create_dataset(
              name = 'image',
              robj = img_arr,
              dtype = GuessDType(x = img_arr)
            )
          } else {
            if (verbose) {
              message("      Warning: Lowres image is not a 3D array, skipping")
            }
          }
        } else if ('hires' %in% image_types) {
          # Fallback to hires if lowres not available
          if (verbose) {
            message("    Adding hires as primary image (lowres not available)")
          }

          img_data <- images_group[['hires']]$read()
          if (length(dim(img_data)) == 3L) {
            # h5ad/anndata stores images as (channels, width, height) in HDF5
            # Seurat expects (height, width, channels)
            # aperm(c(3,2,1)) converts (C,W,H) -> (H,W,C)
            img_arr <- aperm(img_data, c(3L, 2L, 1L))
            if (img_h5$exists(name = 'image')) {
              img_h5$link_delete(name = 'image')
            }
            img_h5$create_dataset(
              name = 'image',
              robj = img_arr,
              dtype = GuessDType(x = img_arr)
            )
          }
        }

        # Add scale factors if available
        if (lib_group$exists(name = 'scalefactors')) {
          if (verbose) {
            message("    Adding scale factors")
          }

          # Create scale.factors group
          if (img_h5$exists(name = 'scale.factors')) {
            img_h5$link_delete(name = 'scale.factors')
          }
          sf_group <- img_h5$create_group(name = 'scale.factors')

          # Map from h5ad names to h5seurat names
          sf_map <- c(
            tissue_hires_scalef = 'hires',
            tissue_lowres_scalef = 'lowres',
            spot_diameter_fullres = 'spot',
            fiducial_diameter_fullres = 'fiducial'
          )

          sf_src <- lib_group[['scalefactors']]
          for (sf_h5ad_name in names(x = sf_src)) {
            # Find corresponding h5seurat name from the map
            # sf_map has h5ad names as keys and h5seurat names as values
            if (sf_h5ad_name %in% names(sf_map)) {
              sf_h5seurat_name <- unname(sf_map[sf_h5ad_name])
            } else {
              # If no mapping found, use the h5ad name directly
              sf_h5seurat_name <- sf_h5ad_name
            }

            sf_value <- sf_src[[sf_h5ad_name]][]

            # Convert diameter to radius for spot and fiducial
            if (sf_h5ad_name %in% c('spot_diameter_fullres', 'fiducial_diameter_fullres')) {
              sf_value <- sf_value / 2.0
            }

            sf_group$create_dataset(
              name = sf_h5seurat_name,
              robj = sf_value,
              dtype = GuessDType(x = sf_value)
            )
          }

          # If lowres scale factor is missing but hires exists, set lowres = hires
          # This is needed when h5ad only has hires image (no lowres)
          # Seurat's SpatialFeaturePlot defaults to using lowres scale
          if (!sf_group$exists('lowres') && sf_group$exists('hires')) {
            hires_val <- sf_group[['hires']][]
            sf_group$create_dataset(
              name = 'lowres',
              robj = hires_val,
              dtype = GuessDType(x = hires_val)
            )
            if (verbose) {
              message("      Setting lowres scale factor = hires (", hires_val, ")")
            }
          }
        }

        # Create boundaries/centroids structure for VisiumV2
        # This requires spatial coordinates from obsm/X_spatial or obsm/spatial
        spatial_key <- if (source$exists(name = 'obsm/X_spatial')) {
          'obsm/X_spatial'
        } else if (source$exists(name = 'obsm/spatial')) {
          'obsm/spatial'
        } else {
          NULL
        }

        if (!is.null(spatial_key)) {
          if (verbose) {
            message("    Creating boundaries/centroids for VisiumV2 (using ", spatial_key, ")")
          }

          # Read spatial coordinates from h5ad
          # Use $read() instead of [] to avoid HDF5 indexing issues
          coords_h5ad <- source[[spatial_key]]$read()

          # Detect coordinate orientation:
          # - Native h5ad (scanpy/CellXGene): shape is (n_cells, 2) where rows are cells
          # - h5ad from Seurat conversion: shape is (2, n_cells) where columns are cells
          # We detect by checking which dimension is 2
          coords_dims <- dim(coords_h5ad)
          is_native_h5ad <- coords_dims[2] == 2  # (n_cells, 2) format

          if (verbose) {
            message("    Coordinate array shape: ", paste(coords_dims, collapse=" x "),
                    " (", if(is_native_h5ad) "native h5ad format" else "Seurat-converted format", ")")
          }

          # Get cell names from the h5Seurat file (already loaded earlier)
          cell_names <- as.character(dfile[['cell.names']][])

          # Get library IDs to filter cells for this library
          # Try common library ID column names
          lib_id_col <- NULL
          lib_ids_vec <- NULL
          for (col_name in c('sangerID', 'library_id', 'sample', 'batch')) {
            if (source$exists(paste0('obs/', col_name))) {
              lib_id_obj <- source[[paste0('obs/', col_name)]]

              # Handle categorical/factor structure
              if (inherits(lib_id_obj, 'H5Group')) {
                if (lib_id_obj$exists('codes') && lib_id_obj$exists('categories')) {
                  codes <- lib_id_obj[['codes']][]
                  categories <- lib_id_obj[['categories']][]
                  lib_ids_vec <- categories[codes + 1]  # Convert to 1-indexed
                  lib_id_col <- col_name
                  break
                }
              } else {
                lib_ids_vec <- as.character(lib_id_obj[])
                lib_id_col <- col_name
                break
              }
            }
          }

          # Filter cells for this library
          if (!is.null(lib_ids_vec)) {
            cell_mask <- lib_ids_vec == lib_id
            if (sum(cell_mask) == 0) {
              if (verbose) {
                message("    Warning: No cells found for library ", lib_id, ", skipping centroids")
              }
            } else {
              # Filter coordinates based on format
              if (is_native_h5ad) {
                # Native format (n_cells, 2): filter rows
                filtered_coords <- coords_h5ad[cell_mask, , drop = FALSE]
              } else {
                # Seurat-converted format (2, n_cells): filter columns
                filtered_coords <- coords_h5ad[, cell_mask, drop = FALSE]
              }
              filtered_cells <- cell_names[cell_mask]

              if (verbose) {
                message("    Filtered to ", length(filtered_cells), " cells for library ", lib_id)
              }

              # Convert to (n, 2) format for centroids
              if (is_native_h5ad) {
                # Native h5ad: already in (n, 2) format
                coords_matrix <- filtered_coords
              } else {
                # Seurat-converted: transpose from (2, n) to (n, 2)
                coords_matrix <- t(filtered_coords)
              }

              # Note: Coordinates are kept in full resolution
              # Seurat's GetTissueCoordinates() will handle scaling via scale.factors
              # when called with scale = "lowres" or scale = "hires"

              # Calculate radius from scale factors (in full resolution)
              radius <- if (img_h5$exists('scale.factors/spot')) {
                img_h5[['scale.factors/spot']][]
              } else {
                # Default Visium spot diameter is ~143 pixels in full res
                # (based on spot_diameter_fullres from the data)
                71.5  # radius = diameter / 2
              }

              # Create boundaries group
              if (!img_h5$exists('boundaries')) {
                boundaries_group <- img_h5$create_group('boundaries')
              } else {
                boundaries_group <- img_h5[['boundaries']]
              }

              # Create centroids group
              if (boundaries_group$exists('centroids')) {
                boundaries_group$link_delete('centroids')
              }
              centroids_group <- boundaries_group$create_group('centroids')

              # Store coords (n x 2 matrix)
              centroids_group$create_dataset(
                name = 'coords',
                robj = coords_matrix,
                dtype = GuessDType(x = coords_matrix)
              )

              # Store cell names
              centroids_group$create_dataset(
                name = 'cells',
                robj = filtered_cells,
                dtype = GuessDType(x = filtered_cells)
              )

              # Store radius (spot radius for visualization)
              centroids_group$create_dataset(
                name = 'radius',
                robj = radius,
                dtype = GuessDType(x = radius)
              )

              # Store nsides (0 = infinite sides/circle, returns centroids not vertices)
              centroids_group$create_dataset(
                name = 'nsides',
                robj = 0L,
                dtype = GuessDType(x = 0L)
              )

              # Store theta (rotation angle, 0 for Visium)
              centroids_group$create_dataset(
                name = 'theta',
                robj = 0.0,
                dtype = GuessDType(x = 0.0)
              )

              if (verbose) {
                message("    Created centroids with ", nrow(coords_matrix), " spots")
              }
            }
          } else {
            if (verbose) {
              message("    Warning: No library ID column found, using all cells for centroids")
            }

            # Use all cells if no library ID found
            # Convert to (n, 2) format for centroids based on detected format
            if (is_native_h5ad) {
              # Native h5ad: already in (n, 2) format
              coords_matrix <- coords_h5ad
            } else {
              # Seurat-converted: transpose from (2, n) to (n, 2)
              coords_matrix <- t(coords_h5ad)
            }

            radius <- if (img_h5$exists('scale.factors/spot')) {
              img_h5[['scale.factors/spot']][]
            } else {
              71.5  # Default radius in full resolution
            }

            if (!img_h5$exists('boundaries')) {
              boundaries_group <- img_h5$create_group('boundaries')
            } else {
              boundaries_group <- img_h5[['boundaries']]
            }

            if (boundaries_group$exists('centroids')) {
              boundaries_group$link_delete('centroids')
            }
            centroids_group <- boundaries_group$create_group('centroids')

            centroids_group$create_dataset(
              name = 'coords',
              robj = coords_matrix,
              dtype = GuessDType(x = coords_matrix)
            )

            centroids_group$create_dataset(
              name = 'cells',
              robj = cell_names,
              dtype = GuessDType(x = cell_names)
            )

            centroids_group$create_dataset(
              name = 'radius',
              robj = radius,
              dtype = GuessDType(x = radius)
            )

            centroids_group$create_dataset(
              name = 'nsides',
              robj = 0L,
              dtype = GuessDType(x = 0L)
            )

            centroids_group$create_dataset(
              name = 'theta',
              robj = 0.0,
              dtype = GuessDType(x = 0.0)
            )
          }
        }
      }
    }
  }
  # Add layers to the RNA assay (V5 style)
  if (Exists(x = source, name = 'layers')) {
    # Get the RNA assay
    rna.assay <- dfile[['assays']][[assay]]

    # For V5, we add layers as named slots within the assay
    # Common layer names: counts, data, scale.data
    for (layer in names(x = source[['layers']])) {
      # Map layer names to appropriate slots
      # Check if this layer should be added (avoid duplicates)
      layer.name <- layer

      # Skip if this slot already exists (was filled from raw/X or X)
      if (rna.assay$exists(name = layer.name)) {
        if (verbose) {
          message("Skipping layer ", layer, " (already exists in RNA assay)")
        }
        next
      }

      if (verbose) {
        message("Adding layer ", layer, " to RNA assay")
      }

      # Copy the layer data to the appropriate slot in RNA assay
      rna.assay$obj_copy_from(
        src_loc = source[['layers']],
        src_name = layer,
        dst_name = layer.name
      )

      # Fix dimensions attributes if needed
      if (isTRUE(x = AttrExists(x = rna.assay[[layer.name]], name = 'shape'))) {
        dims <- rev(x = h5attr(x = rna.assay[[layer.name]], which = 'shape'))
        rna.assay[[layer.name]]$create_attr(
          attr_name = 'dims',
          robj = dims,
          dtype = GuessDType(x = dims)
        )
        rna.assay[[layer.name]]$attr_delete(attr_name = 'shape')
      }
    }
  }
  return(dfile)
}

#' Convert h5Seurat files to H5AD files
#'
#' @inheritParams scConvert
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}}
#' object
#'
#' @section h5Seurat to AnnData/H5AD:
#' The h5Seurat to AnnData/H5AD conversion will try to automatically fill in
#' datasets based on data presence. Data presense is determined by the h5Seurat
#' index (\code{source$index()}). It works in the following manner:
#' \subsection{Assay data}{
#'  \itemize{
#'   \item \code{X} will be filled with \code{scale.data} if \code{scale.data}
#'   is present; otherwise, it will be filled with \code{data}
#'   \item \code{var} will be filled with \code{meta.features} \strong{only} for
#'   the features present in \code{X}; for example, if \code{X} is filled with
#'   \code{scale.data}, then \code{var} will contain only features that have
#'   been scaled
#'   \item \code{raw.X} will be filled with \code{data} if \code{X} is filled
#'   with \code{scale.data}; otherwise, it will be filled with \code{counts}. If
#'   \code{counts} is not present, then \code{raw} will not be filled
#'   \item \code{raw.var} will be filled with \code{meta.features} with the
#'   features present in \code{raw.X}; if \code{raw.X} is not filled, then
#'   \code{raw.var} will not be filled
#'  }
#' }
#' \subsection{Cell-level metadata}{
#'  Cell-level metadata is added to \code{obs}
#' }
#' \subsection{Dimensional reduction information}{
#'  Only dimensional reductions associated with \code{assay} or marked as
#'  \link[SeuratObject:IsGlobal]{global} will be transfered to the H5AD file. For
#'  every reduction \code{reduc}:
#'  \itemize{
#'   \item cell embeddings are placed in \code{obsm} and renamed to
#'   \code{X_reduc}
#'   \item feature loadings, if present, are placed in \code{varm} and renamed
#'   to either \dQuote{PCs} if \code{reduc} is \dQuote{pca} otherwise
#'   \code{reduc} in all caps
#'  }
#'  For example, if \code{reduc} is \dQuote{ica}, then cell embeddings will be
#'  \dQuote{X_ica} in \code{obsm} and feature loaodings, if present, will be
#'  \dQuote{ICA} in \code{varm}
#' }
#' \subsection{Nearest-neighbor graphs}{
#'  If a nearest-neighbor graph is associated with \code{assay}, it will be
#'  added to \code{uns/neighbors/distances}; if more than one graph is present,
#'  then \strong{only} the last graph according to the index will be added.
#' }
#' \subsection{Layers}{
#'  Data from other assays can be added to \code{layers} if they have the same
#'  shape as \code{X} (same number of cells and features). To determine this,
#'  the shape of each alternate assays's \code{scale.data} and \code{data} slots
#'  are determined. If they are the same shape as \code{X}, then that slot
#'  (\code{scale.data} is given priority over \code{data}) will be added as a
#'  layer named the name of the assay (eg. \dQuote{SCT}). In addition, the
#'  features names will be added to \code{var} as \code{assay_features}
#'  (eg. \dQuote{SCT_features}).
#' }
#'
#' @keywords internal
#'
H5SeuratToH5AD <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE
) {
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination H5AD file exists", call. = FALSE)
    }
  }
  rownames <- '_index'
  dfile <- H5File$new(filename = dest, mode = WriteMode(overwrite = FALSE))

  # Define mapping tables for Seurat to scanpy naming conventions
  if (standardize) {
    if (verbose) {
      message("Standardization enabled: will convert Seurat names to scanpy conventions")
    }
    # Cell-level metadata (obs) mappings
    obs_name_map <- c(
      "nCount_RNA" = "n_counts",
      "nFeature_RNA" = "n_genes",
      "percent.mt" = "percent_mito",
      "seurat_clusters" = "clusters"
    )

    # Feature-level metadata (var) mappings
    var_name_map <- c(
      "vst.variable" = "highly_variable",
      "vst.mean" = "means",
      "vst.variance" = "dispersions",
      "vst.variance.standardized" = "dispersions_norm",
      "vst.variance.expected" = "dispersions_expected"
    )
  } else {
    obs_name_map <- character(0)
    var_name_map <- character(0)
  }

  # Transfer data frames from h5Seurat files to H5AD files
  #
  # @param src Source dataset
  # @param dname Name of destination
  # @param index Integer values of rows to take
  #
  # @return Invisibly returns \code{NULL}
  #
  TransferDF <- function(src, dname, index) {
    # Helper function to apply name mapping based on destination type
    apply_name_map <- function(col_name, dest_type) {
      if (!standardize) return(col_name)

      name_map <- if (dest_type == "obs") obs_name_map else var_name_map

      # If column name has a mapping, use it; otherwise keep original
      if (col_name %in% names(name_map)) {
        return(unname(name_map[col_name]))
      }
      return(col_name)
    }

    # Handle V5 environment-wrapped metadata
    if (dname == "obs" && is.environment(x = src)) {
      # Extract from environment if possible
      if (exists("hgroup", envir = src, inherits = FALSE)) {
        src <- get("hgroup", envir = src, inherits = FALSE)
      } else if (source$exists(name = 'meta.data')) {
        # Direct metadata transfer for V5
        if (verbose) {
          message("Handling V5 metadata transfer for obs")
        }

        # Create obs group
        if (!dfile$exists(name = 'obs')) {
          dfile$create_group(name = 'obs')
        }

        # Transfer metadata columns
        meta_group <- source[['meta.data']]
        if (inherits(x = meta_group, what = 'H5Group')) {
          for (col in names(x = meta_group)) {
            if (col == '__categories' || col == '_index') next

            # Apply name mapping if standardize is enabled
            mapped_col <- apply_name_map(col, "obs")

            if (IsFactor(x = meta_group[[col]])) {
              # Use newer anndata format: each categorical is a group with categories/codes
              dfile[['obs']]$create_group(name = mapped_col)

              # Convert codes to 0-based indices, handling NA values
              # R factor codes are 1-based; h5ad uses 0-based with -1 for missing
              raw_codes <- meta_group[[col]][['values']][index]
              codes_values <- ifelse(is.na(raw_codes), -1L, raw_codes - 1L)
              n_categories <- length(meta_group[[col]][['levels']][])
              has_na <- any(is.na(raw_codes))

              # Choose appropriate integer dtype based on number of categories
              # Must use signed integers when NA values present (need -1)
              codes_dtype <- if (has_na || n_categories <= 127) {
                hdf5r::h5types$H5T_NATIVE_INT8
              } else if (n_categories <= 255) {
                hdf5r::h5types$H5T_NATIVE_UINT8
              } else if (n_categories <= 32767) {
                hdf5r::h5types$H5T_NATIVE_INT16
              } else if (n_categories <= 65535) {
                hdf5r::h5types$H5T_NATIVE_UINT16
              } else {
                hdf5r::h5types$H5T_NATIVE_INT32
              }

              # Add codes dataset with proper dtype for scanpy compatibility
              dfile[['obs']][[mapped_col]]$create_dataset(
                name = 'codes',
                robj = codes_values,
                dtype = codes_dtype
              )

              # Add encoding attributes to codes
              dfile[['obs']][[mapped_col]][['codes']]$create_attr(
                attr_name = 'encoding-type',
                robj = 'array',
                dtype = CachedGuessDType(x = 'array'),
                space = ScalarSpace()
              )
              dfile[['obs']][[mapped_col]][['codes']]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = CachedGuessDType(x = '0.2.0'),
                space = ScalarSpace()
              )

              # Add categories dataset
              dfile[['obs']][[mapped_col]]$create_dataset(
                name = 'categories',
                robj = as.character(meta_group[[col]][['levels']][]),
                dtype = StringType('utf8')
              )

              # Add encoding attributes to categories
              AddAnndataEncoding(dfile[['obs']][[mapped_col]][['categories']], encoding_type = 'string-array')

              # Add attributes to the categorical group
              dfile[['obs']][[mapped_col]]$create_attr(
                attr_name = 'encoding-type',
                robj = 'categorical',
                dtype = CachedGuessDType(x = 'categorical'),
                space = ScalarSpace()
              )
              dfile[['obs']][[mapped_col]]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = CachedGuessDType(x = '0.2.0'),
                space = ScalarSpace()
              )
              dfile[['obs']][[mapped_col]]$create_attr(
                attr_name = 'ordered',
                robj = FALSE,
                dtype = GuessDType(x = FALSE),
                space = ScalarSpace()
              )
            } else {
              # Regular column
              dfile[['obs']]$create_dataset(
                name = mapped_col,
                robj = if (is.null(index)) meta_group[[col]][] else meta_group[[col]][index],
                dtype = meta_group[[col]]$get_type()
              )
            }
          }

          # Add column order with mapped names (required by anndata >= 0.11)
          if (meta_group$attr_exists(attr_name = 'colnames')) {
            original_colnames <- h5attr(x = meta_group, which = 'colnames')
            if (standardize) {
              mapped_colnames <- sapply(original_colnames, function(cn) apply_name_map(cn, "obs"))
            } else {
              mapped_colnames <- original_colnames
            }
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = mapped_colnames,
              dtype = GuessDType(x = mapped_colnames)
            )
          } else {
            # Fallback: derive column order from group contents
            obs_col_names <- setdiff(dfile[['obs']]$names, c('_index', '__categories'))
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = if (length(obs_col_names) > 0) obs_col_names else character(0),
              dtype = if (length(obs_col_names) > 0) GuessDType(x = obs_col_names) else StringType('utf8')
            )
          }

          # Add encoding attributes
          encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
          names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
          for (i in seq_along(along.with = encoding.info)) {
            attr.name <- names(x = encoding.info)[i]
            if (!dfile[['obs']]$attr_exists(attr_name = attr.name)) {
              dfile[['obs']]$create_attr(
                attr_name = attr.name,
                robj = encoding.info[i],
                dtype = GuessDType(x = encoding.info[i]),
                space = ScalarSpace()
              )
            }
          }
        }
        return(invisible(x = NULL))
      }
    }

    # Validate source for standard transfer
    if (is.null(x = src) || !inherits(x = src, what = c('H5D', 'H5Group'))) {
      message("TransferDF: invalid source for ", dname)
      return(invisible(x = NULL))
    }

    if (verbose) {
      message("Transfering ", basename(path = src$get_obj_name()), " to ", dname)
    }

    if (inherits(x = src, what = 'H5D')) {
      # Create a name mapping function for this dname
      map_fn <- if (standardize) {
        function(col_name) apply_name_map(col_name, dname)
      } else {
        NULL
      }
      CompoundToGroup(
        src = src,
        dst = dfile,
        dname = dname,
        order = 'column-order',
        index = index,
        name_map_fn = map_fn
      )
    } else if (inherits(x = src, what = 'H5Group')) {
      dfile$create_group(name = dname)
      for (i in src$names) {
        # Apply name mapping if standardize is enabled
        mapped_i <- apply_name_map(i, dname)

        if (IsFactor(x = src[[i]])) {
          # Use newer anndata format: each categorical is a group with categories/codes
          dfile[[dname]]$create_group(name = mapped_i)

          # Convert codes to 0-based indices, handling NA values
          # R factor codes are 1-based; h5ad uses 0-based with -1 for missing
          raw_codes <- src[[i]][['values']][index]
          codes_values <- ifelse(is.na(raw_codes), -1L, raw_codes - 1L)
          n_categories <- length(src[[i]][['levels']][])
          has_na <- any(is.na(raw_codes))

          # Choose appropriate integer dtype based on number of categories
          # Must use signed integers when NA values present (need -1)
          codes_dtype <- if (has_na || n_categories <= 127) {
            hdf5r::h5types$H5T_NATIVE_INT8
          } else if (n_categories <= 255) {
            hdf5r::h5types$H5T_NATIVE_UINT8
          } else if (n_categories <= 32767) {
            hdf5r::h5types$H5T_NATIVE_INT16
          } else if (n_categories <= 65535) {
            hdf5r::h5types$H5T_NATIVE_UINT16
          } else {
            hdf5r::h5types$H5T_NATIVE_INT32
          }

          # Add codes dataset with proper dtype for scanpy compatibility
          dfile[[dname]][[mapped_i]]$create_dataset(
            name = 'codes',
            robj = codes_values,
            dtype = codes_dtype
          )

          # Add encoding attributes to codes
          dfile[[dname]][[mapped_i]][['codes']]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = CachedGuessDType(x = 'array'),
            space = ScalarSpace()
          )
          dfile[[dname]][[mapped_i]][['codes']]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = CachedGuessDType(x = '0.2.0'),
            space = ScalarSpace()
          )

          # Add categories dataset
          dfile[[dname]][[mapped_i]]$create_dataset(
            name = 'categories',
            robj = as.character(src[[i]][['levels']][]),
            dtype = StringType('utf8')
          )

          # Add encoding attributes to categories
          AddAnndataEncoding(dfile[[dname]][[mapped_i]][['categories']], encoding_type = 'string-array')

          # Add attributes to the categorical group
          AddAnndataEncoding(dfile[[dname]][[mapped_i]], encoding_type = 'categorical')
          dfile[[dname]][[mapped_i]]$create_attr(
            attr_name = 'ordered',
            robj = FALSE,
            dtype = GuessDType(x = FALSE),
            space = ScalarSpace()
          )
        } else {
          # Regular dataset
          dfile[[dname]]$create_dataset(
            name = mapped_i,
            robj = src[[i]][index],
            dtype = src[[i]]$get_type()
          )

          # Add encoding attributes
          dfile[[dname]][[mapped_i]]$create_attr(
            attr_name = 'encoding-type',
            robj = 'array',
            dtype = CachedGuessDType(x = 'array'),
            space = ScalarSpace()
          )
          dfile[[dname]][[mapped_i]]$create_attr(
            attr_name = 'encoding-version',
            robj = '0.2.0',
            dtype = CachedGuessDType(x = '0.2.0'),
            space = ScalarSpace()
          )
        }
      }
      # Add column order with mapped names (required by anndata >= 0.11)
      if (src$attr_exists(attr_name = 'colnames')) {
        original_colnames <- h5attr(x = src, which = 'colnames')
        if (standardize) {
          mapped_colnames <- sapply(original_colnames, function(cn) apply_name_map(cn, dname))
        } else {
          mapped_colnames <- original_colnames
        }
        dfile[[dname]]$create_attr(
          attr_name = 'column-order',
          robj = mapped_colnames,
          dtype = GuessDType(x = mapped_colnames)
        )
      } else {
        # Fallback: derive column order from group contents
        col_names <- setdiff(dfile[[dname]]$names, c('_index', '__categories'))
        dfile[[dname]]$create_attr(
          attr_name = 'column-order',
          robj = if (length(col_names) > 0) col_names else character(0),
          dtype = if (length(col_names) > 0) GuessDType(x = col_names) else StringType('utf8')
        )
      }
      encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
      names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        attr.value <- encoding.info[i]
        if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
          dfile[[dname]]$attr_delete(attr_name = attr.name)
        }
        dfile[[dname]]$create_attr(
          attr_name = attr.name,
          robj = attr.value,
          dtype = GuessDType(x = attr.value),
          space = ScalarSpace()
        )
      }
    }
    return(invisible(x = NULL))
  }
  # Because AnnData can't figure out that sparse matrices are stored as groups
  # known_ncols: pass the known number of columns to avoid reading all indices
  AddEncoding <- function(dname, known_ncols = NULL) {
    encoding.info <- c('type' = 'csr_matrix', 'version' = '0.1.0')
    names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
    if (inherits(x = dfile[[dname]], what = 'H5Group')) {
      # Add encoding type and version
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        attr.value <- encoding.info[i]
        if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
          dfile[[dname]]$attr_delete(attr_name = attr.name)
        }
        dfile[[dname]]$create_attr(
          attr_name = attr.name,
          robj = attr.value,
          dtype = CachedGuessDType(x = attr.value),
          space = ScalarSpace()
        )
      }

      # Recompute shape for AnnData convention (n_obs, n_vars).
      # Seurat stores CSC matrices (indptr = column pointers for cells),
      # which we reinterpret as CSR for h5ad (indptr = row pointers for obs).
      if (dfile[[dname]]$exists(name = 'indptr') &&
          dfile[[dname]]$exists(name = 'indices')) {
        indptr_len <- prod(dfile[[dname]][['indptr']]$dims)
        nrows <- as.integer(indptr_len - 1L)

        # Use known_ncols if available; fall back to dims attr from source
        ncols <- if (!is.null(known_ncols)) {
          as.integer(known_ncols)
        } else if (dfile[[dname]]$attr_exists(attr_name = 'dims')) {
          dims <- h5attr(dfile[[dname]], 'dims')
          as.integer(dims[1])  # dims[1] is nrow of original matrix = ncols after CSC->CSR
        } else {
          # Last resort: read indices (slow for large matrices)
          all_indices <- dfile[[dname]][['indices']][]
          if (length(all_indices) > 0) as.integer(max(all_indices) + 1L) else 0L
        }

        shape <- c(nrows, ncols)
        if (dfile[[dname]]$attr_exists(attr_name = 'shape')) {
          dfile[[dname]]$attr_delete(attr_name = 'shape')
        }
        dfile[[dname]]$create_attr(
          attr_name = 'shape',
          robj = shape,
          dtype = GuessDType(x = shape)
        )
      }
    }
    return(invisible(x = NULL))
  }
  # Add assay data
  assay.group <- source[['assays']][[assay]]

  # Check if this is a V5 h5Seurat file with layers structure
  has_layers <- assay.group$exists(name = 'layers') && inherits(assay.group[['layers']], 'H5Group')

  # Determine the appropriate data paths based on structure
  if (has_layers) {
    # V5 structure: data is under layers/
    # Prioritize layers/data (all genes) over scale.data (variable features only)
    if (assay.group[['layers']]$exists(name = 'data')) {
      x.data <- 'layers/data'
      raw.data <- if (assay.group[['layers']]$exists(name = 'counts')) {
        'layers/counts'
      } else {
        NULL
      }
    } else if (assay.group[['layers']]$exists(name = 'counts')) {
      x.data <- 'layers/counts'
      raw.data <- NULL
    } else if (assay.group$exists(name = 'scale.data')) {
      # Fallback: only scale.data available (unusual)
      x.data <- 'scale.data'
      raw.data <- NULL
    } else {
      stop("Cannot find data or counts in V5 h5Seurat file", call. = FALSE)
    }
  } else {
    # Legacy structure
    # Prioritize data (all genes) over scale.data (variable features only)
    if (source$index()[[assay]]$slots[['data']]) {
      x.data <- 'data'
      raw.data <- if (source$index()[[assay]]$slots[['counts']]) {
        'counts'
      } else {
        NULL
      }
    } else if (source$index()[[assay]]$slots[['counts']]) {
      x.data <- 'counts'
      raw.data <- NULL
    } else if (source$index()[[assay]]$slots[['scale.data']]) {
      # Fallback: only scale.data available (unusual)
      x.data <- 'scale.data'
      raw.data <- NULL
    } else {
      stop("Cannot find data or counts in h5Seurat file", call. = FALSE)
    }
  }

  if (verbose) {
    message("Adding ", x.data, " from ", assay, " as X")
  }
  # Get shape from source BEFORE copying (most reliable)
  source_dims <- NULL
  if (grepl("/", x.data)) {
    # Handle nested paths like "layers/data"
    parts <- strsplit(x.data, "/")[[1]]
    src_obj <- assay.group
    for (part in parts) {
      src_obj <- src_obj[[part]]
    }
    if (src_obj$attr_exists(attr_name = 'dims')) {
      source_dims <- h5attr(x = src_obj, which = 'dims')
    }
  } else {
    # Simple path
    if (assay.group$exists(name = x.data) &&
        assay.group[[x.data]]$attr_exists(attr_name = 'dims')) {
      source_dims <- h5attr(x = assay.group[[x.data]], which = 'dims')
    }
  }

  # Copy the data
  CopyH5MatrixData(
    src_group = assay.group,
    src_path = x.data,
    dst_loc = dfile,
    dst_name = 'X',
    verbose = verbose
  )

  # Ensure shape attribute exists (required by anndata)
  if (!dfile[['X']]$attr_exists(attr_name = 'shape')) {
    # Try to get from copied dims attribute first
    if (dfile[['X']]$attr_exists(attr_name = 'dims')) {
      dims <- h5attr(x = dfile[['X']], which = 'dims')
      dfile[['X']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = dims),
        dtype = GuessDType(x = dims)
      )
      dfile[['X']]$attr_delete(attr_name = 'dims')
    } else if (!is.null(source_dims)) {
      # Use dims from source
      dfile[['X']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = source_dims),
        dtype = GuessDType(x = source_dims)
      )
    }
  }

  AddEncoding(dname = 'X')
  # Determine n_features from the SOURCE metadata (reliable).
  # Note: reading shape from dfile[['X']] after AddEncoding is unreliable
  # due to hdf5r attribute caching after attr_delete + create_attr.
  if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
    n_features <- length(SafeH5DRead(assay.group[['meta.data/_index']]))
  } else if (assay.group$exists(name = 'features')) {
    n_features <- prod(assay.group[['features']]$dims)
  } else {
    # Last resort: close and reopen dfile to clear cache, then read shape
    dfile$close_all()
    dfile <- H5File$new(filename = dest, mode = 'r+')
    x_shape <- h5attr(x = dfile[['X']], which = 'shape')
    n_features <- x_shape[2]
  }

  x.features <- switch(
    EXPR = x.data,
    'scale.data' = {
      if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
        # V5: match meta.data/_index against scaled.features
        which(x = SafeH5DRead(assay.group[['meta.data/_index']]) %in% SafeH5DRead(assay.group[['scaled.features']]))
      } else {
        # Legacy: match features against scaled.features
        which(x = SafeH5DRead(assay.group[['features']]) %in% SafeH5DRead(assay.group[['scaled.features']]))
      }
    },
    seq.default(from = 1, to = n_features)
  )
  # Add meta.features with validation
  if (assay.group$exists(name = 'meta.features')) {
    meta.features.src <- assay.group[['meta.features']]
    if (inherits(x = meta.features.src, what = c('H5D', 'H5Group'))) {
      TransferDF(
        src = meta.features.src,
        dname = 'var',
        index = x.features
      )
    } else {
      message("meta.features is not a valid H5D or H5Group, creating empty var")
      dfile$create_group(name = 'var')
    }
  } else {
    dfile$create_group(name = 'var')
  }
  # Add feature names
  if (Exists(x = dfile[['var']], name = rownames)) {
    dfile[['var']]$link_delete(name = rownames)
  }
  # Get feature names from correct location based on structure
  if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
    # V5 structure: feature names are in meta.data/_index
    features_data <- SafeH5DRead(assay.group[['meta.data/_index']])
  } else {
    # Legacy structure: feature names are in features dataset
    features_data <- SafeH5DRead(assay.group[['features']])
  }
  # Ensure features are character strings
  features_subset <- as.character(features_data[x.features])
  dfile[['var']]$create_dataset(
    name = rownames,
    robj = features_subset,
    dtype = StringType('utf8')
  )

  # Add encoding attributes to _index
  dfile[['var']][[rownames]]$create_attr(
    attr_name = 'encoding-type',
    robj = 'string-array',
    dtype = CachedGuessDType(x = 'string-array'),
    space = ScalarSpace()
  )
  dfile[['var']][[rownames]]$create_attr(
    attr_name = 'encoding-version',
    robj = '0.2.0',
    dtype = CachedGuessDType(x = '0.2.0'),
    space = ScalarSpace()
  )

  # Add _index attribute pointing to itself
  dfile[['var']]$create_attr(
    attr_name = rownames,
    robj = rownames,
    dtype = GuessDType(x = rownames),
    space = ScalarSpace()
  )

  # Add highly variable genes information if available
  if (assay.group$exists(name = 'variable.features')) {
    if (verbose) {
      message("Adding highly variable gene information")
    }
    variable.features <- assay.group[['variable.features']][]
    all.features <- dfile[['var']][[rownames]][]

    # Create boolean array for highly variable status
    is.variable <- all.features %in% variable.features

    dfile[['var']]$create_dataset(
      name = 'highly_variable',
      robj = is.variable,
      dtype = GuessDType(x = is.variable)
    )

    # Add number of variable features to uns
    if (!dfile$exists(name = 'uns')) {
      dfile$create_group(name = 'uns')
    }
    dfile[['uns']]$create_dataset(
      name = 'n_variable_features',
      robj = sum(is.variable),
      dtype = GuessDType(x = sum(is.variable))
    )
  }
  # Because AnnData requries meta.features and can't build an empty data frame
  if (!dfile[['var']]$attr_exists(attr_name = 'column-order')) {
    var.cols <- setdiff(
      x = names(x = dfile[['var']]),
      y = c(rownames, '__categories')
    )
    if (!length(x = var.cols)) {
      var.cols <- 'features'
      dfile[['var']]$obj_copy_to(
        dst_loc = dfile[['var']],
        dst_name = var.cols,
        src_name = rownames
      )
    }
    dfile[['var']]$create_attr(
      attr_name = 'column-order',
      robj = var.cols,
      dtype = GuessDType(x = var.cols)
    )
  }
  
  # Add encoding, to ensure compatibility with python's anndata > 0.8.0:
  encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
  names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
  for (i in seq_along(along.with = encoding.info)) {
    attr.name <- names(x = encoding.info)[i]
    attr.value <- encoding.info[i]
    if (dfile[['var']]$attr_exists(attr_name = attr.name)) {
      dfile[['var']]$attr_delete(attr_name = attr.name)
    }
    dfile[['var']]$create_attr(
      attr_name = attr.name,
      robj = attr.value,
      dtype = GuessDType(x = attr.value),
      space = ScalarSpace()
    )
  }
  
  # Add raw
  if (!is.null(x = raw.data)) {
    if (verbose) {
      message("Adding ", raw.data, " from ", assay, " as raw")
    }

    # Get shape from source BEFORE copying (most reliable)
    raw_source_dims <- NULL
    if (grepl("/", raw.data)) {
      # Handle nested paths like "layers/counts"
      parts <- strsplit(raw.data, "/")[[1]]
      src_obj <- assay.group
      for (part in parts) {
        src_obj <- src_obj[[part]]
      }
      if (src_obj$attr_exists(attr_name = 'dims')) {
        raw_source_dims <- h5attr(x = src_obj, which = 'dims')
      }
    } else {
      # Simple path
      if (assay.group$exists(name = raw.data) &&
          assay.group[[raw.data]]$attr_exists(attr_name = 'dims')) {
        raw_source_dims <- h5attr(x = assay.group[[raw.data]], which = 'dims')
      }
    }

    dfile$create_group(name = 'raw')
    CopyH5MatrixData(
      src_group = assay.group,
      src_path = raw.data,
      dst_loc = dfile[['raw']],
      dst_name = 'X',
      verbose = verbose
    )

    # Ensure shape attribute exists for raw/X (required by anndata)
    if (!dfile[['raw/X']]$attr_exists(attr_name = 'shape')) {
      # Try to get from copied dims attribute first
      if (dfile[['raw/X']]$attr_exists(attr_name = 'dims')) {
        dims <- h5attr(x = dfile[['raw/X']], which = 'dims')
        dfile[['raw/X']]$create_attr(
          attr_name = 'shape',
          robj = rev(x = dims),
          dtype = GuessDType(x = dims)
        )
        dfile[['raw/X']]$attr_delete(attr_name = 'dims')
      } else if (!is.null(raw_source_dims)) {
        # Use dims from source
        dfile[['raw/X']]$create_attr(
          attr_name = 'shape',
          robj = rev(x = raw_source_dims),
          dtype = GuessDType(x = raw_source_dims)
        )
      }
    }

    AddEncoding(dname = 'raw/X')
    # Get number of features from the ACTUAL raw/X matrix dimensions
    raw_x_shape <- NULL
    if (dfile[['raw/X']]$attr_exists(attr_name = 'shape')) {
      raw_x_shape <- h5attr(x = dfile[['raw/X']], which = 'shape')
    } else if (dfile[['raw/X']]$attr_exists(attr_name = 'dims')) {
      raw_x_shape <- rev(h5attr(x = dfile[['raw/X']], which = 'dims'))
    }
    if (!is.null(raw_x_shape) && length(raw_x_shape) >= 2) {
      n_raw_features <- raw_x_shape[2]
    } else {
      n_raw_features <- prod(assay.group[['features']]$dims)
    }
    # Add meta.features with validation
    if (assay.group$exists(name = 'meta.features')) {
      meta.features.src <- assay.group[['meta.features']]
      if (inherits(x = meta.features.src, what = c('H5D', 'H5Group'))) {
        TransferDF(
          src = meta.features.src,
          dname = 'raw/var',
          index = seq.default(from = 1, to = n_raw_features)
        )
      } else {
        message("meta.features is not a valid H5D or H5Group in raw, creating empty var")
        dfile[['raw']]$create_group(name = 'var')
      }
    } else {
      dfile[['raw']]$create_group(name = 'var')
    }
    # Add feature names
    if (Exists(x = dfile[['raw/var']], name = rownames)) {
      dfile[['raw/var']]$link_delete(name = rownames)
    }
    # Get feature names from correct location based on structure
    if (has_layers && assay.group$exists(name = 'meta.data/_index')) {
      # V5 structure: feature names are in meta.data/_index
      features_data_raw <- SafeH5DRead(assay.group[['meta.data/_index']])
    } else {
      # Legacy structure: feature names are in features dataset
      features_data_raw <- SafeH5DRead(assay.group[['features']])
    }
    # Subset to match actual raw/X dimensions
    features_data_raw <- as.character(features_data_raw[seq_len(n_raw_features)])
    dfile[['raw/var']]$create_dataset(
      name = rownames,
      robj = features_data_raw,
      dtype = GuessDType(x = features_data_raw[1])
    )
    # Add encoding attributes to raw/var/_index (required for AnnData to read gene names correctly)
    dfile[['raw/var']][[rownames]]$create_attr(
      attr_name = 'encoding-type',
      robj = 'string-array',
      dtype = CachedGuessDType(x = 'string-array'),
      space = ScalarSpace()
    )
    dfile[['raw/var']][[rownames]]$create_attr(
      attr_name = 'encoding-version',
      robj = '0.2.0',
      dtype = CachedGuessDType(x = '0.2.0'),
      space = ScalarSpace()
    )
    dfile[['raw/var']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames),
      space = ScalarSpace()
    )
    # Add dataframe encoding attributes to raw/var group (required for AnnData)
    # Only add column-order if it doesn't already exist (TransferDF may have added it)
    if (!dfile[['raw/var']]$attr_exists(attr_name = 'column-order')) {
      # Get list of columns (everything except _index)
      raw_var_cols <- setdiff(dfile[['raw/var']]$names, rownames)
      if (length(raw_var_cols) > 0) {
        dfile[['raw/var']]$create_attr(
          attr_name = 'column-order',
          robj = raw_var_cols,
          dtype = GuessDType(x = raw_var_cols)
        )
      } else {
        # Empty column-order for dataframes with only _index
        dfile[['raw/var']]$create_attr(
          attr_name = 'column-order',
          robj = character(0),
          dtype = StringType('utf8')
        )
      }
    }
    # Add encoding-type and encoding-version
    encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
    names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
    for (i in seq_along(along.with = encoding.info)) {
      attr.name <- names(x = encoding.info)[i]
      attr.value <- encoding.info[i]
      if (dfile[['raw/var']]$attr_exists(attr_name = attr.name)) {
        dfile[['raw/var']]$attr_delete(attr_name = attr.name)
      }
      dfile[['raw/var']]$create_attr(
        attr_name = attr.name,
        robj = attr.value,
        dtype = GuessDType(x = attr.value),
        space = ScalarSpace()
      )
    }
  }
  # Add cell-level metadata with validation
  if (source$exists(name = 'meta.data')) {
    meta.src <- source[['meta.data']]

    # First attempt normal transfer
    TransferDF(
      src = meta.src,
      dname = 'obs',
      index = seq.default(from = 1, to = length(x = Cells(x = source)))
    )

    # If obs was not created successfully, try direct extraction
    if (!dfile$exists(name = 'obs')) {
      if (verbose) {
        message("Attempting direct metadata extraction for V5 compatibility")
      }

      # Create obs group
      dfile$create_group(name = 'obs')

      # Try to read metadata directly from h5Seurat structure
      if (inherits(x = meta.src, what = 'H5Group')) {
        # Get all column names except internal ones
        meta_cols <- setdiff(names(x = meta.src), c('__categories', '_index'))

        for (col in meta_cols) {
          tryCatch({
            if (IsFactor(x = meta.src[[col]])) {
              # Use newer anndata format: each categorical is a group with categories/codes
              dfile[['obs']]$create_group(name = col)

              # Convert codes to 0-based indices, handling NA values
              # R factor codes are 1-based; h5ad uses 0-based with -1 for missing
              raw_codes <- meta.src[[col]][['values']][]
              codes_values <- ifelse(is.na(raw_codes), -1L, raw_codes - 1L)
              n_categories <- length(meta.src[[col]][['levels']][])
              has_na <- any(is.na(raw_codes))

              # Choose appropriate integer dtype based on number of categories
              # Must use signed integers when NA values present (need -1)
              codes_dtype <- if (has_na || n_categories <= 127) {
                hdf5r::h5types$H5T_NATIVE_INT8
              } else if (n_categories <= 255) {
                hdf5r::h5types$H5T_NATIVE_UINT8
              } else if (n_categories <= 32767) {
                hdf5r::h5types$H5T_NATIVE_INT16
              } else if (n_categories <= 65535) {
                hdf5r::h5types$H5T_NATIVE_UINT16
              } else {
                hdf5r::h5types$H5T_NATIVE_INT32
              }

              # Add codes dataset with proper dtype for scanpy compatibility
              dfile[['obs']][[col]]$create_dataset(
                name = 'codes',
                robj = codes_values,
                dtype = codes_dtype
              )

              # Add encoding attributes to codes
              dfile[['obs']][[col]][['codes']]$create_attr(
                attr_name = 'encoding-type',
                robj = 'array',
                dtype = CachedGuessDType(x = 'array'),
                space = ScalarSpace()
              )
              dfile[['obs']][[col]][['codes']]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = CachedGuessDType(x = '0.2.0'),
                space = ScalarSpace()
              )

              # Add categories dataset
              dfile[['obs']][[col]]$create_dataset(
                name = 'categories',
                robj = as.character(meta.src[[col]][['levels']][]),
                dtype = StringType('utf8')
              )

              # Add encoding attributes to categories
              dfile[['obs']][[col]][['categories']]$create_attr(
                attr_name = 'encoding-type',
                robj = 'string-array',
                dtype = CachedGuessDType(x = 'string-array'),
                space = ScalarSpace()
              )
              dfile[['obs']][[col]][['categories']]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = CachedGuessDType(x = '0.2.0'),
                space = ScalarSpace()
              )

              # Add attributes to the categorical group
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-type',
                robj = 'categorical',
                dtype = CachedGuessDType(x = 'categorical'),
                space = ScalarSpace()
              )
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = CachedGuessDType(x = '0.2.0'),
                space = ScalarSpace()
              )
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'ordered',
                robj = FALSE,
                dtype = GuessDType(x = FALSE),
                space = ScalarSpace()
              )
            } else {
              # Handle regular columns
              dfile[['obs']]$create_dataset(
                name = col,
                robj = meta.src[[col]][],
                dtype = meta.src[[col]]$get_type()
              )

              # Add encoding attributes
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-type',
                robj = 'array',
                dtype = CachedGuessDType(x = 'array'),
                space = ScalarSpace()
              )
              dfile[['obs']][[col]]$create_attr(
                attr_name = 'encoding-version',
                robj = '0.2.0',
                dtype = CachedGuessDType(x = '0.2.0'),
                space = ScalarSpace()
              )
            }
          }, error = function(e) {
            if (verbose) {
              message("Could not transfer metadata column ", col, ": ", e$message)
            }
          })
        }

        # Add column-order attribute (required by anndata >= 0.11)
        if (!dfile[['obs']]$attr_exists(attr_name = 'column-order')) {
          if (length(meta_cols) > 0) {
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = meta_cols,
              dtype = GuessDType(x = meta_cols)
            )
          } else {
            dfile[['obs']]$create_attr(
              attr_name = 'column-order',
              robj = character(0),
              dtype = StringType('utf8')
            )
          }
        }

        # Add encoding attributes
        encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
        names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
        for (i in seq_along(along.with = encoding.info)) {
          attr.name <- names(x = encoding.info)[i]
          if (!dfile[['obs']]$attr_exists(attr_name = attr.name)) {
            dfile[['obs']]$create_attr(
              attr_name = attr.name,
              robj = encoding.info[i],
              dtype = GuessDType(x = encoding.info[i]),
              space = ScalarSpace()
            )
          }
        }
      } else {
        message("meta.data structure not recognized, creating minimal obs")
      }
    }
  } else {
    # Create empty obs group if meta.data doesn't exist
    dfile$create_group(name = 'obs')
  }

  # Add cell names if obs was created
  if (dfile$exists(name = 'obs')) {
    if (Exists(x = dfile[['obs']], name = rownames)) {
      dfile[['obs']]$link_delete(name = rownames)
    }
    cell_names <- as.character(Cells(x = source))
    dfile[['obs']]$create_dataset(
      name = rownames,
      robj = cell_names,
      dtype = StringType('utf8')
    )

    # Add encoding attributes to _index
    dfile[['obs']][[rownames]]$create_attr(
      attr_name = 'encoding-type',
      robj = 'string-array',
      dtype = CachedGuessDType(x = 'string-array'),
      space = ScalarSpace()
    )
    dfile[['obs']][[rownames]]$create_attr(
      attr_name = 'encoding-version',
      robj = '0.2.0',
      dtype = CachedGuessDType(x = '0.2.0'),
      space = ScalarSpace()
    )

    # Add _index attribute pointing to itself
    dfile[['obs']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames),
      space = ScalarSpace()
    )

    # Ensure encoding attributes exist for anndata compatibility
    if (!dfile[['obs']]$attr_exists(attr_name = 'encoding-type')) {
      encoding.info <- c('type' = 'dataframe', 'version' = '0.2.0')
      names(x = encoding.info) <- paste0('encoding-', names(x = encoding.info))
      for (i in seq_along(along.with = encoding.info)) {
        attr.name <- names(x = encoding.info)[i]
        if (!dfile[['obs']]$attr_exists(attr_name = attr.name)) {
          dfile[['obs']]$create_attr(
            attr_name = attr.name,
            robj = encoding.info[i],
            dtype = GuessDType(x = encoding.info[i]),
            space = ScalarSpace()
          )
        }
      }
    }
  }
  # Add dimensional reduction information
  if (!dfile$exists(name = 'obsm')) {
    obsm <- dfile$create_group(name = 'obsm')
  } else {
    obsm <- dfile[['obsm']]
  }
  if (!dfile$exists(name = 'varm')) {
    varm <- dfile$create_group(name = 'varm')
  } else {
    varm <- dfile[['varm']]
  }
  reductions <- source$index()[[assay]]$reductions
  for (reduc in names(x = reductions)) {
    if (verbose) {
      message("Adding dimensional reduction information for ", reduc)
    }
    # hdf5r reads R matrix as (n_cells, n_dims), but when written back,
    # HDF5/h5py sees reversed dims. scTranspose so anndata sees (n_cells, n_dims).
    emb_src <- source[[H5Path('reductions', reduc, 'cell.embeddings')]]
    emb_data <- emb_src$read()
    if (is.matrix(emb_data) && nrow(emb_data) > ncol(emb_data)) {
      emb_data <- t(emb_data)
    }
    obsm$create_dataset(
      name = paste0('X_', reduc),
      robj = emb_data,
      dtype = emb_src$get_type()
    )
    if (reductions[[reduc]]['feature.loadings']) {
      if (verbose) {
        message("Adding feature loadings for ", reduc)
      }
      loadings <- source[['reductions']][[reduc]][['feature.loadings']]
      reduc.features <- loadings$dims[1]
      x.features <- dfile[['var']][[rownames]]$dims
      varm.name <- switch(EXPR = reduc, 'pca' = 'PCs', toupper(x = reduc))
      # Because apparently AnnData requires nPCs == nrow(X)
      if (reduc.features < x.features) {
        pad <- paste0('pad_', varm.name)
        # Get feature indices and filter out NAs
        feature_indices <- match(
          x = SafeH5DRead(source[['reductions']][[reduc]][['features']]),
          table = SafeH5DRead(dfile[['var']][[rownames]])
        )
        # Replace NA values with new rows (append to end)
        na_indices <- which(is.na(feature_indices))
        if (length(na_indices) > 0) {
          # For NA indices, assign new row positions at the end
          feature_indices[na_indices] <- (reduc.features + 1L):(reduc.features + length(na_indices))
        }

        PadMatrix(
          src = loadings,
          dest = dfile[['varm']],
          dname = pad,
          dims = c(x.features, loadings$dims[2]),
          index = list(
            feature_indices,
            seq.default(from = 1, to = loadings$dims[2])
          )
        )
        loadings <- dfile[['varm']][[pad]]
      }
      scTranspose(x = loadings, dest = varm, dname = varm.name, verbose = FALSE)
      if (reduc.features < x.features) {
        dfile$link_delete(name = loadings$get_obj_name())
      }
    }
  }
  # Add global dimensional reduction information
  global.reduc <- source$index()[['global']][['reductions']]
  for (reduc in global.reduc) {
    if (reduc %in% names(x = reductions)) {
      next
    } else if (verbose) {
      message("Adding dimensional reduction information for ", reduc, " (global)")
    }
    emb_src <- source[[H5Path('reductions', reduc, 'cell.embeddings')]]
    emb_data <- emb_src$read()
    if (is.matrix(emb_data) && nrow(emb_data) > ncol(emb_data)) {
      emb_data <- t(emb_data)
    }
    obsm$create_dataset(
      name = paste0('X_', reduc),
      robj = emb_data,
      dtype = emb_src$get_type()
    )
  }
  # Add spatial coordinates if available
  # Check for spatial data in images slot
  if (source$exists(name = 'images')) {
    cell_names <- SafeH5DRead(source[['cell.names']])
    images <- names(x = source[['images']])
    if (length(images) > 0 && verbose) {
      message("Processing spatial data from images")
    }

    # Process first image for now (extend for multiple images later)
    if (length(images) > 0) {
      img_name <- images[1]
      img_group <- source[['images']][[img_name]]

      # Check for coordinates
      spatial_matrix <- NULL
      if (img_group$exists(name = 'coordinates')) {
        if (verbose) {
          message("Adding spatial coordinates to obsm['spatial']")
        }

        # Read coordinates
        coord_group <- img_group[['coordinates']]

        # Get x and y coordinates
        if (coord_group$exists('x') && coord_group$exists('y')) {
          x_coords <- coord_group[['x']][]
          y_coords <- coord_group[['y']][]

          # Create spatial matrix (cells x 2) as [X, Y] for scanpy/squidpy
          spatial_matrix <- cbind(x_coords, y_coords)
        } else if (coord_group$exists('imagerow') && coord_group$exists('imagecol')) {
          # Visium-style coordinates: imagerow=Y, imagecol=X
          row_coords <- coord_group[['imagerow']][]
          col_coords <- coord_group[['imagecol']][]

          # Create spatial matrix (cells x 2) as [X, Y] = [col, row]
          spatial_matrix <- cbind(col_coords, row_coords)
        }
      } else if (img_group$exists(name = 'boundaries') &&
                 img_group[['boundaries']]$exists(name = 'centroids')) {
        centroid_group <- img_group[['boundaries/centroids']]
        if (centroid_group$exists(name = 'coords')) {
          coords_mat <- centroid_group[['coords']]$read()
          row_index <- seq_len(nrow(coords_mat))
          if (centroid_group$exists(name = 'cells')) {
            centroid_cells <- centroid_group[['cells']][]
            if (!is.null(cell_names)) {
              matched <- match(cell_names, centroid_cells)
              if (any(is.na(matched))) {
                message("Unable to align all centroid coordinates with cell names")
              } else {
                row_index <- matched
              }
            }
          }
          coords_mat <- coords_mat[row_index, , drop = FALSE]
          # coords_mat is stored as [X, Y] in Seurat h5Seurat, keep as-is for scanpy
          spatial_matrix <- coords_mat
        }
      }
      if (!is.null(spatial_matrix)) {
        spatial_matrix <- as.matrix(spatial_matrix)
        storage.mode(spatial_matrix) <- 'double'
        obsm$create_dataset(
          name = 'spatial',
          robj = spatial_matrix,
          dtype = h5types$H5T_NATIVE_DOUBLE
        )
        scTranspose(
          x = obsm[['spatial']],
          dest = obsm,
          dname = 'spatial',
          overwrite = TRUE,
          verbose = FALSE
        )
        obsm[['spatial']]$create_attr(
          attr_name = 'encoding-type',
          robj = 'array',
          dtype = CachedGuessDType(x = 'array'),
          space = ScalarSpace()
        )
        obsm[['spatial']]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.2.0',
          dtype = CachedGuessDType(x = '0.2.0'),
          space = ScalarSpace()
        )
      }
    }
  }

  # Create uns if it doesn't exist
  if (!dfile$exists(name = 'uns')) {
    dfile$create_group(name = 'uns')
  }

  # Preserve Seurat-specific metadata in uns['seurat']
  if (verbose) {
    message("Preserving Seurat-specific metadata in uns['seurat']")
  }
  seurat_group <- dfile[['uns']]$create_group(name = 'seurat')

  # Store assay name
  seurat_group$create_dataset(
    name = 'assay_name',
    robj = assay,
    dtype = StringType('utf8')
  )

  # Store h5Seurat format version
  seurat_group$create_dataset(
    name = 'version',
    robj = 'SeuratDisk-V5',
    dtype = StringType('utf8')
  )

  # Preserve graph names if they exist
  if (source$exists(name = 'graphs')) {
    graph_names <- names(x = source[['graphs']])
    if (length(graph_names) > 0) {
      seurat_group$create_dataset(
        name = 'graph_names',
        robj = graph_names,
        dtype = StringType('utf8')
      )
    }
  }

  # Preserve obs column names to track which metadata came from Seurat
  if (dfile$exists(name = 'obs')) {
    if (dfile[['obs']]$attr_exists(attr_name = 'column-order')) {
      obs_cols <- h5attr(x = dfile[['obs']], which = 'column-order')
      seurat_group$create_dataset(
        name = 'obs_columns',
        robj = obs_cols,
        dtype = StringType('utf8')
      )
    }
  }

  # Add spatial metadata to uns if available
  if (source$exists(name = 'images')) {
    images <- names(x = source[['images']])

    if (length(images) > 0) {
      if (verbose) {
        message("Adding spatial metadata to uns['spatial']")
      }

      # Create spatial group in uns
      spatial_uns <- dfile[['uns']]$create_group(name = 'spatial')

      # Each Seurat image becomes a separate library in h5ad format
      # This matches the multi-library structure expected by scanpy/squidpy
      for (img_name in images) {
        img_group <- source[['images']][[img_name]]

        # Use the Seurat image name as the library ID (e.g., 'anterior1')
        lib_id <- img_name
        lib_group <- spatial_uns$create_group(name = lib_id)

        # Add scale factors for this library
        if (img_group$exists(name = 'scale.factors')) {
          sf_group <- lib_group$create_group(name = 'scalefactors')
          sf_src <- img_group[['scale.factors']]
          sf_map <- c(
            hires = 'tissue_hires_scalef',
            lowres = 'tissue_lowres_scalef',
            spot = 'spot_diameter_fullres',
            fiducial = 'fiducial_diameter_fullres'
          )
          for (sf_name in sf_src$names) {
            dst_name <- sf_map[[sf_name]]
            if (is.null(dst_name)) {
              dst_name <- sf_name
            }
            sf_value <- sf_src[[sf_name]][]
            sf_group$create_dataset(
              name = dst_name,
              robj = sf_value,
              dtype = h5types$H5T_NATIVE_DOUBLE,
              space = ScalarSpace(),
              chunk_dims = NULL
            )
          }
        }

        # Add metadata
        meta_group <- lib_group$create_group(name = 'metadata')

        # Create images group with standard h5ad keys ('hires', 'lowres')
        images_group <- lib_group$create_group(name = 'images')

        # Add primary image as 'lowres' and 'hires'
        if (img_group$exists(name = 'image')) {
          img_data <- img_group[['image']]$read()
          if (length(dim(img_data)) == 3L) {
            img_arr <- aperm(img_data, c(3L, 2L, 1L))

            # Check for high-resolution image in attributes
            has_hires <- FALSE
            if (!is.null(img_group[['image']]$attr_names) &&
                'hires.image' %in% img_group[['image']]$attr_names) {
              hires_data <- h5attr(x = img_group[['image']], which = 'hires.image')
              if (!is.null(hires_data) && length(dim(hires_data)) == 3L) {
                images_group$create_dataset(
                  name = 'hires',
                  robj = aperm(hires_data, c(3L, 2L, 1L)),
                  dtype = h5types$H5T_NATIVE_DOUBLE
                )
                has_hires <- TRUE
              }
            }

            # If no separate hires image, use the primary image for both hires and lowres
            # This ensures compatibility with squidpy which defaults to 'hires'
            if (!has_hires) {
              images_group$create_dataset(
                name = 'hires',
                robj = img_arr,
                dtype = h5types$H5T_NATIVE_DOUBLE
              )
            }

            # Always create lowres
            images_group$create_dataset(
              name = 'lowres',
              robj = img_arr,
              dtype = h5types$H5T_NATIVE_DOUBLE
            )
          }
        }
      }
    }
  }

  # Add graphs to obsp (pairwise observations)
  # Include graphs from both the assay index and no.assay (untagged graphs)
  graphs.available <- unique(c(
    source$index()[[assay]]$graphs,
    source$index()[['no.assay']]$graphs
  ))
  if (length(graphs.available) > 0 && source$exists('graphs')) {
    if (verbose) {
      message("Adding graph information to obsp")
    }
    if (!dfile$exists(name = 'obsp')) {
      dfile$create_group(name = 'obsp')
    }

    # Transfer all available graphs
    for (graph.name in graphs.available) {
      if (source[['graphs']]$exists(name = graph.name)) {
        # Map Seurat graph names to anndata conventions
        # RNA_nn (kNN distances) -> distances, RNA_snn (shared NN weights) -> connectivities
        obsp.name <- if (grepl("_nn$", graph.name)) {
          "distances"
        } else if (grepl("_snn$", graph.name)) {
          "connectivities"
        } else {
          gsub(paste0("^", assay, "_"), "", graph.name)
        }

        if (verbose) {
          message("  - Adding ", graph.name, " as obsp/", obsp.name)
        }

        source[['graphs']]$obj_copy_to(
          dst_loc = dfile[['obsp']],
          dst_name = obsp.name,
          src_name = graph.name
        )

        # Add shape attribute for anndata compatibility
        if (source[['graphs']][[graph.name]]$attr_exists(attr_name = 'dims')) {
          dims <- h5attr(x = source[['graphs']][[graph.name]], which = 'dims')
          dfile[['obsp']][[obsp.name]]$create_attr(
            attr_name = 'shape',
            robj = rev(x = dims),
            dtype = GuessDType(x = dims)
          )
        }

        # Add encoding type
        AddEncoding(dname = paste0('obsp/', obsp.name))
      }
    }
  }

  # Add graph (legacy - keep for backward compatibility)
  graph <- source$index()[[assay]]$graphs
  graph <- graph[length(x = graph)]
  if (!is.null(x = graph)) {
    if (verbose) {
      message("Adding ", graph, " as neighbors")
    }
    dgraph <- dfile[['uns']]$create_group(name = 'neighbors')
    source[['graphs']]$obj_copy_to(
      dst_loc = dgraph,
      dst_name = 'distances',
      src_name = graph
    )
    if (source[['graphs']][[graph]]$attr_exists(attr_name = 'dims')) {
      dims <- h5attr(x = source[['graphs']][[graph]], which = 'dims')
      dgraph[['distances']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = dims),
        dtype = GuessDType(x = dims)
      )
    }
    AddEncoding(dname = 'uns/neighbors/distances')
    # Add parameters
    dgraph$create_group(name = 'params')
    dgraph[['params']]$create_dataset(
      name = 'method',
      robj = gsub(pattern = paste0('^', assay, '_'), replacement = '', x = graph),
      dtype = GuessDType(x = graph)
    )
    cmdlog <- paste(
      paste0('FindNeighbors.', assay),
      unique(x = c(names(x = reductions), source$index()$global$reductions)),
      sep = '.',
      collapse = '|'
    )
    cmdlog <- grep(
      pattern = cmdlog,
      x = names(x = source[['commands']]),
      value = TRUE
    )
    if (length(x = cmdlog) > 1) {
      timestamps <- sapply(
        X = cmdlog,
        FUN = function(cmd) {
          ts <- if (source[['commands']][[cmd]]$attr_exists(attr_name = 'time.stamp')) {
            h5attr(x = source[['commands']][[cmd]], which = 'time.stamp')
          } else {
            NULL
          }
          return(ts)
        },
        simplify = TRUE,
        USE.NAMES = FALSE
      )
      timestamps <- Filter(f = Negate(f = is.null), x = timestamps)
      cmdlog <- cmdlog[order(timestamps, decreasing = TRUE)][1]
    }
    if (length(x = cmdlog) && !is.na(x = cmdlog)) {
      cmdlog <- source[['commands']][[cmdlog]]
      if ('k.param' %in% names(x = cmdlog)) {
        dgraph[['params']]$obj_copy_from(
          src_loc = cmdlog,
          src_name = 'k.param',
          dst_name = 'n_neighbors'
        )
      }
    }
  }
  # Add layers
  other.assays <- setdiff(
    x = names(x = source$index()),
    y = c(assay, 'global', 'no.assay')
  )
  if (length(x = other.assays)) {
    x.dims <- Dims(x = dfile[['X']])
    layers <- dfile$create_group(name = 'layers')
    for (other in other.assays) {
      layer.slot <- NULL
      other.group <- source[['assays']][[other]]

      # Check if this assay has V5 structure with layers
      other.has_layers <- other.group$exists(name = 'layers') && inherits(other.group[['layers']], 'H5Group')

      for (slot in c('scale.data', 'data')) {
        # Determine the actual path for the slot
        actual_path <- if (other.has_layers && slot %in% c('data', 'counts')) {
          paste0('layers/', slot)
        } else {
          slot
        }

        # Check if the slot exists at the determined path
        slot.exists <- if (other.has_layers && slot %in% c('data', 'counts')) {
          other.group[['layers']]$exists(name = slot)
        } else {
          other.group$exists(name = slot)
        }

        if (slot.exists) {
          slot.dims <- Dims(x = other.group[[actual_path]])
          if (isTRUE(all.equal(slot.dims, x.dims))) {
            layer.slot <- actual_path
            break
          }
        }
      }

      if (!is.null(x = layer.slot)) {
        if (verbose) {
          message("Adding ", layer.slot, " from ", other, " as a layer")
        }
        layers$obj_copy_from(
          src_loc = source[['assays']][[other]],
          src_name = layer.slot,
          dst_name = other
        )
        if (layers[[other]]$attr_exists(attr_name = 'dims')) {
          dims <- h5attr(x = layers[[other]], which = 'dims')
          layers[[other]]$create_attr(
            attr_name = 'shape',
            robj = rev(x = dims),
            dtype = GuessDType(x = dims)
          )
          layers[[other]]$attr_delete(attr_name = 'dims')
        }
        AddEncoding(dname = paste('layers', other, sep = '/'))
        layer.features <- switch(
          EXPR = layer.slot,
          'scale.data' = 'scaled.features',
          'features'
        )
        var.name <- paste0(other, '_features')
        dfile[['var']]$obj_copy_from(
          src_loc = source[['assays']][[other]],
          src_name = layer.features,
          dst_name = var.name
        )
        col.order <- h5attr(x = dfile[['var']], which = 'column-order')
        col.order <- c(col.order, var.name)
        dfile[['var']]$attr_rename(
          old_attr_name = 'column-order',
          new_attr_name = 'old-column-order'
        )
        dfile[['var']]$create_attr(
          attr_name = 'column-order',
          robj = col.order,
          dtype = GuessDType(x = col.order)
        )
        dfile[['var']]$attr_delete(attr_name = 'old-column-order')
      }
    }
  }

  # Ensure all required anndata groups exist (even if empty)
  required_groups <- c('obsm', 'obsp', 'varm', 'varp', 'layers', 'uns')
  for (group_name in required_groups) {
    if (!dfile$exists(name = group_name)) {
      dfile$create_group(name = group_name)
    }
  }

  # Add encoding attributes to all top-level groups
  groups_to_encode <- c('obsm', 'obsp', 'varm', 'varp', 'layers', 'uns')
  for (group_name in groups_to_encode) {
    if (dfile$exists(name = group_name) && inherits(dfile[[group_name]], 'H5Group')) {
      # Add encoding-type
      if (!dfile[[group_name]]$attr_exists(attr_name = 'encoding-type')) {
        dfile[[group_name]]$create_attr(
          attr_name = 'encoding-type',
          robj = 'dict',
          dtype = CachedGuessDType(x = 'dict'),
          space = ScalarSpace()
        )
      }
      # Add encoding-version
      if (!dfile[[group_name]]$attr_exists(attr_name = 'encoding-version')) {
        dfile[[group_name]]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.1.0',
          dtype = CachedGuessDType(x = '0.1.0'),
          space = ScalarSpace()
        )
      }
    }
  }

  # Add encoding attributes to dimensional reductions in obsm
  if (dfile$exists(name = 'obsm')) {
    for (reduc_name in names(dfile[['obsm']])) {
      if (!dfile[['obsm']][[reduc_name]]$attr_exists(attr_name = 'encoding-type')) {
        dfile[['obsm']][[reduc_name]]$create_attr(
          attr_name = 'encoding-type',
          robj = 'array',
          dtype = CachedGuessDType(x = 'array'),
          space = ScalarSpace()
        )
      }
      if (!dfile[['obsm']][[reduc_name]]$attr_exists(attr_name = 'encoding-version')) {
        dfile[['obsm']][[reduc_name]]$create_attr(
          attr_name = 'encoding-version',
          robj = '0.2.0',
          dtype = CachedGuessDType(x = '0.2.0'),
          space = ScalarSpace()
        )
      }
    }
  }

  dfile$flush()
  return(dfile)
}



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
#' Exploits the zero-copy CSC↔CSR reinterpretation: dgCMatrix(genes×cells) CSC
#' is identical to h5ad's cells×genes CSR with relabeled arrays.
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
  mat <- NULL
  for (ln in c("counts", "data")) {
    tryCatch({
      m <- GetAssayData(assay_obj, layer = ln)
      if (inherits(m, "dgCMatrix") && length(m@x) > 0) { mat <- m; break }
    }, error = function(e) NULL)
  }
  if (is.null(mat)) return(FALSE)

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

  result <- .Call("C_write_h5ad",
    filename, mat_list, meta, reductions, graphs,
    assay_name, as.integer(gzip_level),
    PACKAGE = "scConvert"
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
  data_mat <- tryCatch({
    d <- GetAssayData(object, assay = assay, layer = 'data')
    if (!is.null(d) && prod(dim(d)) > 0) d else NULL
  }, error = function(e) NULL)

  counts_mat <- tryCatch({
    d <- GetAssayData(object, assay = assay, layer = 'counts')
    if (!is.null(d) && prod(dim(d)) > 0) d else NULL
  }, error = function(e) NULL)

  if (!is.null(data_mat)) {
    write_csr_group(dfile, "X", data_mat)
    if (!is.null(counts_mat)) {
      raw_grp <- dfile$create_group("raw")
      write_csr_group(raw_grp, "X", counts_mat)
    }
  } else if (!is.null(counts_mat)) {
    write_csr_group(dfile, "X", counts_mat)
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
  if (!is.null(data_mat) && !is.null(counts_mat)) {
    raw_gene_names <- rownames(counts_mat)
    raw_var_df <- data.frame(row.names = raw_gene_names)
    write_df_group(dfile[["raw"]], "var", raw_var_df, raw_gene_names)
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

    # Write misc items (exclude __varp__ which is written to the varp group)
    for (item_name in names(misc)) {
      if (item_name == "__varp__") next
      item <- misc[[item_name]]
      if (is.null(item)) next
      if (is.character(item) && length(item) == 1) {
        uns_grp$create_dataset(item_name, robj = item, dtype = CachedUtf8Type())
      } else if (is.numeric(item) && length(item) == 1) {
        uns_grp$create_dataset(item_name, robj = item)
      } else if (is.data.frame(item)) {
        write_df_group(uns_grp, item_name, item, rownames(item) %||% as.character(seq_len(nrow(item))))
      }
    }

    # Delegate spatial data
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# H5MU Conversion Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert H5MU files to h5Seurat files
#'
#' @inheritParams scConvert
#'
#' @return Returns a handle to \code{dest} as an \code{\link{h5Seurat}} object
#'
#' @keywords internal
#'
H5MUToH5Seurat <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Converting H5MU to h5Seurat via Seurat object...")
  }

  # Load h5mu file as Seurat object
  seurat_obj <- readH5MU(
    file = source$filename,
    verbose = verbose
  )

  # Save as h5Seurat
  h5seurat_file <- writeH5Seurat(
    object = seurat_obj,
    filename = dest,
    overwrite = overwrite,
    verbose = verbose
  )

  # Return h5Seurat connection
  dfile <- h5Seurat$new(filename = dest, mode = 'r')
  return(dfile)
}


#' Convert h5Seurat files to H5MU files
#'
#' @inheritParams scConvert
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}} object
#'
#' @keywords internal
#'
H5SeuratToH5MU <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Converting h5Seurat to H5MU via Seurat object...")
  }

  # Load h5Seurat as Seurat object
  seurat_obj <- readH5Seurat(
    file = source$filename,
    verbose = verbose
  )

  # Save as h5mu
  h5mu_file <- writeH5MU(
    object = seurat_obj,
    filename = dest,
    overwrite = overwrite,
    verbose = verbose
  )

  # Return H5File connection to h5mu
  dfile <- H5File$new(filename = dest, mode = 'r')
  return(dfile)
}


#' Convert H5MU files to H5AD files (extract single modality)
#'
#' @inheritParams scConvert
#' @param modality Name of modality to extract from h5mu file
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}} object
#'
#' @keywords internal
#'
H5MUToH5AD <- function(
  source,
  dest,
  modality = 'rna',
  overwrite = FALSE,
  verbose = TRUE
) {
  if (verbose) {
    message("Extracting modality '", modality, "' from H5MU to H5AD...")
  }

  # Load h5mu file
  seurat_obj <- readH5MU(
    file = source$filename,
    modalities = modality,
    verbose = verbose
  )

  # Get the corresponding assay name
  assay_names <- Assays(seurat_obj)
  if (length(assay_names) == 0) {
    stop("No assays found in converted object", call. = FALSE)
  }

  # Use first assay (should be the only one if modality was specified)
  target_assay <- assay_names[1]

  # Save as h5Seurat first
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit(file.remove(temp_h5seurat), add = TRUE)

  writeH5Seurat(
    object = seurat_obj,
    filename = temp_h5seurat,
    overwrite = TRUE,
    verbose = FALSE
  )

  # Convert h5Seurat to h5ad
  temp_h5seurat_conn <- scConnect(filename = temp_h5seurat, force = TRUE)
  dfile <- H5SeuratToH5AD(
    source = temp_h5seurat_conn,
    dest = dest,
    assay = target_assay,
    overwrite = overwrite,
    verbose = verbose
  )

  temp_h5seurat_conn$close_all()

  return(dfile)
}


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
        # Strip assay prefix if present (e.g. "RNA_snn" → "snn")
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Direct format converters: zarr
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an h5ad file to zarr format
#'
#' Converts an AnnData h5ad file to zarr format. By default uses streaming
#' conversion that copies fields directly without loading into R memory.
#' Both formats follow the same AnnData on-disk specification.
#'
#' @param source Path to input .h5ad file
#' @param dest Path for output .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly without creating
#'   a Seurat intermediate. If FALSE, load as Seurat then save as zarr.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5ADToZarr <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_h5ad_to_zarr(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5AD(source, verbose = verbose)
    writeZarr(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_h5ad_to_zarr <- function(h5ad_path, zarr_path, gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for H5ADToZarr", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for H5ADToZarr", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
  on.exit(h5$close_all())

  compressor <- list(id = "zlib", level = as.integer(gzip))

  # Root group
  .zarr_create_group(zarr_path, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # Helper: stream a DataFrame group (obs or var)
  .stream_h5_df_to_zarr <- function(h5_group, zarr_group_path) {
    if (!inherits(h5_group, "H5Group")) return()

    # Read index
    index <- NULL
    index_name <- "_index"
    if (h5_group$attr_exists("_index")) {
      index_name <- hdf5r::h5attr(h5_group, "_index")
    }
    if (h5_group$exists(index_name)) {
      index <- as.character(h5_group[[index_name]][])
    } else if (h5_group$exists("_index")) {
      index <- as.character(h5_group[["_index"]][])
    } else if (h5_group$exists("index")) {
      index <- as.character(h5_group[["index"]][])
    }
    if (is.null(index)) return()

    # Get column-order
    col_order <- if (h5_group$attr_exists("column-order")) {
      hdf5r::h5attr(h5_group, "column-order")
    } else {
      setdiff(names(h5_group), c("_index", "index", "__categories"))
    }

    # Create group with DataFrame attrs
    .zarr_create_group(zarr_group_path, attrs = list(
      `encoding-type` = "dataframe",
      `encoding-version` = "0.2.0",
      `_index` = "_index",
      `column-order` = col_order
    ))

    # Write index
    .zarr_write_strings(
      dir = file.path(zarr_group_path, "_index"),
      strings = index,
      compressor = compressor,
      attrs = list(`encoding-type` = "string-array", `encoding-version` = "0.2.0")
    )

    # Cache legacy categories
    has_cats <- h5_group$exists("__categories")
    cats <- if (has_cats) h5_group[["__categories"]] else NULL
    cat_names <- if (has_cats) names(cats) else character(0)

    # Stream each column
    for (col in col_order) {
      if (!h5_group$exists(col)) next
      col_obj <- h5_group[[col]]
      col_dir <- file.path(zarr_group_path, col)

      if (inherits(col_obj, "H5Group")) {
        # Modern categorical
        enc_type <- tryCatch(hdf5r::h5attr(col_obj, "encoding-type"),
                             error = function(e) "")
        if (enc_type == "categorical" &&
            col_obj$exists("categories") && col_obj$exists("codes")) {
          codes <- col_obj[["codes"]]$read()
          categories <- as.character(col_obj[["categories"]]$read())
          .zarr_create_group(col_dir, attrs = list(
            `encoding-type` = "categorical",
            `encoding-version` = "0.2.0",
            ordered = FALSE
          ))
          .zarr_write_numeric(
            dir = file.path(col_dir, "codes"),
            data = as.integer(codes),
            dtype = "|i1",
            compressor = compressor
          )
          .zarr_write_strings(
            dir = file.path(col_dir, "categories"),
            strings = categories,
            compressor = compressor,
            attrs = list(`encoding-type` = "string-array",
                         `encoding-version` = "0.2.0")
          )
        }
      } else if (inherits(col_obj, "H5D")) {
        if (has_cats && col %in% cat_names) {
          # Legacy categorical
          codes <- col_obj$read()
          categories <- as.character(cats[[col]]$read())
          .zarr_create_group(col_dir, attrs = list(
            `encoding-type` = "categorical",
            `encoding-version` = "0.2.0",
            ordered = FALSE
          ))
          .zarr_write_numeric(
            dir = file.path(col_dir, "codes"),
            data = as.integer(codes),
            dtype = "|i1",
            compressor = compressor
          )
          .zarr_write_strings(
            dir = file.path(col_dir, "categories"),
            strings = categories,
            compressor = compressor,
            attrs = list(`encoding-type` = "string-array",
                         `encoding-version` = "0.2.0")
          )
        } else {
          vals <- col_obj$read()
          if (is.character(vals)) {
            .zarr_write_strings(
              dir = col_dir, strings = vals, compressor = compressor,
              attrs = list(`encoding-type` = "string-array",
                           `encoding-version` = "0.2.0")
            )
          } else if (is.logical(vals)) {
            .zarr_write_numeric(dir = col_dir, data = as.integer(vals),
                                dtype = "|b1", compressor = compressor)
          } else if (is.integer(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<i4", compressor = compressor)
          } else if (is.numeric(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<f8", compressor = compressor)
          }
        }
      }
    }
  }

  # Helper: stream a matrix (sparse or dense) from h5ad to zarr
  .stream_h5_matrix_to_zarr <- function(h5_obj, zarr_mat_path, transpose) {
    if (inherits(h5_obj, "H5Group") &&
        h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
      # Sparse
      encoding <- tryCatch(hdf5r::h5attr(h5_obj, "encoding-type"),
                           error = function(e) "csr_matrix")
      shape <- tryCatch(hdf5r::h5attr(h5_obj, "shape"), error = function(e) NULL)

      data_vals <- h5_obj[["data"]][]
      indices <- h5_obj[["indices"]][]
      indptr <- h5_obj[["indptr"]][]

      is_csr <- identical(encoding, "csr_matrix")
      is_csc <- identical(encoding, "csc_matrix")

      if (is.null(shape)) {
        if (is_csr) {
          shape <- c(length(indptr) - 1L,
                     if (length(indices) > 0) max(indices) + 1L else 0L)
        } else {
          shape <- c(if (length(indices) > 0) max(indices) + 1L else 0L,
                     length(indptr) - 1L)
        }
      }

      if (transpose && is_csr) {
        # CSR cells×genes -> reinterpret as CSC genes×cells -> write back as CSR cells×genes
        # Just copy directly — the data is already in cells×genes CSR
        out_encoding <- "csr_matrix"
        out_shape <- as.integer(shape)
      } else if (!transpose) {
        # Keep original encoding (obsp graphs)
        out_encoding <- encoding
        out_shape <- as.integer(shape)
      } else {
        # CSC + transpose: need to re-encode as CSR
        # Build dgCMatrix from CSC, then DeconstructSparseCSR for transposed CSR
        mat <- Matrix::sparseMatrix(
          i = as.integer(indices) + 1L,
          p = as.integer(indptr),
          x = as.numeric(data_vals),
          dims = as.integer(shape),
          repr = "C"
        )
        csr <- DeconstructSparseCSR(mat)
        data_vals <- csr$data
        indices <- csr$indices
        indptr <- csr$indptr
        out_shape <- as.integer(csr$shape)
        out_encoding <- "csr_matrix"
      }

      .zarr_create_group(zarr_mat_path, attrs = list(
        `encoding-type` = out_encoding,
        `encoding-version` = "0.1.0",
        shape = as.list(out_shape)
      ))
      .zarr_parallel_write(
        write_specs = list(
          list(dir = file.path(zarr_mat_path, "data"),
               data = as.numeric(data_vals), dtype = "<f8"),
          list(dir = file.path(zarr_mat_path, "indices"),
               data = as.integer(indices), dtype = "<i4"),
          list(dir = file.path(zarr_mat_path, "indptr"),
               data = as.integer(indptr), dtype = "<i4")
        ),
        compressor = compressor
      )
    } else if (inherits(h5_obj, "H5D")) {
      # Dense — check dimensionality: 1D datasets need [] not [,]
      ndims <- length(h5_obj$dims)
      mat <- if (ndims == 1L) {
        matrix(h5_obj[], ncol = 1L)
      } else {
        h5_obj[,]
      }
      if (transpose) mat <- t(mat)
      .zarr_write_numeric(
        dir = zarr_mat_path, data = mat, dtype = "<f8",
        compressor = compressor,
        attrs = list(`encoding-type` = "array", `encoding-version` = "0.2.0")
      )
    }
  }

  # 1. obs
  if (verbose) message("Streaming obs...")
  if (h5$exists("obs")) {
    .stream_h5_df_to_zarr(h5[["obs"]], file.path(zarr_path, "obs"))
  } else {
    .zarr_create_group(file.path(zarr_path, "obs"))
  }

  # 2. var
  if (verbose) message("Streaming var...")
  if (h5$exists("var")) {
    .stream_h5_df_to_zarr(h5[["var"]], file.path(zarr_path, "var"))
  } else {
    .zarr_create_group(file.path(zarr_path, "var"))
  }

  # 3. X
  if (h5$exists("X")) {
    if (verbose) message("Streaming X...")
    .stream_h5_matrix_to_zarr(h5[["X"]], file.path(zarr_path, "X"),
                               transpose = FALSE)
  }

  # 4. layers
  .zarr_create_group(file.path(zarr_path, "layers"))
  if (h5$exists("layers")) {
    for (ln in names(h5[["layers"]])) {
      if (verbose) message("  Streaming layer: ", ln)
      tryCatch({
        .stream_h5_matrix_to_zarr(
          h5[["layers"]][[ln]],
          file.path(zarr_path, "layers", ln),
          transpose = FALSE
        )
      }, error = function(e) {
        if (verbose) warning("Could not stream layer ", ln, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 5. obsm
  .zarr_create_group(file.path(zarr_path, "obsm"))
  if (h5$exists("obsm") && inherits(h5[["obsm"]], "H5Group")) {
    for (rn in names(h5[["obsm"]])) {
      if (verbose) message("  Streaming obsm: ", rn)
      tryCatch({
        item <- h5[["obsm"]][[rn]]
        if (inherits(item, "H5Group")) {
          # Sparse embedding (rare)
          .stream_h5_matrix_to_zarr(item, file.path(zarr_path, "obsm", rn),
                                     transpose = FALSE)
        } else {
          # Dense: read and write as-is (already cells × dims in h5ad)
          ndims_item <- length(item$dims)
          emb <- if (ndims_item == 1L) {
            matrix(item[], ncol = 1L)
          } else {
            item[,]
          }
          # h5ad stores (n_obs, n_dims) but hdf5r may return transposed
          if (ndims_item > 1L && ncol(emb) > nrow(emb) && item$dims[1] > item$dims[2]) {
            emb <- t(emb)
          }
          .zarr_write_numeric(
            dir = file.path(zarr_path, "obsm", rn),
            data = emb, dtype = "<f8", compressor = compressor
          )
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream obsm ", rn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 6. obsp
  .zarr_create_group(file.path(zarr_path, "obsp"))
  if (h5$exists("obsp") && inherits(h5[["obsp"]], "H5Group")) {
    for (gn in names(h5[["obsp"]])) {
      if (verbose) message("  Streaming obsp: ", gn)
      tryCatch({
        .stream_h5_matrix_to_zarr(h5[["obsp"]][[gn]],
                                   file.path(zarr_path, "obsp", gn),
                                   transpose = FALSE)
      }, error = function(e) {
        if (verbose) warning("Could not stream obsp ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 7. uns
  .zarr_create_group(file.path(zarr_path, "uns"))
  if (h5$exists("uns") && inherits(h5[["uns"]], "H5Group")) {
    for (un in names(h5[["uns"]])) {
      tryCatch({
        item <- h5[["uns"]][[un]]
        if (inherits(item, "H5D")) {
          vals <- item[]
          if (is.character(vals)) {
            .zarr_write_strings(
              dir = file.path(zarr_path, "uns", un),
              strings = vals, compressor = compressor
            )
          } else if (is.numeric(vals) || is.logical(vals)) {
            .zarr_write_numeric(
              dir = file.path(zarr_path, "uns", un),
              data = vals, compressor = compressor
            )
          }
        }
        # Skip complex groups silently
      }, error = function(e) {
        if (verbose) warning("Could not stream uns ", un, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 8. varp
  if (h5$exists("varp") && inherits(h5[["varp"]], "H5Group")) {
    if (verbose) message("Streaming varp...")
    .zarr_create_group(file.path(zarr_path, "varp"))
    for (vn in names(h5[["varp"]])) {
      tryCatch({
        .stream_h5_matrix_to_zarr(h5[["varp"]][[vn]],
                                   file.path(zarr_path, "varp", vn),
                                   transpose = FALSE)
      }, error = function(e) {
        if (verbose) warning("Could not stream varp ", vn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # 9. raw (if present, write as raw/ group with X and var)
  if (h5$exists("raw")) {
    if (verbose) message("Streaming raw...")
    raw_grp <- h5[["raw"]]
    .zarr_create_group(file.path(zarr_path, "raw"))
    if (raw_grp$exists("X")) {
      .stream_h5_matrix_to_zarr(raw_grp[["X"]],
                                 file.path(zarr_path, "raw", "X"),
                                 transpose = FALSE)
    }
    if (raw_grp$exists("var")) {
      .stream_h5_df_to_zarr(raw_grp[["var"]], file.path(zarr_path, "raw", "var"))
    }
  }

  if (verbose) message("Streaming h5ad -> zarr complete: ", zarr_path)
  invisible(zarr_path)
}

#' Convert a zarr store to h5ad format
#'
#' Converts an AnnData zarr store to h5ad (HDF5) format. By default uses
#' streaming conversion that copies fields directly without loading into R
#' memory.
#'
#' @param source Path to input .zarr directory
#' @param dest Path for output .h5ad file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly without creating
#'   a Seurat intermediate. If FALSE, load as Seurat then save as h5ad.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToH5AD <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!dir.exists(source))
    stop("Zarr store not found: ", source, call. = FALSE)
  if (file.exists(dest) || dir.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_zarr_to_h5ad(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readZarr(source, verbose = verbose)
    writeH5AD(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_zarr_to_h5ad <- function(zarr_path, h5ad_path, gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for ZarrToH5AD", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for ZarrToH5AD", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5ad_path, mode = "w")
  on.exit(h5$close_all())

  # Root attrs
  h5$create_attr(attr_name = "encoding-type", robj = "anndata",
                 dtype = CachedGuessDType("anndata"), space = ScalarSpace())
  h5$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                 dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())

  # Helper: stream a zarr DataFrame to h5ad
  .stream_zarr_df_to_h5 <- function(store_path, zarr_group, h5_parent, h5_name) {
    grp_path <- zarr_group
    if (.zarr_node_type(store_path, grp_path) != "group") {
      # Write empty group
      g <- h5_parent$create_group(h5_name)
      g$create_attr(attr_name = "encoding-type", robj = "dataframe",
                    dtype = CachedGuessDType("dataframe"), space = ScalarSpace())
      g$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                    dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
      g$create_attr(attr_name = "_index", robj = "_index",
                    dtype = CachedGuessDType("_index"), space = ScalarSpace())
      g$create_attr(attr_name = "column-order", robj = character(0),
                    dtype = CachedUtf8Type())
      return()
    }

    attrs <- .zarr_read_attrs(store_path, grp_path)
    index <- .zarr_read_anndata_index(store_path, zarr_group)
    col_order <- attrs[["column-order"]]
    if (is.null(col_order)) {
      children <- .zarr_list_children(store_path, grp_path)
      col_order <- setdiff(children, c("_index", "index"))
    }

    # Build a data.frame for WriteDFGroup
    df <- data.frame(row.names = if (!is.null(index)) index else seq_along(index))
    for (col in col_order) {
      col_path <- file.path(grp_path, col)
      if (.zarr_node_type(store_path, col_path) == "missing") next
      tryCatch({
        vals <- .zarr_read_anndata_column(store_path, col_path)
        if (!is.null(vals)) df[[col]] <- vals
      }, error = function(e) {
        if (verbose) warning("Could not read column ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }

    idx <- if (!is.null(index)) index else rownames(df)
    WriteDFGroup(h5_parent, h5_name, df, idx, gzip = gzip)
  }

  # Helper: stream a zarr matrix to h5ad
  .stream_zarr_matrix_to_h5 <- function(store_path, zarr_mat_path,
                                         h5_parent, h5_name) {
    node_type <- .zarr_node_type(store_path, zarr_mat_path)
    if (node_type == "missing") return()

    if (node_type == "group") {
      attrs <- .zarr_read_attrs(store_path, zarr_mat_path)
      enc <- attrs[["encoding-type"]] %||% ""

      if (enc %in% c("csr_matrix", "csc_matrix")) {
        data_vals <- .zarr_read_numeric(store_path,
                                         file.path(zarr_mat_path, "data"))
        indices <- .zarr_read_numeric(store_path,
                                       file.path(zarr_mat_path, "indices"))
        indptr <- .zarr_read_numeric(store_path,
                                      file.path(zarr_mat_path, "indptr"))
        shape <- attrs[["shape"]]
        if (is.list(shape)) shape <- unlist(shape)

        chunk_size <- 65536L
        grp <- h5_parent$create_group(h5_name)
        grp$create_dataset("data", robj = as.numeric(data_vals),
                           chunk_dims = min(length(data_vals), chunk_size),
                           gzip_level = gzip)
        grp$create_dataset("indices", robj = as.integer(indices),
                           chunk_dims = min(length(indices), chunk_size),
                           gzip_level = gzip)
        grp$create_dataset("indptr", robj = as.integer(indptr),
                           chunk_dims = min(length(indptr), chunk_size),
                           gzip_level = gzip)
        grp$create_attr(attr_name = "encoding-type", robj = enc,
                        dtype = CachedGuessDType(enc), space = ScalarSpace())
        grp$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                        dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())
        grp$create_attr(attr_name = "shape", robj = as.integer(shape),
                        dtype = GuessDType(as.integer(shape)))
      }
    } else if (node_type == "array") {
      mat <- .zarr_read_numeric(store_path, zarr_mat_path)
      if (is.matrix(mat)) {
        h5_parent$create_dataset(h5_name, robj = mat,
                                 chunk_dims = dim(mat), gzip_level = gzip)
      } else {
        h5_parent$create_dataset(h5_name, robj = mat,
                                 chunk_dims = length(mat), gzip_level = gzip)
      }
      ds <- h5_parent[[h5_name]]
      ds$create_attr(attr_name = "encoding-type", robj = "array",
                     dtype = CachedGuessDType("array"), space = ScalarSpace())
      ds$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                     dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
    }
  }

  # 1. obs
  if (verbose) message("Streaming obs...")
  .stream_zarr_df_to_h5(zarr_path, "obs", h5, "obs")

  # 2. var
  if (verbose) message("Streaming var...")
  .stream_zarr_df_to_h5(zarr_path, "var", h5, "var")

  # 3. X
  if (.zarr_node_type(zarr_path, "X") != "missing") {
    if (verbose) message("Streaming X...")
    .stream_zarr_matrix_to_h5(zarr_path, "X", h5, "X")
  }

  # 4. layers
  if (.zarr_node_type(zarr_path, "layers") == "group") {
    layer_names <- .zarr_list_children(zarr_path, "layers")
    if (length(layer_names) > 0) {
      layers_grp <- h5$create_group("layers")
      for (ln in layer_names) {
        if (verbose) message("  Streaming layer: ", ln)
        tryCatch({
          .stream_zarr_matrix_to_h5(zarr_path, file.path("layers", ln),
                                     layers_grp, ln)
        }, error = function(e) {
          if (verbose) warning("Could not stream layer ", ln, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 5. obsm
  if (.zarr_node_type(zarr_path, "obsm") == "group") {
    obsm_names <- .zarr_list_children(zarr_path, "obsm")
    if (length(obsm_names) > 0) {
      obsm_grp <- h5$create_group("obsm")
      for (rn in obsm_names) {
        if (verbose) message("  Streaming obsm: ", rn)
        tryCatch({
          node <- .zarr_node_type(zarr_path, file.path("obsm", rn))
          if (node == "group") {
            .stream_zarr_matrix_to_h5(zarr_path, file.path("obsm", rn),
                                       obsm_grp, rn)
          } else if (node == "array") {
            emb <- .zarr_read_numeric(zarr_path, file.path("obsm", rn))
            if (is.matrix(emb)) {
              obsm_grp$create_dataset(rn, robj = emb,
                                      chunk_dims = dim(emb), gzip_level = gzip)
            }
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream obsm ", rn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 6. obsp
  if (.zarr_node_type(zarr_path, "obsp") == "group") {
    obsp_names <- .zarr_list_children(zarr_path, "obsp")
    if (length(obsp_names) > 0) {
      obsp_grp <- h5$create_group("obsp")
      for (gn in obsp_names) {
        if (verbose) message("  Streaming obsp: ", gn)
        tryCatch({
          .stream_zarr_matrix_to_h5(zarr_path, file.path("obsp", gn),
                                     obsp_grp, gn)
        }, error = function(e) {
          if (verbose) warning("Could not stream obsp ", gn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 7. uns
  if (.zarr_node_type(zarr_path, "uns") == "group") {
    uns_names <- .zarr_list_children(zarr_path, "uns")
    if (length(uns_names) > 0) {
      uns_grp <- h5$create_group("uns")
      for (un in uns_names) {
        tryCatch({
          node <- .zarr_node_type(zarr_path, file.path("uns", un))
          if (node == "array") {
            meta <- .zarr_read_json(file.path(zarr_path, "uns", un, ".zarray"))
            if (!is.null(meta$dtype) && meta$dtype == "|O") {
              vals <- .zarr_read_strings(zarr_path, file.path("uns", un))
              uns_grp$create_dataset(un, robj = vals,
                                     dtype = CachedUtf8Type(),
                                     chunk_dims = length(vals), gzip_level = gzip)
            } else {
              vals <- .zarr_read_numeric(zarr_path, file.path("uns", un))
              uns_grp$create_dataset(un, robj = vals,
                                     chunk_dims = length(vals), gzip_level = gzip)
            }
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream uns ", un, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 8. varp
  if (.zarr_node_type(zarr_path, "varp") == "group") {
    varp_names <- .zarr_list_children(zarr_path, "varp")
    if (length(varp_names) > 0) {
      varp_grp <- h5$create_group("varp")
      for (vn in varp_names) {
        if (verbose) message("  Streaming varp: ", vn)
        tryCatch({
          .stream_zarr_matrix_to_h5(zarr_path, file.path("varp", vn),
                                     varp_grp, vn)
        }, error = function(e) {
          if (verbose) warning("Could not stream varp ", vn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # 9. raw
  if (.zarr_node_type(zarr_path, "raw") == "group") {
    if (verbose) message("Streaming raw...")
    raw_grp <- h5$create_group("raw")
    if (.zarr_node_type(zarr_path, "raw/X") != "missing") {
      .stream_zarr_matrix_to_h5(zarr_path, "raw/X", raw_grp, "X")
    }
    if (.zarr_node_type(zarr_path, "raw/var") == "group") {
      .stream_zarr_df_to_h5(zarr_path, "raw/var", raw_grp, "var")
    }
  }

  if (verbose) message("Streaming zarr -> h5ad complete: ", h5ad_path)
  invisible(h5ad_path)
}

#' Convert an h5seurat file to zarr (AnnData) format
#'
#' Converts an h5seurat file to AnnData zarr format, translating the h5seurat
#' layout (levels/values factors, genes x cells CSC, reductions as
#' cell.embeddings) to AnnData conventions (categories/codes, cells x genes
#' CSR, obsm arrays). By default uses streaming that avoids a Seurat
#' intermediate.
#'
#' @param source Path to input .h5seurat file
#' @param dest Path for output .zarr directory
#' @param assay Assay name (default "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as zarr.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5SeuratToZarr <- function(source, dest, assay = "RNA", overwrite = FALSE,
                            gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_h5seurat_to_zarr(source, dest, assay = assay, gzip = gzip,
                              verbose = verbose)
  } else {
    obj <- readH5Seurat(source, verbose = verbose)
    writeZarr(obj, dest, assay = assay, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_h5seurat_to_zarr <- function(h5seurat_path, zarr_path, assay = "RNA",
                                      gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for H5SeuratToZarr", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for H5SeuratToZarr", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5seurat_path, mode = "r")
  on.exit(h5$close_all())

  compressor <- list(id = "zlib", level = as.integer(gzip))

  # Detect active assay
  if (h5$attr_exists("active.assay")) {
    assay <- hdf5r::h5attr(h5, "active.assay")
  }

  # Root zarr group
  .zarr_create_group(zarr_path, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # Read cell names (SeuratDisk may store as Nx1 2D dataset)
  cell_names <- if (h5$exists("cell.names")) {
    cn_ds <- h5[["cell.names"]]
    if (length(cn_ds$dims) == 2L) as.character(cn_ds[, 1]) else as.character(cn_ds[])
  } else {
    NULL
  }

  # Read feature names (Seurat v5 then v4 paths)
  assay_path <- paste0("assays/", assay)
  feature_names <- NULL
  if (h5$exists(assay_path)) {
    ag <- h5[[assay_path]]
    if (ag$exists("features")) {
      feat_ds <- ag[["features"]]
      # SeuratDisk may store features as Nx1 (2D), handle both 1D and 2D
      feature_names <- if (length(feat_ds$dims) == 2L) {
        as.character(feat_ds[, 1])
      } else {
        as.character(feat_ds[])
      }
    }
  }

  n_cells <- length(cell_names)
  n_features <- length(feature_names)

  # --- obs (meta.data → AnnData DataFrame) ---
  if (verbose) message("Streaming obs (meta.data)...")
  if (h5$exists("meta.data")) {
    md <- h5[["meta.data"]]
    md_cols <- character(0)
    # Get column order from colnames attr
    if (md$attr_exists("colnames")) {
      md_cols <- hdf5r::h5attr(md, "colnames")
    } else {
      md_cols <- names(md)
    }

    .zarr_create_group(file.path(zarr_path, "obs"), attrs = list(
      `encoding-type` = "dataframe",
      `encoding-version` = "0.2.0",
      `_index` = "_index",
      `column-order` = md_cols
    ))

    # Write _index
    if (!is.null(cell_names)) {
      .zarr_write_strings(
        dir = file.path(zarr_path, "obs", "_index"),
        strings = cell_names,
        compressor = compressor,
        attrs = list(`encoding-type` = "string-array",
                     `encoding-version` = "0.2.0")
      )
    }

    # Write each column
    for (col in md_cols) {
      if (!md$exists(col)) next
      col_obj <- md[[col]]
      col_dir <- file.path(zarr_path, "obs", col)

      tryCatch({
        if (inherits(col_obj, "H5Group")) {
          # Factor: levels/values (1-based) → categories/codes (0-based)
          if (col_obj$exists("levels") && col_obj$exists("values")) {
            levs <- as.character(col_obj[["levels"]][])
            vals <- col_obj[["values"]][]
            codes <- as.integer(vals) - 1L  # 1-based → 0-based

            .zarr_create_group(col_dir, attrs = list(
              `encoding-type` = "categorical",
              `encoding-version` = "0.2.0",
              ordered = FALSE
            ))
            .zarr_write_numeric(
              dir = file.path(col_dir, "codes"),
              data = codes, dtype = "|i1", compressor = compressor
            )
            .zarr_write_strings(
              dir = file.path(col_dir, "categories"),
              strings = levs, compressor = compressor,
              attrs = list(`encoding-type` = "string-array",
                           `encoding-version` = "0.2.0")
            )
          }
        } else if (inherits(col_obj, "H5D")) {
          vals <- col_obj$read()
          if (is.character(vals)) {
            .zarr_write_strings(
              dir = col_dir, strings = vals, compressor = compressor,
              attrs = list(`encoding-type` = "string-array",
                           `encoding-version` = "0.2.0")
            )
          } else if (is.logical(vals)) {
            .zarr_write_numeric(dir = col_dir, data = as.integer(vals),
                                dtype = "|b1", compressor = compressor)
          } else if (is.integer(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<i4", compressor = compressor)
          } else if (is.numeric(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<f8", compressor = compressor)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream obs column ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  } else {
    .zarr_create_group(file.path(zarr_path, "obs"), attrs = list(
      `encoding-type` = "dataframe", `encoding-version` = "0.2.0",
      `_index` = "_index", `column-order` = character(0)
    ))
    if (!is.null(cell_names)) {
      .zarr_write_strings(
        dir = file.path(zarr_path, "obs", "_index"),
        strings = cell_names, compressor = compressor,
        attrs = list(`encoding-type` = "string-array",
                     `encoding-version` = "0.2.0")
      )
    }
  }

  # --- var (features → AnnData DataFrame) ---
  if (verbose) message("Streaming var...")
  var_cols <- character(0)
  if (h5$exists("var")) {
    var_grp <- h5[["var"]]
    if (var_grp$attr_exists("colnames")) {
      var_cols <- hdf5r::h5attr(var_grp, "colnames")
    } else {
      var_cols <- names(var_grp)
    }
  }
  .zarr_create_group(file.path(zarr_path, "var"), attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = var_cols
  ))
  if (!is.null(feature_names)) {
    .zarr_write_strings(
      dir = file.path(zarr_path, "var", "_index"),
      strings = feature_names, compressor = compressor,
      attrs = list(`encoding-type` = "string-array",
                   `encoding-version` = "0.2.0")
    )
  }
  # Write var metadata columns
  if (h5$exists("var") && length(var_cols) > 0) {
    var_grp <- h5[["var"]]
    for (col in var_cols) {
      if (!var_grp$exists(col)) next
      col_obj <- var_grp[[col]]
      col_dir <- file.path(zarr_path, "var", col)
      tryCatch({
        if (inherits(col_obj, "H5Group") &&
            col_obj$exists("levels") && col_obj$exists("values")) {
          levs <- as.character(col_obj[["levels"]][])
          vals <- col_obj[["values"]][]
          codes <- as.integer(vals) - 1L
          .zarr_create_group(col_dir, attrs = list(
            `encoding-type` = "categorical", `encoding-version` = "0.2.0",
            ordered = FALSE))
          .zarr_write_numeric(dir = file.path(col_dir, "codes"),
                              data = codes, dtype = "|i1", compressor = compressor)
          .zarr_write_strings(dir = file.path(col_dir, "categories"),
                              strings = levs, compressor = compressor,
                              attrs = list(`encoding-type` = "string-array",
                                           `encoding-version` = "0.2.0"))
        } else if (inherits(col_obj, "H5D")) {
          vals <- col_obj$read()
          if (is.character(vals)) {
            .zarr_write_strings(dir = col_dir, strings = vals,
                                compressor = compressor,
                                attrs = list(`encoding-type` = "string-array",
                                             `encoding-version` = "0.2.0"))
          } else if (is.logical(vals)) {
            .zarr_write_numeric(dir = col_dir, data = as.integer(vals),
                                dtype = "|b1", compressor = compressor)
          } else if (is.integer(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<i4", compressor = compressor)
          } else if (is.numeric(vals)) {
            .zarr_write_numeric(dir = col_dir, data = vals,
                                dtype = "<f8", compressor = compressor)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream var column ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- X (expression matrix: genes×cells CSC → cells×genes CSR) ---
  # Try v5 path first, then v4
  x_path <- NULL
  for (p in c(paste0(assay_path, "/layers/data"),
              paste0(assay_path, "/data"))) {
    if (h5$exists(p)) { x_path <- p; break }
  }
  if (!is.null(x_path)) {
    if (verbose) message("Streaming X (", x_path, ")...")
    x_obj <- h5[[x_path]]
    .stream_h5seurat_matrix_to_zarr(x_obj, file.path(zarr_path, "X"),
                                     compressor, transpose = TRUE)
  }

  # --- layers (counts → raw/X if separate from data) ---
  counts_path <- NULL
  for (p in c(paste0(assay_path, "/layers/counts"),
              paste0(assay_path, "/counts"))) {
    if (h5$exists(p)) { counts_path <- p; break }
  }

  # If no data layer exists but counts does, also write counts as X
  # so that anndata.read_zarr() puts data in adata.X (not just adata.raw.X)
  if (is.null(x_path) && !is.null(counts_path)) {
    if (verbose) message("Streaming X (from counts: ", counts_path, ")...")
    counts_obj <- h5[[counts_path]]
    .stream_h5seurat_matrix_to_zarr(counts_obj, file.path(zarr_path, "X"),
                                     compressor, transpose = TRUE)
  }

  .zarr_create_group(file.path(zarr_path, "layers"))
  if (!is.null(counts_path) && !identical(counts_path, x_path)) {
    if (verbose) message("Streaming raw/X (", counts_path, ")...")
    .zarr_create_group(file.path(zarr_path, "raw"))
    counts_obj <- h5[[counts_path]]
    .stream_h5seurat_matrix_to_zarr(counts_obj, file.path(zarr_path, "raw", "X"),
                                     compressor, transpose = TRUE)
    # Write raw/var with just _index
    .zarr_create_group(file.path(zarr_path, "raw", "var"), attrs = list(
      `encoding-type` = "dataframe", `encoding-version` = "0.2.0",
      `_index` = "_index", `column-order` = character(0)
    ))
    if (!is.null(feature_names)) {
      .zarr_write_strings(
        dir = file.path(zarr_path, "raw", "var", "_index"),
        strings = feature_names, compressor = compressor,
        attrs = list(`encoding-type` = "string-array",
                     `encoding-version` = "0.2.0")
      )
    }
  }

  # --- obsm (reductions → obsm/X_{name}) ---
  .zarr_create_group(file.path(zarr_path, "obsm"))
  if (h5$exists("reductions")) {
    reduc_grp <- h5[["reductions"]]
    for (rn in names(reduc_grp)) {
      if (verbose) message("  Streaming obsm: X_", rn)
      tryCatch({
        r_obj <- reduc_grp[[rn]]
        if (r_obj$exists("cell.embeddings")) {
          emb <- r_obj[["cell.embeddings"]]
          if (inherits(emb, "H5D")) {
            ndims_emb <- length(emb$dims)
            mat <- if (ndims_emb == 1L) {
              matrix(emb[], ncol = 1L)
            } else {
              emb[,]
            }
            # hdf5r reads h5seurat cell.embeddings as [n_cells, n_comp]
            # zarr AnnData also needs [n_cells, n_comp] — no transpose needed
            # Only fix if dimensions are unexpectedly swapped
            if (nrow(mat) != n_cells && ncol(mat) == n_cells) {
              mat <- t(mat)
            }
            .zarr_write_numeric(
              dir = file.path(zarr_path, "obsm", paste0("X_", rn)),
              data = mat, dtype = "<f8", compressor = compressor
            )
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream reduction ", rn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- obsp (graphs → obsp/{name}) ---
  .zarr_create_group(file.path(zarr_path, "obsp"))
  if (h5$exists("graphs")) {
    graph_grp <- h5[["graphs"]]
    for (gn in names(graph_grp)) {
      if (verbose) message("  Streaming obsp: ", gn)
      tryCatch({
        g_obj <- graph_grp[[gn]]
        if (inherits(g_obj, "H5Group") &&
            g_obj$exists("data") && g_obj$exists("indices") && g_obj$exists("indptr")) {
          # h5seurat graphs are CSC (dgCMatrix layout)
          data_vals <- g_obj[["data"]][]
          indices <- g_obj[["indices"]][]
          indptr <- g_obj[["indptr"]][]
          shape <- if (g_obj$attr_exists("dims")) {
            hdf5r::h5attr(g_obj, "dims")
          } else {
            c(length(indptr) - 1L,
              if (length(indices) > 0) max(indices) + 1L else 0L)
          }

          .zarr_create_group(file.path(zarr_path, "obsp", gn), attrs = list(
            `encoding-type` = "csc_matrix",
            `encoding-version` = "0.1.0",
            shape = as.list(as.integer(shape))
          ))
          .zarr_parallel_write(
            write_specs = list(
              list(dir = file.path(zarr_path, "obsp", gn, "data"),
                   data = as.numeric(data_vals), dtype = "<f8"),
              list(dir = file.path(zarr_path, "obsp", gn, "indices"),
                   data = as.integer(indices), dtype = "<i4"),
              list(dir = file.path(zarr_path, "obsp", gn, "indptr"),
                   data = as.integer(indptr), dtype = "<i4")
            ),
            compressor = compressor
          )
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream graph ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- uns (misc → uns) ---
  .zarr_create_group(file.path(zarr_path, "uns"))
  if (h5$exists("misc")) {
    misc_grp <- h5[["misc"]]
    for (un in names(misc_grp)) {
      tryCatch({
        item <- misc_grp[[un]]
        if (inherits(item, "H5D")) {
          vals <- item[]
          if (is.character(vals)) {
            .zarr_write_strings(dir = file.path(zarr_path, "uns", un),
                                strings = vals, compressor = compressor)
          } else if (is.numeric(vals) || is.logical(vals)) {
            .zarr_write_numeric(dir = file.path(zarr_path, "uns", un),
                                data = vals, compressor = compressor)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream misc ", un, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  if (verbose) message("Streaming h5seurat -> zarr complete: ", zarr_path)
  invisible(zarr_path)
}

#' Stream h5seurat sparse matrix to zarr
#'
#' h5seurat stores CSC (genes x cells). If transpose=TRUE, reinterprets as
#' CSR (cells x genes) using the zero-copy CSC↔CSR duality.
#'
#' @keywords internal
#' @noRd
.stream_h5seurat_matrix_to_zarr <- function(h5_obj, zarr_mat_path,
                                             compressor, transpose = TRUE) {
  if (inherits(h5_obj, "H5Group") &&
      h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
    # Sparse CSC: i=row indices (0-based), p=col pointers, x=values
    data_vals <- h5_obj[["data"]][]
    indices <- h5_obj[["indices"]][]
    indptr <- h5_obj[["indptr"]][]
    shape <- if (h5_obj$attr_exists("dims")) {
      hdf5r::h5attr(h5_obj, "dims")
    } else {
      c(if (length(indices) > 0) max(indices) + 1L else 0L,
        length(indptr) - 1L)
    }

    if (transpose) {
      # CSC of (genes × cells) == CSR of (cells × genes) — zero-copy
      out_shape <- as.integer(rev(shape))  # [n_cells, n_genes]
      .zarr_create_group(zarr_mat_path, attrs = list(
        `encoding-type` = "csr_matrix",
        `encoding-version` = "0.1.0",
        shape = as.list(out_shape)
      ))
    } else {
      out_shape <- as.integer(shape)
      .zarr_create_group(zarr_mat_path, attrs = list(
        `encoding-type` = "csc_matrix",
        `encoding-version` = "0.1.0",
        shape = as.list(out_shape)
      ))
    }
    .zarr_parallel_write(
      write_specs = list(
        list(dir = file.path(zarr_mat_path, "data"),
             data = as.numeric(data_vals), dtype = "<f8"),
        list(dir = file.path(zarr_mat_path, "indices"),
             data = as.integer(indices), dtype = "<i4"),
        list(dir = file.path(zarr_mat_path, "indptr"),
             data = as.integer(indptr), dtype = "<i4")
      ),
      compressor = compressor
    )
  } else if (inherits(h5_obj, "H5D")) {
    ndims <- length(h5_obj$dims)
    if (ndims == 1L) {
      # 1D dense: likely an embedding or a single feature.
      # Cannot reliably determine (rows, cols) from a 1D vector alone,
      # so skip rather than write a wrong-shape matrix.
      warning("Skipping 1D dense dataset in H5SeuratToZarr: ",
              "cannot determine matrix dimensions from 1D data")
    } else {
      mat <- h5_obj[,]
      if (transpose) mat <- t(mat)
      .zarr_write_numeric(
        dir = zarr_mat_path, data = mat, dtype = "<f8",
        compressor = compressor,
        attrs = list(`encoding-type` = "array", `encoding-version` = "0.2.0")
      )
    }
  }
}

#' Convert a zarr (AnnData) store to h5seurat format
#'
#' Converts an AnnData zarr store to h5seurat format, translating AnnData
#' conventions (categories/codes 0-based, cells x genes CSR, obsm arrays)
#' to h5seurat layout (levels/values 1-based, genes x cells CSC,
#' reductions/cell.embeddings). By default uses streaming that avoids a
#' Seurat intermediate.
#'
#' @param source Path to input .zarr directory
#' @param dest Path for output .h5seurat file
#' @param assay Assay name (default "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as h5seurat.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToH5Seurat <- function(source, dest, assay = "RNA", overwrite = FALSE,
                            gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!dir.exists(source))
    stop("Zarr store not found: ", source, call. = FALSE)
  if (file.exists(dest) || dir.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    .stream_zarr_to_h5seurat(source, dest, assay = assay, gzip = gzip,
                              verbose = verbose)
  } else {
    obj <- readZarr(source, verbose = verbose)
    writeH5Seurat(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' @keywords internal
.stream_zarr_to_h5seurat <- function(zarr_path, h5seurat_path, assay = "RNA",
                                      gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for ZarrToH5Seurat", call. = FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("jsonlite required for ZarrToH5Seurat", call. = FALSE)

  h5 <- hdf5r::H5File$new(h5seurat_path, mode = "w")
  on.exit(h5$close_all())

  # Read cell and feature names from zarr
  cell_names <- .zarr_read_anndata_index(zarr_path, "obs")
  feature_names <- .zarr_read_anndata_index(zarr_path, "var")
  n_cells <- length(cell_names)
  n_features <- length(feature_names)

  # Root attributes
  h5$create_attr(attr_name = "active.assay", robj = assay,
                 dtype = CachedGuessDType(assay), space = ScalarSpace())
  h5$create_attr(attr_name = "project", robj = "scConvert",
                 dtype = CachedGuessDType("scConvert"), space = ScalarSpace())
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("SeuratObject")),
    error = function(e) "5.0.0"
  )
  h5$create_attr(attr_name = "version", robj = pkg_ver,
                 dtype = CachedGuessDType(pkg_ver), space = ScalarSpace())

  # --- cell.names ---
  if (!is.null(cell_names)) {
    h5$create_dataset("cell.names", robj = cell_names,
                      dtype = CachedUtf8Type(),
                      chunk_dims = length(cell_names), gzip_level = gzip)
  }

  # --- meta.data (obs → meta.data with levels/values factors) ---
  if (verbose) message("Streaming meta.data (obs)...")
  obs_attrs <- .zarr_read_attrs(zarr_path, "obs")
  obs_cols <- obs_attrs[["column-order"]]
  if (is.null(obs_cols)) {
    obs_cols <- setdiff(.zarr_list_children(zarr_path, "obs"),
                        c("_index", "index"))
  }

  md_grp <- h5$create_group("meta.data")
  col_written <- character(0)
  for (col in obs_cols) {
    col_path <- file.path("obs", col)
    if (.zarr_node_type(zarr_path, col_path) == "missing") next
    tryCatch({
      node_type <- .zarr_node_type(zarr_path, col_path)
      if (node_type == "group") {
        attrs <- .zarr_read_attrs(zarr_path, col_path)
        enc <- attrs[["encoding-type"]] %||% ""
        if (enc == "categorical") {
          codes <- .zarr_read_numeric(zarr_path, file.path(col_path, "codes"))
          # Try reading categories as strings first
          cats_path <- file.path(col_path, "categories")
          cats_meta <- .zarr_read_json(file.path(zarr_path, cats_path, ".zarray"))
          cats <- if (!is.null(cats_meta$dtype) && cats_meta$dtype == "|O") {
            .zarr_read_strings(zarr_path, cats_path)
          } else {
            tryCatch(.zarr_read_strings(zarr_path, cats_path),
                     error = function(e) {
                       as.character(.zarr_read_numeric(zarr_path, cats_path))
                     })
          }
          # Write as h5seurat factor: levels + values (1-based)
          fac_grp <- md_grp$create_group(col)
          fac_grp$create_dataset("levels", robj = cats,
                                 dtype = CachedUtf8Type(),
                                 chunk_dims = length(cats), gzip_level = gzip)
          values_1based <- as.integer(codes) + 1L
          fac_grp$create_dataset("values", robj = values_1based,
                                 chunk_dims = length(values_1based),
                                 gzip_level = gzip)
          col_written <- c(col_written, col)
        }
      } else if (node_type == "array") {
        meta <- .zarr_read_json(file.path(zarr_path, col_path, ".zarray"))
        if (!is.null(meta$dtype) && meta$dtype == "|O") {
          vals <- .zarr_read_strings(zarr_path, col_path)
          md_grp$create_dataset(col, robj = vals, dtype = CachedUtf8Type(),
                                chunk_dims = length(vals), gzip_level = gzip)
        } else {
          vals <- .zarr_read_numeric(zarr_path, col_path)
          if (is.logical(vals)) vals <- as.integer(vals)
          md_grp$create_dataset(col, robj = vals,
                                chunk_dims = length(vals), gzip_level = gzip)
        }
        col_written <- c(col_written, col)
      }
    }, error = function(e) {
      if (verbose) warning("Could not stream obs column ", col, ": ",
                           e$message, immediate. = TRUE)
    })
  }
  if (length(col_written) > 0) {
    md_grp$create_attr(attr_name = "colnames", robj = col_written,
                       dtype = GuessDType(col_written))
  }

  # --- assays/{assay} ---
  assay_grp <- h5$create_group("assays")
  a_grp <- assay_grp$create_group(assay)
  a_grp$create_attr(attr_name = "key", robj = paste0(tolower(assay), "_"),
                    dtype = CachedGuessDType(paste0(tolower(assay), "_")),
                    space = ScalarSpace())
  a_grp$create_attr(attr_name = "s4class", robj = "SeuratObject::Assay5",
                    dtype = CachedGuessDType("SeuratObject::Assay5"),
                    space = ScalarSpace())

  # features
  if (!is.null(feature_names)) {
    a_grp$create_dataset("features", robj = feature_names,
                         dtype = CachedUtf8Type(),
                         chunk_dims = length(feature_names), gzip_level = gzip)
  }

  # --- X → layers/data (transpose CSR cells×genes → CSC genes×cells) ---
  layers_grp <- a_grp$create_group("layers")
  if (.zarr_node_type(zarr_path, "X") != "missing") {
    if (verbose) message("Streaming X -> layers/data...")
    .stream_zarr_matrix_to_h5seurat(zarr_path, "X", layers_grp, "data",
                                     gzip, transpose = TRUE)
  }

  # --- raw/X or layers/counts → layers/counts ---
  counts_src <- NULL
  if (.zarr_node_type(zarr_path, "raw/X") != "missing") {
    counts_src <- "raw/X"
  } else if (.zarr_node_type(zarr_path, "layers") == "group") {
    lc <- .zarr_list_children(zarr_path, "layers")
    if ("counts" %in% lc) counts_src <- "layers/counts"
  }
  if (!is.null(counts_src)) {
    if (verbose) message("Streaming ", counts_src, " -> layers/counts...")
    .stream_zarr_matrix_to_h5seurat(zarr_path, counts_src, layers_grp,
                                     "counts", gzip, transpose = TRUE)
  }

  # --- reductions (obsm → reductions) ---
  reductions_grp <- h5$create_group("reductions")
  if (.zarr_node_type(zarr_path, "obsm") == "group") {
    for (rn in .zarr_list_children(zarr_path, "obsm")) {
      clean_name <- gsub("^X_", "", rn)
      if (verbose) message("  Streaming reduction: ", clean_name)
      tryCatch({
        r_grp <- reductions_grp$create_group(clean_name)
        r_grp$create_attr(attr_name = "active.assay", robj = assay,
                          dtype = CachedGuessDType(assay), space = ScalarSpace())
        key <- AnnDataReductionKey(clean_name)
        r_grp$create_attr(attr_name = "key", robj = key,
                          dtype = CachedGuessDType(key), space = ScalarSpace())
        r_grp$create_attr(attr_name = "global", robj = 0L,
                          dtype = GuessDType(0L), space = ScalarSpace())
        r_grp$create_group("misc")

        node_type <- .zarr_node_type(zarr_path, file.path("obsm", rn))
        if (node_type == "array") {
          emb <- .zarr_read_numeric(zarr_path, file.path("obsm", rn))
          # zarr stores [n_cells, n_comp] in C order.
          # hdf5r writes R matrices in Fortran order, so writing [n_cells, n_comp]
          # from R produces the same HDF5 layout as writeH5Seurat (which writes
          # Embeddings() directly — also [n_cells, n_comp] in R).
          # No transpose needed.
          if (is.matrix(emb) && nrow(emb) != n_cells && ncol(emb) == n_cells) {
            emb <- t(emb)  # fix only if dimensions are swapped
          }
          r_grp$create_dataset("cell.embeddings", robj = emb,
                               chunk_dims = dim(emb), gzip_level = gzip)
        } else if (node_type == "group") {
          # Sparse embedding — read and densify
          mat <- .zarr_read_anndata_matrix(zarr_path, file.path("obsm", rn),
                                            transpose = FALSE)
          mat <- as.matrix(mat)
          if (nrow(mat) != n_cells && ncol(mat) == n_cells) mat <- t(mat)
          r_grp$create_dataset("cell.embeddings", robj = mat,
                               chunk_dims = dim(mat), gzip_level = gzip)
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream reduction ", clean_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- graphs (obsp → graphs) ---
  graphs_grp <- h5$create_group("graphs")
  if (.zarr_node_type(zarr_path, "obsp") == "group") {
    for (gn in .zarr_list_children(zarr_path, "obsp")) {
      if (verbose) message("  Streaming graph: ", gn)
      tryCatch({
        gp <- file.path("obsp", gn)
        attrs <- .zarr_read_attrs(zarr_path, gp)
        enc <- attrs[["encoding-type"]] %||% ""

        if (enc %in% c("csr_matrix", "csc_matrix")) {
          data_vals <- .zarr_read_numeric(zarr_path, file.path(gp, "data"))
          indices <- .zarr_read_numeric(zarr_path, file.path(gp, "indices"))
          indptr <- .zarr_read_numeric(zarr_path, file.path(gp, "indptr"))
          shape <- attrs[["shape"]]
          if (is.list(shape)) shape <- unlist(shape)

          if (enc == "csr_matrix") {
            # CSR → CSC: need actual conversion via dgCMatrix
            mat <- Matrix::sparseMatrix(
              i = rep(seq_len(shape[1]), diff(as.integer(indptr))),
              j = as.integer(indices) + 1L,
              x = as.numeric(data_vals),
              dims = as.integer(shape)
            )
            mat <- as(mat, "dgCMatrix")
            g_grp <- graphs_grp$create_group(gn)
            g_grp$create_dataset("data", robj = mat@x,
                                 chunk_dims = min(length(mat@x), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indices", robj = mat@i,
                                 chunk_dims = min(length(mat@i), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indptr", robj = mat@p,
                                 chunk_dims = min(length(mat@p), 65536L),
                                 gzip_level = gzip)
            g_grp$create_attr(attr_name = "dims", robj = dim(mat),
                              dtype = GuessDType(dim(mat)))
          } else {
            # CSC — direct copy (same as h5seurat native format)
            g_grp <- graphs_grp$create_group(gn)
            g_grp$create_dataset("data", robj = as.numeric(data_vals),
                                 chunk_dims = min(length(data_vals), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indices", robj = as.integer(indices),
                                 chunk_dims = min(length(indices), 65536L),
                                 gzip_level = gzip)
            g_grp$create_dataset("indptr", robj = as.integer(indptr),
                                 chunk_dims = min(length(indptr), 65536L),
                                 gzip_level = gzip)
            g_grp$create_attr(attr_name = "dims", robj = as.integer(shape),
                              dtype = GuessDType(as.integer(shape)))
          }
          g_grp$create_attr(attr_name = "assay.used", robj = assay,
                            dtype = CachedGuessDType(assay), space = ScalarSpace())
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream graph ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- misc (uns → misc) ---
  misc_grp <- h5$create_group("misc")
  if (.zarr_node_type(zarr_path, "uns") == "group") {
    for (un in .zarr_list_children(zarr_path, "uns")) {
      tryCatch({
        node <- .zarr_node_type(zarr_path, file.path("uns", un))
        if (node == "array") {
          meta <- .zarr_read_json(file.path(zarr_path, "uns", un, ".zarray"))
          if (!is.null(meta$dtype) && meta$dtype == "|O") {
            vals <- .zarr_read_strings(zarr_path, file.path("uns", un))
            misc_grp$create_dataset(un, robj = vals, dtype = CachedUtf8Type(),
                                    chunk_dims = length(vals), gzip_level = gzip)
          } else {
            vals <- .zarr_read_numeric(zarr_path, file.path("uns", un))
            misc_grp$create_dataset(un, robj = vals,
                                    chunk_dims = length(vals), gzip_level = gzip)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream uns ", un, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- Required empty groups ---
  for (g in c("tools", "commands", "images", "neighbors")) {
    h5$create_group(g)
  }

  # --- active.ident (default: all cells = single identity) ---
  ident_grp <- h5$create_group("active.ident")
  ident_grp$create_dataset("levels", robj = "scConvert",
                           dtype = CachedUtf8Type(),
                           chunk_dims = 1L, gzip_level = gzip)
  ident_grp$create_dataset("values", robj = rep(1L, n_cells),
                           chunk_dims = n_cells, gzip_level = gzip)

  if (verbose) message("Streaming zarr -> h5seurat complete: ", h5seurat_path)
  invisible(h5seurat_path)
}

#' Stream zarr sparse or dense matrix to h5seurat format
#'
#' AnnData stores cells x genes (CSR). h5seurat needs genes x cells (CSC).
#' CSR of (cells x genes) == CSC of (genes x cells) — zero-copy reinterpretation.
#'
#' @keywords internal
#' @noRd
.stream_zarr_matrix_to_h5seurat <- function(store_path, zarr_mat_path,
                                             h5_parent, h5_name, gzip,
                                             transpose = TRUE) {
  node_type <- .zarr_node_type(store_path, zarr_mat_path)
  if (node_type == "missing") return()

  chunk_size <- 65536L

  if (node_type == "group") {
    attrs <- .zarr_read_attrs(store_path, zarr_mat_path)
    enc <- attrs[["encoding-type"]] %||% ""

    if (enc %in% c("csr_matrix", "csc_matrix")) {
      data_vals <- .zarr_read_numeric(store_path,
                                       file.path(zarr_mat_path, "data"))
      indices <- .zarr_read_numeric(store_path,
                                     file.path(zarr_mat_path, "indices"))
      indptr <- .zarr_read_numeric(store_path,
                                    file.path(zarr_mat_path, "indptr"))
      shape <- attrs[["shape"]]
      if (is.list(shape)) shape <- unlist(shape)

      if (transpose && enc == "csr_matrix") {
        # CSR(cells × genes) == CSC(genes × cells) — zero-copy
        dims <- as.integer(rev(shape))
      } else if (transpose && enc == "csc_matrix") {
        # CSC already, but we need to transpose → build dgCMatrix + reinterpret
        mat <- Matrix::sparseMatrix(
          i = as.integer(indices) + 1L, p = as.integer(indptr),
          x = as.numeric(data_vals), dims = as.integer(shape), repr = "C"
        )
        mat <- Matrix::t(mat)
        data_vals <- mat@x; indices <- mat@i; indptr <- mat@p
        dims <- dim(mat)
      } else {
        dims <- as.integer(shape)
      }

      grp <- h5_parent$create_group(h5_name)
      grp$create_dataset("data", robj = as.numeric(data_vals),
                         chunk_dims = min(length(data_vals), chunk_size),
                         gzip_level = gzip)
      grp$create_dataset("indices", robj = as.integer(indices),
                         chunk_dims = min(length(indices), chunk_size),
                         gzip_level = gzip)
      grp$create_dataset("indptr", robj = as.integer(indptr),
                         chunk_dims = min(length(indptr), chunk_size),
                         gzip_level = gzip)
      grp$create_attr(attr_name = "dims", robj = as.integer(dims),
                      dtype = GuessDType(as.integer(dims)))
    }
  } else if (node_type == "array") {
    mat <- .zarr_read_numeric(store_path, zarr_mat_path)
    if (is.matrix(mat) && transpose) mat <- t(mat)
    h5_parent$create_dataset(h5_name, robj = mat,
                             chunk_dims = if (is.matrix(mat)) dim(mat) else length(mat),
                             gzip_level = gzip)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CLI-accelerated conversion helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Find the scconvert CLI binary
#'
#' Searches for the compiled scconvert binary in the package source tree
#' and on the system PATH.
#'
#' @return Path to the scconvert binary, or NULL if not found
#'
#' @keywords internal
#'
sc_find_cli <- function() {
  # 1. Check package src/ directory (development / devtools::load_all)
  pkg_src <- system.file("src", "scconvert", package = "scConvert")
  if (nzchar(pkg_src) && file.exists(pkg_src)) {
    return(pkg_src)
  }

  # 2. Check relative to package root (devtools::load_all from source tree)
  pkg_dir <- find.package("scConvert", quiet = TRUE)
  if (length(pkg_dir) > 0) {
    src_bin <- file.path(pkg_dir, "src", "scconvert")
    if (file.exists(src_bin)) {
      return(src_bin)
    }
  }

  # 3. Check system PATH
  on_path <- Sys.which("scconvert")
  if (nzchar(on_path)) {
    return(unname(on_path))
  }

  return(NULL)
}

#' Run the scconvert CLI for file-to-file conversion
#'
#' Performs file-to-file format conversion with three tiers:
#' \enumerate{
#'   \item \strong{C binary}: HDF5 format pairs (h5ad, h5seurat, h5mu)
#'   \item \strong{Streaming}: zarr pairs (h5ad/h5seurat \eqn{\leftrightarrow}{<->} zarr) — copies
#'     fields directly without creating a Seurat intermediate
#'   \item \strong{R hub}: all other format pairs via \code{scConvert()}
#' }
#'
#' Supported formats: h5ad, h5seurat, h5mu, loom, rds, zarr.
#'
#' @param input Input file path
#' @param output Output file path
#' @param assay Assay name (passed as --assay)
#' @param gzip Gzip compression level (0-9, passed as --gzip)
#' @param overwrite If TRUE, remove output file before conversion
#' @param verbose If TRUE, print CLI output
#'
#' @return TRUE on success, FALSE on failure
#'
#' @export
#'
scConvert_cli <- function(
  input,
  output,
  assay = "RNA",
  gzip = 4L,
  overwrite = FALSE,
  verbose = TRUE
) {
  is_uri <- grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", input)
  if (!is_uri && !file.exists(input) && !dir.exists(input)) {
    stop("Input file not found: ", input, call. = FALSE)
  }

  stype <- FileType(file = input)
  dtype <- FileType(file = output)

  # Tier 1: Try C binary first for all supported HDF5/loom/zarr format pairs
  cli_formats <- c("h5ad", "h5seurat", "h5mu", "loom", "zarr")
  use_c_cli <- stype %in% cli_formats && dtype %in% cli_formats

  if (use_c_cli) {
    cli_bin <- sc_find_cli()
    if (!is.null(cli_bin)) {
      if (file.exists(output) || dir.exists(output)) {
        if (overwrite) {
          unlink(output, recursive = TRUE)
        } else {
          stop("Output exists: ", output, ". Use overwrite = TRUE.", call. = FALSE)
        }
      }

      args <- c(input, output, "--assay", assay, "--gzip", as.character(gzip))
      if (verbose) message("Using scconvert CLI for fast on-disk conversion")

      result <- tryCatch({
        ret <- system2(cli_bin, args = args, stdout = TRUE, stderr = TRUE)
        status <- attr(ret, "status")
        if (!is.null(status) && status != 0L) {
          if (verbose) {
            message("CLI conversion failed (exit ", status, "): ",
                    paste(ret, collapse = "\n"))
          }
          FALSE
        } else {
          if (verbose && length(ret) > 0) message(paste(ret, collapse = "\n"))
          file.exists(output)
        }
      }, error = function(e) {
        if (verbose) message("CLI binary error: ", e$message)
        FALSE
      })

      if (!result && (file.exists(output) || dir.exists(output))) {
        try(unlink(output, recursive = TRUE), silent = TRUE)
      }
      if (result) return(TRUE)
      if (verbose) message("C binary failed, falling back to R streaming")
    }
  }

  # Tier 2: R streaming converters (zarr pairs + h5mu/loom pairs)
  streaming_converters <- list(
    # zarr pairs
    "h5ad|zarr"      = function() H5ADToZarr(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|h5ad"      = function() ZarrToH5AD(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5seurat|zarr"  = function() H5SeuratToZarr(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|h5seurat"  = function() ZarrToH5Seurat(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # h5mu ↔ zarr
    "h5mu|zarr"      = function() H5MUToZarr(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|h5mu"      = function() ZarrToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom ↔ zarr
    "loom|zarr"      = function() LoomToZarr(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|loom"      = function() ZarrToLoom(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # h5mu ↔ h5seurat (R streaming)
    "h5mu|h5seurat"  = function() H5MUToH5Seurat(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5seurat|h5mu"  = function() H5SeuratToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom ↔ h5seurat (R streaming)
    "loom|h5seurat"  = function() LoomToH5Seurat(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5seurat|loom"  = function() H5SeuratToLoom(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # h5mu ↔ h5ad (R streaming via h5seurat)
    "h5mu|h5ad"      = function() H5MUToH5AD(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5ad|h5mu"      = function() H5ADToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom ↔ h5ad (R streaming via h5seurat)
    "loom|h5ad"      = function() LoomToH5AD(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5ad|loom"      = function() H5ADToLoom(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom ↔ h5mu (R streaming via h5seurat)
    "loom|h5mu"      = function() LoomToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5mu|loom"      = function() H5MUToLoom(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose)
  )

  pair_key <- paste0(stype, "|", dtype)
  if (pair_key %in% names(streaming_converters)) {
    result <- tryCatch({
      streaming_converters[[pair_key]]()
      TRUE
    }, error = function(e) {
      if (verbose) message("R streaming failed: ", e$message,
                           " — falling back to R hub")
      FALSE
    })
    if (result) return(TRUE)
  }

  # R-level conversion: supports all format pairs
  if (verbose) message("Converting ", stype, " -> ", dtype, " via R")
  tryCatch({
    scConvert(source = input, dest = output, assay = assay,
            overwrite = overwrite, verbose = verbose)
    TRUE
  }, error = function(e) {
    if (verbose) message("Conversion failed: ", e$message)
    FALSE
  })
}
