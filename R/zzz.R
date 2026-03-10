#' @importFrom rlang %||%
#' @importFrom hdf5r H5T_COMPOUND
#' @importFrom methods setOldClass is
#' @importFrom stats na.omit
#' @importFrom utils data read.csv
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
#'
#' scConvert uses the following options to control behavior, users can configure
#' these with \code{\link[base]{options}}:
#'
#' \describe{
#'  \item{\code{scConvert.dtypes.logical_to_int}}{
#'   When writing \link[base]{logical} vectors, coerce to integer types to
#'   ensure compatibility across languages (see \code{\link{BoolToInt}} for
#'   more details)
#'  }
#'  \item{\code{scConvert.dtypes.dataframe_as_group}}{
#'   When writing \link[base]{data.frame}s, always write out as a group
#'   regardless of factor presence
#'  }
#'  \item{\code{scConvert.chunking.MARGIN}}{
#'   Default direction for chunking datasets; choose from:
#'   \describe{
#'    \item{largest}{Chunk along the largest dimension of a dataset}
#'    \item{smallest}{Chunk along the smallest dimension}
#'    \item{first}{Chunk along the first dimension}
#'    \item{last}{Chunk along the last dimension}
#'   }
#'  }
#'  \item{\code{scConvert.dimreducs.allglobal}}{
#'   Treat all DimReducs as global, regardless of actual global status
#'  }
#'  \item{\code{scConvert.compression.level}}{
#'   Gzip compression level for HDF5 dataset writes (0-9). Level 0 means
#'   no compression, level 4 (default) is a good balance of speed vs size.
#'   Only applies to newly created datasets, not to data copied via
#'   \code{obj_copy_from}.
#'  }
#' }
#'
#' @aliases scConvert
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

modes <- list(
  'new' = c('w', 'w-', 'x'),
  'existing' = c('r', 'r+')
)

version.regex <- '^\\d+(\\.\\d+){2}(\\.9\\d{3})?$'

scdisk.types <- new.env()

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

#' Convert a logical to an integer
#'
#' Unlike most programming languages, R has three possible \link[base]{logical}
#' (boolean) values: \code{TRUE}, \code{FALSE}, and \code{\link[base]{NA}};
#' moreover, the \code{NA} value has representations in other data types, such
#' as \code{NA_integer_}, \code{NA_real_}, and \code{NA_character_}. Simply
#' writing out the logical values to an HDF5 file would cause issues when trying
#' to read the data in to another language, such as Python. To encode these three
#' logical values for other languages, we can encode the logicals as integers:
#' \itemize{
#'  \item \code{FALSE} becomes \code{0L}
#'  \item \code{TRUE} becomes \code{1L}
#'  \item \code{NA} becomes \code{2L}
#' }
#' This encoding scheme allows other languages to handle \code{NA}s in their own
#' manner while preserving all three logicals for R
#'
#' @param x A logical vector
#'
#' @return An integer vector
#'
#' @seealso \link[base]{integer} \link[base]{logical} \code{\link[base]{NA}}
#' \code{\link{scSaveH5Seurat}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::BoolToInt(x = c(TRUE, FALSE, NA))
#' }
#'
BoolToInt <- function(x) {
  x <- as.integer(x = x)
  x[which(x = is.na(x = x))] <- 2L
  return(x)
}

#' Generate chunk points
#'
#' @param dsize Size of data being chunked
#' @param csize Size of chunk; if \code{NA}, assumes single chunk
#'
#' @return A matrix where each row is a chunk, column 1 is start points, column
#' 2 is end points
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::ChunkPoints(100, 3)
#' scConvert:::ChunkPoints(100, NA)
#' }
#'
ChunkPoints <- function(dsize, csize) {
  if (is.na(x = csize)) {
    return(matrix(
      data = c(1, dsize),
      ncol = 2,
      dimnames = list(NULL, c('start', 'end'))
    ))
  }
  return(t(x = vapply(
    X = seq.default(from = 1L, to = ceiling(dsize / csize)),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  )))
}

#' Find the closest version
#'
#' API changes happen at set versions, and knowing how a current running version
#' relates to versions introducing API changes is important.
#' \code{ClosestVersion} approximages both \dQuote{rounding down} (eg. to
#' determine minimum version with new API addition) and \dQuote{rounding up}
#' (eg. to determine maximum version before API deletion) for semantic versions.
#'
#' @param query A query version (\code{\link[base]{character}} or
#' \code{\link[base]{numeric_version}})
#' @param targets A vector of target versions (\code{\link[base]{character}} or
#' \code{\link[base]{numeric_version}})
#' @param direction Which way should we check for closest version? Choose from:
#' \describe{
#'  \item{min}{Closest version less than or equal to \code{query}}
#'  \item{max}{Closest version greater than or equal to \code{query}}
#' }
#' @param inclusive Perform an inclusive comparison (eg. \code{>=} or \code{<=}
#' versus to \code{>} or \code{<}) for \dQuote{rounding}
#'
#' @return The version from \code{targets} that is closest to \code{query} as a
#' \code{\link[base]{character}} vector
#'
#' @seealso \code{\link[base]{numeric_version}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::ClosestVersion('3.1.0', targets = c('3.0.0', '1.4.9', '4.3.2'))
#' scConvert:::ClosestVersion('3.1.0', targets = c('3.0.0', '1.4.9', '4.3.2'), direction = 'max')
#' }
#'
ClosestVersion <- function(
  query,
  targets,
  direction = c('min', 'max'),
  inclusive = direction == 'min'
) {
  direction <- match.arg(arg = direction)
  query <- numeric_version(x = query)
  targets <- sort(x = numeric_version(x = targets))
  switch(
    EXPR = direction,
    'min' = {
      compare <- ifelse(test = inclusive, yes = `<=`, no = `<`)
      collapse <- max
    },
    'max' = {
      compare <- ifelse(test = inclusive, yes = `>=`, no = `>`)
      collapse <- min
    }
  )
  index <- suppressWarnings(expr = collapse(which(x = compare(
    e1 = targets,
    e2 = query
  ))))
  if (is.infinite(x = index)) {
    stop(
      "All target versions ",
      switch(EXPR = direction, 'min' = 'greater', 'max' = 'less'),
      " than query version (",
      as.character(x = query),
      ")",
      call. = FALSE
    )
  }
  return(as.character(x = targets[index]))
}

#' Convert an HDF5 compound dataset to a group (legacy version)
#'
#' @param src An HDF5 dataset
#' @param dest An HDF5 file or group
#' @param dst.name Name of group in \code{dest}
# CompoundToGroup moved to ReadH5.R (single implementation with name_map_fn support)

#' Determine a filetype based on its extension
#'
#' @param file Name of file
#'
#' @return The extension, all lowercase
#'
#' @importFrom tools file_ext
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::FileType('pbmc3k.h5Seurat')
#' scConvert:::FileType('h5ad')
#' }
#'
FileType <- function(file) {
  ext <- file_ext(x = file)
  ext <- ifelse(test = nchar(x = ext), yes = ext, no = basename(path = file))
  return(tolower(x = ext))
}

#' Fix Feature Names
#'
#' @param features A vector of feature names
#'
#' @return Fixed features
#'
#' @keywords internal
#'
FixFeatures <- function(features) {
  # Replace underscores with dashes FIRST (before dedup, to avoid
  # creating new duplicates after make.unique)
  if (any(grepl(pattern = '_', x = features))) {
    message("Feature names cannot have underscores ('_'), replacing with dashes ('-')")
    features <- gsub(pattern = '_', replacement = '-', x = features)
  }
  # Deduplicate: sort duplicate groups by value, then add .1, .2 suffixes
  if (anyDuplicated(x = features)) {
    message("Non-unique features (rownames) present, making unique")
    features <- make.unique(names = features)
  }
  return(features)
}

#' Get a class string with package information
#'
#' S4 classes are useful in the context of their defining package (benefits of
#' stricter typing). In order to ensure class information is properly retained
#' in HDF5 files, S4 class names are written as \dQuote{package:classname} with
#' certain exceptions (eg. S4 classes defined by
#' \link[Seurat:Seurat-package]{Seurat})
#'
#' @param class Class name
#' @param packages A vector of packages to exclude from resulting class
#' information
#'
#' @return A character vector with the class
#'
#' @importFrom methods getClass slot
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::GetClass('Seurat')
#' scConvert:::GetClass('Matrix')
#' }
#'
GetClass <- function(class, packages = 'Seurat') {
  class <- class[1]
  classdef <- getClass(Class = class)
  classpkg <- slot(object = classdef, name = 'package')
  if (classpkg %in% packages) {
    classpkg <- NULL
  }
  class <- paste(classpkg, class, sep = ':')
  return(gsub(pattern = '^:', replacement = '', x = class))
}

#' Determine the margin to use for a dataset
#'
#' @param dims Dimensions of a dataset
#' @param MARGIN Either an integer value contained within
#' \code{1:length(x = dims)} or one of the possible values of
#' the \code{scConvert.chunking.MARGIN} option
#'
#' @return An integer value with the \code{MARGIN}
#'
#' @seealso \code{\link{scConvert-package}} for package options
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::GetMargin(c(4, 10))
#' }
#'
GetMargin <- function(dims, MARGIN = getOption(x = 'scConvert.chunking.MARGIN')) {
  if (isFALSE(x = is.numeric(x = MARGIN))) {
    MARGIN <- tryCatch(
      expr = match.arg(
        arg = MARGIN,
        choices = default.options[['scConvert.chunking.MARGIN']]
      ),
      error = function(err) {
        warning(err$message, call. = FALSE, immediate. = TRUE)
        return(default.options[['scConvert.chunking.MARGIN']][1])
      }
    )
    MARGIN <- switch(
      EXPR = MARGIN,
      'largest' = which.max(x = dims),
      'smallest' = which.min(x = dims),
      'first' = 1L,
      'last' = length(x = dims)
    )
  }
  if (isFALSE(x = MARGIN %in% seq.int(from = 1, to = length(x = dims)))) {
    stop("'MARGIN' must be within the dimensions of the dataset", call. = FALSE)
  }
  return(MARGIN)
}

#' Get the parent of an HDF5 dataset or group
#'
#' @param x An HDF5 dataset or group
#'
#' @return An \code{\link[hdf5r]{H5File}} or \code{\link[hdf5r]{H5Group}} object
#'
#' @keywords internal
#'
GetParent <- function(x) {
  dname <- dirname(path = x$get_obj_name())
  dest <- if (dname == '/') {
    x$get_file_id()
  } else {
    x$get_file_id()[[dname]]
  }
  return(dest)
}

#' Guess an HDF5 Datatype
#'
#' Wrapper around \code{\link[hdf5r:guess_dtype]{hdf5r::guess_dtype}}, allowing
#' for the customization of string types rather than defaulting to
#' variable-length ASCII-encoded strings. Also encodes logicals as
#' \code{\link[hdf5r]{H5T_INTEGER}} instead of \code{\link[hdf5r]{H5T_LOGICAL}}
#' to ensure cross-language compatibility (controlled via
#' \link[=scConvert-package]{package options})
#'
#' @inheritParams StringType
#' @inheritParams hdf5r::guess_dtype
#' @inheritDotParams hdf5r::guess_dtype
#'
#' @return An object of class \code{\link[hdf5r]{H5T}}
#'
#' @importFrom hdf5r guess_dtype
#'
#' @seealso \code{\link[hdf5r]{guess_dtype}} \code{\link{BoolToInt}}
#' \code{\link{StringType}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Characters can either be variable-width UTF8-encoded or
#' # fixed-width ASCII-encoded
#' scConvert:::GuessDType(x = 'hello')
#' scConvert:::GuessDType(x = 'hello', stype = 'ascii7')
#'
#' # Data frames are a compound type; character columns follow the same rules
#' # as character vectors
#' df <- data.frame(x = c('g1', 'g2', 'g3'), y = 1, 2, 3, stringsAsFactors = FALSE)
#' scConvert:::GuessDType(x = df)
#' scConvert:::GuessDType(x = df, stype = 'ascii7')
#'
#' # Logicals are turned into integers to ensure compatibility with Python
#' # TRUE evaluates to 1, FALSE to 0, and NA to 2
#' scConvert:::GuessDType(x = c(TRUE, FALSE, NA))
#' }
#'
GuessDType <- function(x, stype = 'utf8', ...) {
  dtype <- guess_dtype(x = x, ...)
  if (inherits(x = dtype, what = 'H5T_STRING')) {
    dtype <- StringType(stype = stype)
  } else if (inherits(x = dtype, what = 'H5T_COMPOUND')) {
    cpd.dtypes <- dtype$get_cpd_types()
    for (i in seq_along(along.with = cpd.dtypes)) {
      if (inherits(x = cpd.dtypes[[i]], what = 'H5T_STRING')) {
        cpd.dtypes[[i]] <- StringType(stype = stype)
      }
    }
    dtype <- H5T_COMPOUND$new(
      labels = dtype$get_cpd_labels(),
      dtypes = cpd.dtypes,
      size = dtype$get_size()
    )
  } else if (inherits(x = dtype, what = 'H5T_LOGICAL')) {
    if (getOption(x = "scConvert.dtypes.logical_to_int", default = TRUE)) {
      dtype <- guess_dtype(x = BoolToInt(x = x), ...)
    }
  }
  return(dtype)
}

#' Check the datatype of an HDF5 dataset
#'
#' Effectively, an implementation of \code{\link[methods]{is}} for HDF5 datasets;
#' useful to ensure HDF5 validity for specific file structures
#'
#' @param x An HDF5 dataset (object of type \code{\link[hdf5r]{H5D}})
#' @param dtype A character vector of HDF5 datatype names, must be present in
#' \code{\link[hdf5r]{h5types}}
#'
#' @return A logical
#'
#' @importFrom hdf5r h5types
#'
#' @seealso \code{\link[hdf5r]{h5types}}
#'
#' @keywords internal
#'
IsDType <- function(x, dtype) {
  if (!inherits(x = x, what = 'H5D')) {
    stop("'IsDType' only works on HDF5 dataset", call. = FALSE)
  }
  dtypes <- unique(x = sapply(
    X = grep(pattern = '^H5T_', x = names(x = h5types), value = TRUE),
    FUN = function(i) {
      return(class(x = h5types[[i]])[1])
    },
    USE.NAMES = FALSE
  ))
  dtypes <- unique(x = c(dtypes, 'H5T_COMPOUND'))
  match.arg(arg = dtype, choices = dtypes, several.ok = TRUE)
  missing.dtypes <- setdiff(x = dtype, y = dtypes)
  if (length(x = missing.dtypes)) {
    dtype <- setdiff(x = dtype, y = missing.dtypes)
    if (!length(x = dtype)) {
      stop("None of the requested dtypes are valid HDF5 datatypes", call. = FALSE)
    } else {
      warning(
        "The following requested dtypes are not valid HDF5 datatypes: ",
        paste(missing.dtypes, sep = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  return(inherits(x = x$get_type(), what = dtype))
}

#' Check to see if a matrix is empty
#'
#' Determine if a matrix is empty or not. A matrix is considered empty if it
#' satisfies one of the following conditions:
#' \itemize{
#'  \item The dimensions of the matrix are 0-by-0 (\code{all(dim(x) == 0)})
#'  \item The dimensions of the matrix are 1-by-1 (\code{all(dim(x) == 1)}) and
#'  the sole vlaue is \code{NA}
#' }
#' These two situations correspond to matrices generated with either
#' \code{new('matrix')} or \code{matrix()}
#'
#' @param x A matrix
#'
#' @return \code{TRUE} if the matrix is empty otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::IsMatrixEmpty(new('matrix'))
#' scConvert:::IsMatrixEmpty(matrix())
#' scConvert:::IsMatrixEmpty(matrix(1:9, nrow = 3))
#' }
#'
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

#' Make a space
#'
#' Generate a blank space \code{n} characters long; useful for aligning text to
#' be printed to console
#'
#' @param n Length space should be
#'
#' @return A space (' ') of length \code{n}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::MakeSpace(n = 10)
#' cat('hello', scConvert:::MakeSpace(n = 10), 'world\n', sep = '')
#' }
#'
MakeSpace <- function(n) {
  return(paste(rep_len(x = ' ', length.out = n), collapse = ''))
}

#' Add names for unnamed or partially named objects
#'
#' @param x An object that can be named
#' @param prefix A prefix to be added to each name
#'
#' @return \code{x} with unnamed values named
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' a <- list(1, b = 2, 3)
#' scConvert:::PadNames(a)
#' }
#'
PadNames <- function(x, prefix = 'index') {
  if (length(x = x)) {
    xnames <- names(x = x) %||% paste0(
      prefix,
      seq.int(from = 1L, to = length(x = x))
    )
    missing <- which(x = !nchar(x = xnames))
    if (length(x = missing)) {
      xnames[missing] <- paste0(prefix, missing)
    }
    names(x = x) <- xnames
  }
  return(x)
}

#' Create a progress bar
#'
#' Progress bars are useful ways of getting updates on how close a task is to
#' completion. However, they can get in the way of RMarkdown documents with
#' lots of unnecesssary printing. \code{PB} is a convenience function that
#' creates progress bars with the following defaults
#' \itemize{
#'  \item \code{char = '='}
#'  \item \code{style = 3}
#'  \item \code{file = stderr()}
#' }
#'
#' @return An object of class \code{\link[utils]{txtProgressBar}}
#'
#' @importFrom utils txtProgressBar
#'
#' @seealso \code{\link[utils]{txtProgressBar}} \code{\link[base]{stderr}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' pb <- scConvert:::PB()
#' for (i in 1:10) {
#'   utils::setTxtProgressBar(pb, i / 10)
#' }
#' close(pb)
#' }
PB <- function() {
  return(txtProgressBar(char = '=', style = 3, file = stderr()))
}

#' Generate a random string of characters
#'
#' @param length Length (\code{\link[base]{nchar}}) of string to generate
#' @param ... Extra parameters passed to \code{\link[base]{sample}}
#'
#' @return A random string of characters of length (\code{\link[base]{nchar}})
#' of \code{length}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::RandomName()
#' }
#'
RandomName <- function(length = 5L, ...) {
  paste(sample(letters, size = length, ...), collapse = "")
}

#' Generate an HDF5 string dtype
#'
#' Presets for encoding variations of \code{\link[hdf5r]{H5T_STRING}}; used to
#' generate HDF5 datatype specifications with specific string encodings
#'
#' @param stype Type of string encoding to use, choose from:
#' \describe{
#'  \item{utf8}{Variable-width, UTF-8}
#'  \item{ascii7}{Fixed-width (7 bits), ASCII}
#' }
#'
#' @return An \code{\link[hdf5r]{H5T_STRING}} object
#'
#' @importFrom hdf5r h5const H5T_STRING
#'
#' @seealso \code{\link[hdf5r]{H5T_STRING}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::StringType()
#' scConvert:::StringType('ascii7')
#' }
#'
StringType <- function(stype = c('utf8', 'ascii7')) {
  stype <- match.arg(arg = stype)
  return(switch(
    EXPR = stype,
    'utf8' = H5T_STRING$new(size = Inf)$set_cset(cset = h5const$H5T_CSET_UTF8)$set_strpad(strpad = h5const$H5T_STR_NULLTERM),
    'ascii7' = H5T_STRING$new(size = 7L)
  ))
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

#' Update a Seurat key
#'
#' Attempts to validate a string to use as a Seurat key. Valid keys must match
#' the regular expression \code{^[[:alnum:]]+_$}; if \code{key} fails this
#' regular expression, an attempt to modify it to said key will be made by
#' removing all non-alphanumeric characters, collapsing the resulting vector,
#' and appending \dQuote{_}. If this stil fails, a random string of lowercase
#' characters will be generated, followed by \dQuote{_}, to be used as the key
#'
#' @param key A key to validate and update
#'
#' @return \code{key}, updated if invalid
#'
#' @seealso \code{\link[Seurat]{Key}} \code{\link{RandomName}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::UpdateKey("RNA_")
#' scConvert:::UpdateKey("potato")
#' scConvert:::UpdateKey("*@)")
#' }
#'
UpdateKey <- function(key) {
  if (grepl(pattern = "^[[:alnum:]]+_$", x = key)) {
    return(key)
  } else {
    new.key <- regmatches(x = key, m = gregexpr(pattern = "[[:alnum:]]+",text = key))
    new.key <- paste0(paste(unlist(x = new.key), collapse = ""), "_")
    if (new.key == "_") {
      new.key <- paste0(RandomName(length = 3), "_")
    }
    return(new.key)
  }
}

#' Update slots in an object
#'
#' @param object An object to update
#'
#' @return \code{object} with the latest slot definitions
#'
#' @importFrom methods slotNames slot slot<-
#'
#' @keywords internal
#'
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(
        mode = class(x = xobj),
        length = 1L
      )
    }
  }
  return(object)
}

#' Write an attribute to an HDF5 file, group, or dataset (legacy version)
#'
#' @param h5 An HDF5 file, group, or dataset
#' @param name Name to store attribute as
#' @param robj An object to write out
#' @param dtype Data type of attribute
#' @param scalar Is this a scalar or simple (vectorized) attribute?
#' @param overwrite Overwrite the attribute if it already exists
#' @param ... Extra parameters passed to \code{H5S$new}
#'
#' @return Invisibly returns \code{NULL}
#'
#' @importFrom hdf5r H5S
#'
#' @keywords internal
#' @noRd
WriteAttribute <- function(
  h5,
  name,
  robj,
  dtype = GuessDType(x = robj),
  scalar = length(x = robj) == 1,
  overwrite = FALSE,
  ...
) {
  if (!inherits(x = h5, what = c('H5File', 'H5Group', 'H5D'))) {
    stop("'h5' must be an HDF5 file, group, or dataset", call. = FALSE)
  }
  if (h5$attr_exists(attr_name = name)) {
    if (overwrite) {
      h5$attr_delete(attr_name = name)
    } else {
      stop("Attribute ", name, " already exists", call. = FALSE)
    }
  }
  if (is.logical(x = robj) && getOption(x = "scConvert.dtypes.logical_to_int", default = TRUE)) {
    robj <- BoolToInt(x = robj)
  }
  space.type <- ifelse(test = isTRUE(x = scalar), yes = 'scalar', no = 'simple')
  dims <- if (space.type == 'scalar') {
    NULL
  } else {
    dim(x = robj) %||% length(x = robj)
  }
  h5$create_attr(
    attr_name = name,
    robj = robj,
    dtype = dtype,
    space = H5S$new(type = space.type, dims = dims, ...)
  )
  return(invisible(x = NULL))
}

#' Get the proper HDF5 connection mode for writing depending on overwrite status
#'
#' @param overwrite Overwrite a file
#'
#' @return \code{w} if \code{overwrite} else \code{w-}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' scConvert:::WriteMode(TRUE)
#' scConvert:::WriteMode(FALSE)
#' }
#'
WriteMode <- function(overwrite = FALSE) {
  return(ifelse(test = overwrite, yes = 'w', no = 'w-'))
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
#' Checks an internal cache before calling \code{\link{GuessDType}}.
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

#' Safely set layer data in a V5 Assay
#'
#' Handles the case where LayerData<- might not be exported from SeuratObject
#'
#' @param object An Assay or Assay5 object
#' @param layer Layer name
#' @param value Data to set
#' @return Modified assay object
#' @keywords internal
#'
SafeSetLayerData <- function(object, layer, value) {
  # Try using SeuratObject::LayerData<- if available
  tryCatch({
    SeuratObject::LayerData(object, layer = layer) <- value
    return(object)
  }, error = function(e) {
    # Fall back to direct slot assignment
  })

  # Fall back to direct slot assignment
  tryCatch({
    if (inherits(object, "Assay5")) {
      # For Assay5, try to directly modify the layers slot
      if (is.null(object@layers)) {
        object@layers <- list()
      }
      object@layers[[layer]] <- value
    } else {
      # For standard Assay, use standard slots
      slot(object, layer) <- value
    }
    return(object)
  }, error = function(e) {
    warning("Could not set layer data: ", e$message)
    return(object)  # Return unchanged object as last resort
  })
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Format converter functions (defined outside .onLoad to avoid R CMD check NOTE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.h5seurat_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  scLoadH5Seurat(file = file, verbose = verbose, ...)
}

.h5seurat_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  scSaveH5Seurat(object = object, filename = filename,
               overwrite = overwrite, verbose = verbose, ...)
  invisible(filename)
}

.h5ad_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  temp <- tempfile(fileext = '.h5seurat')
  on.exit(unlink(temp), add = TRUE)
  hfile <- scConnect(filename = file, force = TRUE)
  h5s <- H5ADToH5Seurat(source = hfile, dest = temp, assay = assay,
                         overwrite = TRUE, verbose = verbose)
  h5s$close_all()
  hfile$close_all()
  scLoadH5Seurat(file = temp, verbose = verbose)
}

.h5ad_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  assay_name <- DefaultAssay(object = object)

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

    # Standard pipeline with placeholders: Seurat -> h5seurat -> h5ad
    temp <- tempfile(fileext = '.h5Seurat')
    on.exit(unlink(temp), add = TRUE)
    scSaveH5Seurat(object = obj_write, filename = temp, overwrite = TRUE, verbose = FALSE)
    h5s <- scConnect(filename = temp, force = TRUE)
    H5SeuratToH5AD(source = h5s, dest = filename, assay = assay_name,
                    overwrite = overwrite, verbose = verbose, ...)
    h5s$close_all()
    rm(obj_write)

    gzip <- GetCompressionLevel()
    has_data_layer <- "data" %in% names(bpcells_layers)
    for (ln in names(bpcells_layers)) {
      if (ln == "data") {
        group_name <- "X"
      } else if (ln == "counts" && has_data_layer) {
        group_name <- "raw/X"
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

      if (group_name == "raw/X") {
        h5 <- hdf5r::H5File$new(filename, mode = "r+")
        if (h5$exists("raw/var")) {
          h5$link_delete("raw/var")
        }
        h5$create_group("raw/var")
        gene_names <- rownames(bpcells_layers[[ln]])
        h5[["raw/var"]]$create_dataset(
          name = "_index",
          robj = gene_names,
          dtype = StringType('utf8')
        )
        h5[["raw/var"]]$create_attr(
          attr_name = "encoding-type", robj = "dataframe",
          dtype = CachedGuessDType("dataframe"), space = ScalarSpace()
        )
        h5[["raw/var"]]$create_attr(
          attr_name = "encoding-version", robj = "0.2.0",
          dtype = CachedGuessDType("0.2.0"), space = ScalarSpace()
        )
        h5[["raw/var"]]$create_attr(
          attr_name = "_index", robj = "_index",
          dtype = CachedGuessDType("_index"), space = ScalarSpace()
        )
        h5$close_all()
      }

      if (verbose) message("  Streamed '", ln, "' -> ", group_name, " via BPCells (no materialization)")
    }
  } else {
    if (file.exists(filename) && !overwrite) {
      stop("File '", filename, "' already exists; set overwrite = TRUE", call. = FALSE)
    }
    tryCatch({
      DirectSeuratToH5AD(object = object, filename = filename,
                         assay = assay_name, overwrite = overwrite,
                         verbose = verbose, ...)
    }, error = function(e) {
      if (verbose) message("Direct pipeline failed (", e$message, "), falling back to two-pass")
      if (file.exists(filename)) try(file.remove(filename), silent = TRUE)
      temp <- tempfile(fileext = '.h5Seurat')
      on.exit(unlink(temp), add = TRUE)
      scSaveH5Seurat(object = object, filename = temp, overwrite = TRUE, verbose = FALSE)
      h5s <- scConnect(filename = temp, force = TRUE)
      H5SeuratToH5AD(source = h5s, dest = filename, assay = assay_name,
                      overwrite = overwrite, verbose = verbose, ...)
      h5s$close_all()
    })
  }
  invisible(filename)
}

.loom_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  scLoadLoom(file = file, verbose = verbose, ...)
}

.loom_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  scSaveLoom(object = object, filename = filename,
           overwrite = overwrite, verbose = verbose, ...)
  invisible(filename)
}

.h5mu_loader <- function(file, assay = 'RNA', verbose = TRUE, ...) {
  LoadH5MU(file = file, verbose = verbose)
}

.h5mu_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  if (file.exists(filename) && !overwrite) {
    stop("File '", filename, "' already exists; set overwrite = TRUE", call. = FALSE)
  }
  SeuratToH5MU(object = object, filename = filename,
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
  LoadZarr(file = file, assay.name = assay, verbose = verbose, ...)
}

.zarr_saver <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  SaveZarr(object = object, filename = filename,
           overwrite = overwrite, verbose = verbose, ...)
  invisible(filename)
}

.h5ad_to_h5seurat_direct <- function(source, dest, assay = 'RNA',
                                      overwrite = FALSE, verbose = TRUE, ...) {
  if (isTRUE(getOption("scConvert.use_cli"))) {
    source_file <- if (inherits(source, c('H5File', 'h5Seurat'))) source$filename else source
    if (scConvert_cli(input = source_file, output = dest, assay = assay,
                       overwrite = overwrite, verbose = verbose)) {
      return(dest)
    }
    if (verbose) message("CLI unavailable, using R-level conversion")
  }
  H5ADToH5Seurat(source = source, dest = dest, assay = assay,
                 overwrite = overwrite, verbose = verbose)
}

.h5seurat_to_h5ad_direct <- function(source, dest, assay = 'RNA',
                                      overwrite = FALSE, verbose = TRUE,
                                      standardize = FALSE, ...) {
  if (isTRUE(getOption("scConvert.use_cli")) && !standardize) {
    source_file <- if (inherits(source, c('H5File', 'h5Seurat'))) source$filename else source
    if (scConvert_cli(input = source_file, output = dest, assay = assay,
                       overwrite = overwrite, verbose = verbose)) {
      return(dest)
    }
    if (verbose) message("CLI unavailable, using R-level conversion")
  }
  H5SeuratToH5AD(source = source, dest = dest, assay = assay,
                 overwrite = overwrite, verbose = verbose,
                 standardize = standardize)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Loading handler
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  # Make the classes defined in scConvert compatible with S4 generics/methods
  setOldClass(Classes = c('scdisk', 'h5Seurat'))
  RegisterSCDisk(r6class = h5Seurat)
  RegisterSCDisk(r6class = loom)

  # Set some default options
  op <- options()
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }

  # Pre-populate HDF5 caches for common AnnData encoding strings
  .h5_scalar_space <<- H5S$new(type = "scalar")
  .h5_utf8_dtype <<- StringType('utf8')
  for (s in c('string-array', '0.2.0', '0.1.0', 'categorical', 'array',
              'csr_matrix', 'csc_matrix', 'dataframe', 'numeric-scalar',
              '_index')) {
    .h5_guessdtype_cache[[s]] <<- GuessDType(x = s)
  }

  # Register file format loaders and savers
  RegisterFormat(ext = 'h5seurat', loader = .h5seurat_loader, saver = .h5seurat_saver)
  RegisterFormat(ext = 'h5ad', loader = .h5ad_loader, saver = .h5ad_saver)
  RegisterFormat(ext = 'loom', loader = .loom_loader, saver = .loom_saver)
  RegisterFormat(ext = 'h5mu', loader = .h5mu_loader, saver = .h5mu_saver)
  RegisterFormat(ext = 'rds', loader = .rds_loader, saver = .rds_saver)
  RegisterFormat(ext = 'zarr', loader = .zarr_loader, saver = .zarr_saver)

  # Register direct HDF5-level conversion paths
  RegisterDirectPath('h5ad', 'h5seurat', .h5ad_to_h5seurat_direct)
  RegisterDirectPath('h5seurat', 'h5ad', .h5seurat_to_h5ad_direct)

  invisible(x = NULL)
}
