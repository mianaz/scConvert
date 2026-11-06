#' @include H5Utils.R
#' @importFrom hdf5r h5attr
#' @importFrom Matrix sparseMatrix
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Coercing hdf5r objects to R objects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Read HDF5 datasets and groups into R objects
#'
#' Coercion methods for \pkg{hdf5r} handles, following the h5Seurat
#' conventions: data frames are compound datasets or groups of equal-length
#' columns (factors as \code{values}/\code{levels} groups, logicals as
#' integers listed in a \code{logicals} attribute), and sparse matrices are
#' CSC groups (\code{data}, \code{indices}, \code{indptr}) with a \code{dims}
#' attribute.
#'
#' @param x An \code{H5D} dataset or \code{H5Group}
#' @param row.names Row names for the data frame
#' @param optional Keep the stored column names as they are instead of
#'   passing them through \code{\link[base]{make.names}}
#' @param recursive Also convert nested groups (otherwise they are returned
#'   as \code{H5Group} handles)
#' @param transpose Transpose the matrix after reading
#' @param ... Passed on (\code{as.array}) or ignored
#'
#' @return The corresponding array, data frame, list, matrix or
#'   \code{dgCMatrix}
#'
#' @name ReadH5
#' @rdname ReadH5
#' @keywords internal
#'
NULL

# Read a dataset without dropping length-1 axes
.read_keep_dims <- function(x) {
  if (length(x$dims) < 2L) return(x$read())
  tryCatch(x[, , drop = FALSE], error = function(e) x$read())
}

# Decode h5Seurat logical columns (written as 0 = FALSE, 1 = TRUE,
# 2 = NA); as.logical() alone would turn the NA code into TRUE
IntToBool <- function(x) {
  if (is.logical(x)) return(x)
  codes <- as.integer(x)
  ifelse(codes %in% 0:1, codes == 1L, NA)
}

# Restore the logical columns listed in the `logicals` attribute of `x`
.restore_logicals <- function(df, x) {
  if (x$attr_exists(attr_name = 'logicals')) {
    flagged <- intersect(as.character(h5attr(x = x, which = 'logicals')), names(df))
    for (col in flagged) df[[col]] <- IntToBool(df[[col]])
  }
  df
}

# Is the `ordered` attribute of an h5Seurat factor group set?
.is_ordered_group <- function(g) {
  g$attr_exists(attr_name = 'ordered') && isTRUE(as.logical(h5attr(x = g, which = 'ordered')))
}

# Factor from an h5Seurat factor group (1-based `values` + `levels`);
# out-of-range codes become NA
.factor_from_group <- function(g, order_by = NULL) {
  lv <- g[['levels']][]
  codes <- as.integer(g[['values']][])
  codes[!is.na(codes) & (codes < 1L | codes > length(lv))] <- NA_integer_
  if (!is.null(order_by)) codes <- codes[order(order_by)]
  factor(lv[codes], levels = lv, ordered = .is_ordered_group(g))
}

#' @rdname ReadH5
#' @method as.array H5D
#' @export
#'
as.array.H5D <- function(x, ...) {
  as.array(.read_keep_dims(x), ...)
}

#' @rdname ReadH5
#' @method as.data.frame H5D
#' @export
#'
as.data.frame.H5D <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- if (length(x$dims) >= 2L) .read_keep_dims(x) else tryCatch(x[], error = function(e) x$read())
  if (!is.null(row.names)) row.names(df) <- row.names
  if (!optional) {
    labels <- tryCatch(x$get_type()$get_cpd_labels(), error = function(e) NULL)
    if (!is.null(labels)) try(colnames(df) <- make.names(labels), silent = TRUE)
  }
  .restore_logicals(df, x)
}

#' @rdname ReadH5
#' @method as.data.frame H5Group
#' @export
#'
as.data.frame.H5Group <- function(x, row.names = NULL, optional = FALSE, ...) {
  members <- names(x)
  # an `index` dataset gives the order the rows were written in
  order_by <- if ('index' %in% members) x[['index']][] else NULL
  columns <- list()
  for (col in setdiff(members, 'index')) {
    obj <- x[[col]]
    value <- tryCatch({
      if (inherits(obj, 'H5Group')) {
        if (IsFactor(obj)) .factor_from_group(obj, order_by) else NULL
      } else {
        v <- as.vector(obj[])
        if (!is.null(order_by)) v <- v[order(order_by)]
        v
      }
    }, error = function(e) NULL)   # unreadable columns are left out
    if (!is.null(value)) columns[[col]] <- value
  }
  n <- if (length(columns)) length(columns[[1]]) else 0L
  df <- data.frame(row.names = row.names %||% if (n) seq_len(n) else character(0))
  for (col in names(columns)) {
    # a column whose length does not match the first one is left out
    df <- tryCatch({ df[[col]] <- columns[[col]]; df }, error = function(e) df)
  }
  df <- .restore_logicals(df, x)
  # `colnames` records the original column order
  if (x$attr_exists(attr_name = 'colnames') && ncol(df) > 0L) {
    keep <- intersect(h5attr(x = x, which = 'colnames'), names(df))
    if (length(keep)) df <- df[, keep, drop = FALSE]
  }
  df
}

#' @rdname ReadH5
#' @method as.list H5D
#' @export
#'
as.list.H5D <- function(x, ...) {
  as.list(x[])
}

#' @rdname ReadH5
#' @method as.list H5Group
#' @export
#'
as.list.H5Group <- function(x, recursive = TRUE, ...) {
  out <- list()
  for (nm in names(x)) {
    obj <- x[[nm]]
    if (inherits(obj, 'H5Group')) {
      out[[nm]] <- if (isTRUE(recursive)) as.list(obj, recursive = TRUE, ...) else obj
      next
    }
    kind <- GetStringType(obj)
    out[[nm]] <- if (identical(kind, 'H5T_STRING')) {
      as.character(obj[])
    } else if (identical(kind, 'numeric')) {
      as.vector(obj[])
    } else {
      obj[]
    }
  }
  out
}

#' @rdname ReadH5
#' @method as.matrix H5D
#' @export
#'
as.matrix.H5D <- function(x, transpose = FALSE, ...) {
  m <- .read_keep_dims(x)
  as.matrix(if (isTRUE(transpose)) t(m) else m)
}

#' @rdname ReadH5
#' @method as.matrix H5Group
#' @export
#'
as.matrix.H5Group <- function(x, ...) {
  as.sparse(x = x, ...)
}

# Matrix dimensions recorded on a sparse dataset/group, in R order
.sparse_dims <- function(x) {
  key <- intersect(c('dims', 'shape'), hdf5r::h5attr_names(x))
  d <- if (length(key)) as.integer(h5attr(x = x, which = key[1])) else NULL
  if (length(d) < 2L || anyNA(d) || any(d <= 0L)) {
    stop("A sparse matrix needs valid dimensions (got: ", paste(d, collapse = " x "), ")",
         call. = FALSE)
  }
  rev(d)
}

#' @importFrom SeuratObject as.sparse
#' @rdname ReadH5
#' @method as.sparse H5D
#' @export
#'
as.sparse.H5D <- function(x, ...) {
  # a coordinate table: one row per entry, (row, column, value)
  d <- .sparse_dims(x)
  tryCatch({
    triplets <- x$read()
    sparseMatrix(i = triplets[, 1], j = triplets[, 2], x = triplets[, 3], dims = d)
  }, error = function(e) {
    stop("Could not build a sparse matrix: ", conditionMessage(e), call. = FALSE)
  })
}

#' @rdname ReadH5
#' @method as.sparse H5Group
#' @export
#'
as.sparse.H5Group <- function(x, ...) {
  if (!all(c('data', 'indices', 'indptr') %in% names(x))) {
    stop("Not a sparse matrix group", call. = FALSE)
  }
  d <- .sparse_dims(x)
  tryCatch(
    sparseMatrix(i = as.integer(x[['indices']][] + 1), p = as.integer(x[['indptr']][]),
                 x = as.vector(x[['data']][]), dims = d),
    error = function(e) {
      stop("Could not build a sparse matrix: ", conditionMessage(e), call. = FALSE)
    }
  )
}

# Coarse type of an HDF5 type, dataset or group: "factor" for h5Seurat
# factor groups, "numeric" for integer/float data, otherwise the HDF5 type
# class (e.g. "H5T_STRING")
GetStringType <- function(stype) {
  if (inherits(stype, c('H5D', 'H5Group'))) {
    if (IsFactor(stype)) return('factor')
    stype <- stype$get_type()
  }
  if (inherits(stype, c('H5T', 'H5T_COMPOUND'))) stype <- stype$get_class()
  stype <- as.character(stype)
  if (stype %in% c('H5T_INTEGER', 'H5T_FLOAT')) 'numeric' else stype
}
