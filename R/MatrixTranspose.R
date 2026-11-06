#' @include zzz.R
#' @importFrom Matrix t sparseMatrix
NULL

#' Transpose a matrix held in memory or in an HDF5 file
#'
#' For an in-memory \code{dgCMatrix} the transpose is returned. For an HDF5
#' dense dataset or a compressed-sparse group, the transposed matrix is
#' written next to the original (or into \code{dest}) under \code{dname}.
#'
#' @param x A \code{dgCMatrix}, a two-dimensional \code{H5D}, or an
#'   \code{H5Group} holding \code{data}/\code{indices}/\code{indptr}
#' @param dest HDF5 group or file that receives the result
#' @param dname Name of the result inside \code{dest}
#' @param overwrite Replace an existing \code{dname} in \code{dest}
#' @param verbose Report progress
#' @param ... Unused
#'
#' @return The transposed \code{dgCMatrix}, or (invisibly) the new HDF5
#'   dataset/group
#'
#' @name scTranspose
#' @rdname scTranspose
#' @keywords internal
#' @export
#'
scTranspose <- function(x, ...) {
  UseMethod(generic = 'scTranspose', object = x)
}

#' @rdname scTranspose
#' @method scTranspose dgCMatrix
#' @export
#'
scTranspose.dgCMatrix <- function(x, ...) {
  Matrix::t(x)
}

# Resolve the target location shared by the HDF5 methods
.transpose_target <- function(x, dest, dname, overwrite) {
  if (is.null(dest)) {
    parent <- dirname(x$get_obj_name())
    dest <- x$get_file_id()[[if (identical(parent, ".")) "/" else parent]]
  }
  if (dest$exists(dname)) {
    if (!isTRUE(overwrite)) {
      stop("'", dname, "' already exists; set overwrite = TRUE", call. = FALSE)
    }
    dest$link_delete(dname)
  }
  dest
}

.default_tname <- function(x) paste0("t_", basename(x$get_obj_name()))

#' @rdname scTranspose
#' @method scTranspose H5D
#' @export
#'
scTranspose.H5D <- function(x, dest = NULL, dname = .default_tname(x),
                            overwrite = FALSE, verbose = TRUE, ...) {
  d <- x$dims
  if (length(d) != 2L) {
    stop("Only two-dimensional datasets can be transposed", call. = FALSE)
  }
  # when writing over the source itself, materialise it first
  same_place <- is.null(dest) && identical(dname, basename(x$get_obj_name()))
  values <- if (same_place) x$read() else NULL
  dtype <- x$get_type()
  dest <- .transpose_target(x, dest, dname, overwrite || same_place)
  out <- dest$create_dataset(name = dname, dtype = dtype,
                             space = hdf5r::H5S$new(dims = rev(d), maxdims = rev(d)))
  if (!is.null(values)) {
    out[, ] <- t(values)
    return(invisible(out))
  }
  # copy in column blocks of the source so memory stays bounded
  step <- max(1L, as.integer(floor(1e7 / max(1, d[1]))))
  for (start in seq.int(1L, d[2], by = step)) {
    idx <- start:min(d[2], start + step - 1L)
    if (verbose && d[2] > step) message("Transposing columns ", start, "-", max(idx))
    out[idx, ] <- t(x[, idx, drop = FALSE])
  }
  invisible(out)
}

#' @rdname scTranspose
#' @method scTranspose H5Group
#' @export
#'
scTranspose.H5Group <- function(x, dest = NULL, dname = .default_tname(x),
                                overwrite = FALSE, ...) {
  for (part in c("data", "indices", "indptr")) {
    if (!x$exists(part)) {
      stop("Group is not a compressed sparse matrix (missing '", part, "')",
           call. = FALSE)
    }
  }
  dim_attr <- intersect(c("dims", "shape", "h5sparse_shape"), hdf5r::h5attr_names(x))
  ptr <- x[["indptr"]][]
  idx <- x[["indices"]][]
  dims <- if (length(dim_attr)) as.integer(hdf5r::h5attr(x, dim_attr[1])) else
    c(if (length(idx)) max(idx) + 1L else 0L, length(ptr) - 1L)
  m <- Matrix::sparseMatrix(i = idx + 1L, p = ptr, x = as.numeric(x[["data"]][]),
                            dims = dims)
  tm <- methods::as(Matrix::t(m), "CsparseMatrix")
  same_place <- is.null(dest) && identical(dname, basename(x$get_obj_name()))
  dest <- .transpose_target(x, dest, dname, overwrite || same_place)
  g <- dest$create_group(dname)
  g[["data"]] <- tm@x
  g[["indices"]] <- tm@i
  g[["indptr"]] <- tm@p
  g$create_attr(attr_name = "dims", robj = dim(tm), dtype = GuessDType(x = dim(tm)))
  invisible(g)
}
