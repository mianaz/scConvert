
#' @param p A vector of sparse matrix pointers
#'
#' @return \code{PointerToIndex}: A vector of column (j) indices
#'
#' @rdname SparsePointers
#'
#' @keywords internal
#'
#' @source \code{PointerToIndex} came from
#' \href{https://stackoverflow.com/questions/20008200/r-constructing-sparse-matrix}{StackOverflow}
#' @author \code{PointerToIndex} was written by
#' \href{https://stackoverflow.com/users/980833/josh-obrien}{Josh O'Brien on StackOverflow}
#'
PointerToIndex <- function(p) {
  dp <- diff(x = p)
  j <- rep.int(x = seq_along(along.with = dp), times = dp)
  return(j)
}
