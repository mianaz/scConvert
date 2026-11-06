#' @include Convert.R
NULL

#' Convert many files in parallel
#'
#' Runs \code{\link{scConvert}} on each source file, several at a time. Each
#' conversion runs in its own worker process, so files are converted
#' independently (an HDF5 file can only be written by one process, so
#' parallelism is across files rather than within one). Combine with
#' \code{stream = TRUE} to keep each worker's memory bounded for large
#' files.
#'
#' @param sources Paths of the files to convert
#' @param dest One of: a format name (e.g. \code{"h5seurat"}), which writes
#'   each output next to its source with that extension; a directory, which
#'   requires \code{format}; or a character vector of output paths, one per
#'   source
#' @param format Output format when \code{dest} is a directory
#' @param workers Number of files converted at the same time
#' @param overwrite,verbose,stream,... Passed to \code{\link{scConvert}}
#'
#' @return A data frame with one row per source: \code{source}, \code{dest},
#'   \code{ok}, \code{seconds} and \code{error} (the error message of failed
#'   conversions). Failures do not stop the other conversions.
#'
#' @details Workers are forked processes on Linux and macOS and a socket
#' cluster on Windows. Peak memory is roughly \code{workers} times that of
#' one conversion.
#'
#' @examples
#' \dontrun{
#' files <- list.files("samples", pattern = "\\.h5ad$", full.names = TRUE)
#' res <- scConvertBatch(files, dest = "h5seurat", workers = 4, stream = TRUE)
#' res[!res$ok, ]
#' }
#'
#' @export
#'
scConvertBatch <- function(sources, dest, format = NULL, workers = 1L, overwrite = FALSE,
                           verbose = FALSE, stream = NULL, ...) {
  sources <- as.character(sources)
  if (!length(sources)) stop("No source files given", call. = FALSE)
  outputs <- .batch_outputs(sources, dest, format)
  workers <- max(1L, min(as.integer(workers), length(sources)))
  extra <- list(...)
  one <- function(i) {
    t0 <- Sys.time()
    err <- tryCatch({
      do.call(scConvert, c(list(source = sources[i], dest = outputs[i], overwrite = overwrite,
                                verbose = verbose, stream = stream), extra))
      NA_character_
    }, error = function(e) conditionMessage(e))
    data.frame(source = sources[i], dest = outputs[i], ok = is.na(err),
               seconds = round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 2),
               error = err, stringsAsFactors = FALSE)
  }
  rows <- if (workers == 1L) {
    lapply(seq_along(sources), one)
  } else if (.Platform$OS.type == 'unix') {
    parallel::mclapply(seq_along(sources), one, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(scConvert)))
    parallel::parLapply(cl, seq_along(sources), one)
  }
  # a worker that died (e.g. killed for memory) returns an error object
  rows <- lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    if (is.data.frame(r)) return(r)
    data.frame(source = sources[i], dest = outputs[i], ok = FALSE, seconds = NA_real_,
               error = paste("worker failed:", conditionMessage(attr(r, 'condition') %||%
                                                                  simpleError(as.character(r)))),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (verbose) message(sum(out$ok), " of ", nrow(out), " conversions succeeded")
  out
}

# Output path for every source
.batch_outputs <- function(sources, dest, format) {
  if (length(dest) == length(sources) && length(sources) > 1L) return(as.character(dest))
  if (length(dest) != 1L) stop("'dest' must be a format, a directory or one path per source", call. = FALSE)
  stem <- tools::file_path_sans_ext(basename(sources))
  if (dir.exists(dest)) {
    if (is.null(format)) stop("Give 'format' when 'dest' is a directory", call. = FALSE)
    return(file.path(dest, paste0(stem, '.', tolower(format))))
  }
  if (length(sources) == 1L && nzchar(tools::file_ext(dest))) return(dest)
  file.path(dirname(sources), paste0(stem, '.', tolower(dest)))
}
