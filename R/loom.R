#' @include scdisk.R
#' @importFrom R6 R6Class
#' @importFrom hdf5r h5attr H5S
#' @importFrom utils setTxtProgressBar
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# loom file connection
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Connection to a loom file
#'
#' An \code{\link{scdisk}} subclass for \href{http://linnarssonlab.org/loompy/format/}{loom}
#' files. A loom file holds a main \code{/matrix} (genes x cells as seen by
#' loompy), extra matrices under \code{/layers}, per-gene and per-cell
#' vectors under \code{/row_attrs} and \code{/col_attrs}, sparse graphs under
#' \code{/row_graphs} and \code{/col_graphs}, and global values under
#' \code{/attrs}. New files are created with the loom 3.0.0 layout.
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @aliases loom
#' @format An \code{\link[R6]{R6Class}} generator
#' @seealso \code{\link{scdisk}}, \code{\link{readLoom}}, \code{\link{writeLoom}}
#'
#' @export
#'
loom <- R6Class(
  classname = 'loom',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    #' @description Add a global, per-gene (\code{"row"}) or per-cell
    #'   (\code{"col"}) attribute
    #' @param x A vector (or, for global attributes, a matrix)
    #' @param name Attribute name
    #' @param type \code{"global"}, \code{"row"} or \code{"col"}
    #' @return Invisibly, the object
    add_attribute = function(x, name, type = c('global', 'row', 'col')) {
      if (is.list(x)) stop("Lists cannot be stored as loom attributes", call. = FALSE)
      type <- match.arg(type)
      name <- basename(name)
      if (type == 'global') {
        .loom_write_global(lfile = self, x = x, name = name, stype = private$stype())
      } else {
        .loom_write_vector(lfile = self, x = x, name = name, axis = type,
                           stype = private$stype())
      }
      invisible(self)
    },
    #' @description Add a sparse graph over genes (\code{"row"}) or cells
    #'   (\code{"col"})
    #' @param x A \code{dgCMatrix}
    #' @param name Graph name
    #' @param type \code{"col"} or \code{"row"}
    #' @param verbose Unused, kept for compatibility
    #' @return Invisibly, the object
    add_graph = function(x, name, type = c('col', 'row'), verbose = TRUE) {
      .loom_write_graph(lfile = self, x = x, name = basename(name), axis = match.arg(type))
      invisible(self)
    },
    #' @description Add a matrix under \code{/layers}
    #' @param x A matrix with the same shape as the object (genes x cells), or
    #'   its transpose when \code{transpose = FALSE}
    #' @param name Layer name
    #' @param transpose Transpose \code{x} before writing
    #' @param verbose Show progress
    #' @return Invisibly, the object
    add_layer = function(x, name, transpose = TRUE, verbose = TRUE) {
      expected <- if (isTRUE(transpose)) dim(self) else rev(dim(self))
      if (!identical(dim(x), expected)) {
        stop(sprintf("A layer must be %d x %d", expected[1], expected[2]), call. = FALSE)
      }
      parent <- dirname(name)
      if (parent %in% c('.', '')) {
        name <- paste0('/layers/', name)
      } else if (!parent %in% c('layers', '/layers')) {
        stop("Layers must sit directly under /layers", call. = FALSE)
      }
      if (isTRUE(verbose)) message("Adding layer ", basename(name))
      .loom_write_matrix(lfile = self, x = x, name = name, transpose = transpose,
                         verbose = verbose)
      invisible(self)
    },
    #' @description Loom specification version of the file
    #' @return A \code{\link[base]{numeric_version}}
    version = function() {
      key <- 'LOOM_SPEC_VERSION'
      if (AttrExists(x = self, name = key)) {
        return(numeric_version(h5attr(x = self, which = key)))
      }
      if (Exists(x = self, name = paste0('attrs/', key))) {
        return(numeric_version(self[[paste0('attrs/', key)]][]))
      }
      guess <- if (Exists(x = self, name = 'attrs')) '3.0.0' else '0.1.0'
      warning("No loom version recorded in this file; assuming ", guess,
              call. = FALSE, immediate. = TRUE)
      numeric_version(guess)
    },
    #' @description Record the modification time in \code{last_modified}
    #' @param name Object to stamp; \code{NULL} stamps the file itself
    #' @return Invisibly, the object
    timestamp = function(name = NULL) {
      super$timestamp(name = name, attr = 'last_modified', tz = 'UTC',
                      format = TSFormats(type = 'loom'))
      invisible(self)
    },
    #' @description Modification time recorded in \code{last_modified}
    #' @param name Object to read from; \code{NULL} reads the file itself
    #' @param locale Show the time in the local time zone
    #' @return The timestamp, or \code{NA} when there is none
    last.modified = function(name = NULL, locale = FALSE) {
      target <- private$member(name)
      if (!AttrExists(x = target, name = 'last_modified')) {
        warning("No timestamp found", call. = FALSE, immediate. = TRUE)
        return(NA_character_)
      }
      # loompy writes fractional seconds; keep whole seconds and the UTC marker
      stamp <- sub('\\..*$', '', h5attr(x = target, which = 'last_modified'))
      if (!endsWith(stamp, 'Z')) stamp <- paste0(stamp, 'Z')
      FormatTime(time = stamp, locale = locale)
    }
  ),
  private = list(
    # Lay out an empty loom 3.0.0 file
    create = function() {
      if (identical(self$mode, 'r')) stop(private$errors(type = 'mode'), call. = FALSE)
      for (g in c('layers', 'attrs', 'row_attrs', 'col_attrs', 'row_graphs', 'col_graphs')) {
        self$create_group(name = g)
      }
      suppressWarnings(self$add_attribute(x = '3.0.0', name = 'LOOM_SPEC_VERSION'))
      invisible(NULL)
    },
    # Strings are fixed 7-byte ASCII in loom 0.1 files, UTF-8 from 3.0 on
    stype = function() {
      v <- tryCatch(ClosestVersion(query = self$version(), targets = c('0.1.0', '3.0.0')),
                    error = function(e) NA_character_)
      if (identical(v, '0.1.0')) return('ascii7')
      if (!identical(v, '3.0.0')) {
        warning("Unrecognised loom version, writing strings as UTF-8",
                call. = FALSE, immediate. = TRUE)
      }
      'utf8'
    },
    validate = function(validate = TRUE) {
      if (self$mode %in% modes$new) private$create()
      invisible(NULL)
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Writers used by the loom class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check_loom <- function(lfile) {
  if (!inherits(lfile, 'loom')) stop("'lfile' must be a loom connection", call. = FALSE)
}

# Global value under /attrs; matrices are stored transposed (loompy layout)
.loom_write_global <- function(lfile, x, name, stype, transpose = TRUE) {
  .check_loom(lfile)
  if (is.matrix(x)) {
    if (isTRUE(transpose)) x <- t(x)
    dtype <- GuessDType(x = x[1, 1], stype = stype)
  } else {
    if (is.list(x)) stop("Lists cannot be stored as loom attributes", call. = FALSE)
    dtype <- GuessDType(x = x, stype = stype)
  }
  lfile$create_dataset(name = paste0('/attrs/', name), robj = x, dtype = dtype)
  invisible(NULL)
}

# Per-gene (axis = "row") or per-cell (axis = "col") vector
.loom_write_vector <- function(lfile, x, name, axis, stype, verbose = TRUE) {
  .check_loom(lfile)
  if (is.list(x)) stop("Lists cannot be stored as loom attributes", call. = FALSE)
  n <- if (axis == 'row') nrow(lfile) else ncol(lfile)
  if (length(x) != n) {
    stop(sprintf("A %s attribute needs %d values, got %d", axis, n, length(x)),
         call. = FALSE)
  }
  if (isTRUE(verbose)) message("Adding ", axis, " attribute ", name)
  lfile$create_dataset(name = paste0(axis, '_attrs/', name), robj = x,
                       dtype = GuessDType(x = x, stype = stype))
  invisible(NULL)
}

# Sparse graph as loom's coordinate triplets: a (row), b (col), w (weight),
# all 0-based
.loom_write_graph <- function(lfile, x, name, axis) {
  .check_loom(lfile)
  if (!inherits(x, 'dgCMatrix')) x <- methods::as(x, 'dgCMatrix')
  g <- lfile$create_group(name = paste0(axis, '_graphs/', name))
  g$create_dataset(name = 'a', robj = x@i, dtype = GuessDType(x = 1L))
  g$create_dataset(name = 'b', robj = PointerToIndex(p = x@p) - 1, dtype = GuessDType(x = 1L))
  g$create_dataset(name = 'w', robj = x@x, dtype = GuessDType(x = 1))
  invisible(NULL)
}

# Dense dataset from a dense or sparse matrix. Sparse input is written in
# blocks along the larger axis so it is never fully densified.
.loom_write_matrix <- function(lfile, x, name, transpose = TRUE, verbose = TRUE) {
  .check_loom(lfile)
  if (!inherits(x, 'sparseMatrix')) {
    x <- as.matrix(x)
    if (isTRUE(transpose)) x <- t(x)
    lfile$create_dataset(name = name, robj = x, dtype = GuessDType(x = x[1, 1]))
    return(invisible(NULL))
  }
  if (isTRUE(transpose)) x <- Matrix::t(x)
  lfile$create_dataset(name = name, dtype = GuessDType(x = x[1, 1]),
                       space = H5S$new(dims = dim(x)))
  ds <- lfile[[name]]
  axis <- GetMargin(dims = ds$dims)
  blocks <- ChunkPoints(dsize = ds$dims[axis], csize = ds$chunk_dims[axis])
  index <- list(NULL, NULL)
  index[[3L - axis]] <- seq_len(ds$dims[3L - axis])
  bar <- if (isTRUE(verbose)) PB() else NULL
  for (b in seq_len(nrow(blocks))) {
    index[[axis]] <- seq.int(blocks[b, 'start'], blocks[b, 'end'])
    ds$write(args = index, value = as.matrix(x[index[[1]], index[[2]]]))
    if (!is.null(bar)) setTxtProgressBar(pb = bar, value = b / nrow(blocks))
  }
  if (!is.null(bar)) close(bar)
  invisible(NULL)
}
