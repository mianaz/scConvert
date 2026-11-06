#' @include zzz.R
#' @include H5Utils.R
#' @importFrom R6 R6Class is.R6Class
#' @importFrom hdf5r H5File h5attr
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Base class for file connections
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base class for single-cell file connections
#'
#' \code{scdisk} extends \code{\link[hdf5r]{H5File}} with validation on
#' connect, chunk iteration and timestamps. It is abstract: open files
#' through a format class such as \code{\link{h5Seurat}} or \code{\link{loom}}
#' (usually via \code{\link{scConnect}}).
#'
#' @name scdisk-class
#' @aliases scdisk
#' @rdname scdisk-class
#' @docType class
#' @format An \code{\link[R6]{R6Class}} generator
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @export
#'
scdisk <- R6Class(
  'scdisk',
  inherit = H5File,
  lock_class = TRUE,
  cloneable = FALSE,
  public = list(
    #' @description Open or create a file and validate it
    #' @param filename Path to the file
    #' @param mode \code{"r"} / \code{"r+"} open an existing file (read-only /
    #'   read-write); \code{"w"} creates or truncates; \code{"w-"} and
    #'   \code{"x"} create and fail if the file exists; \code{"a"} opens or
    #'   creates read-write
    #' @param validate Run the format's validation after opening
    #' @param ... Passed to the validation routine
    initialize = function(filename = NULL,
                          mode = c('a', 'r', 'r+', 'w', 'w-', 'x'),
                          validate = TRUE, ...) {
      mode <- match.arg(mode)
      if (mode %in% c('r', 'r+') && !file.exists(filename)) {
        stop("File not found: ", filename, call. = FALSE)
      }
      super$initialize(filename = filename, mode = mode)
      private$validate(validate = validate, ...)
    },
    #' @description Release every HDF5 handle when the object is collected
    finalizer = function() {
      self$close_all(close_self = TRUE)
    },
    #' @description Start/end indices for iterating over a dataset in blocks
    #' @param dataset Path of a dataset in the file
    #' @param MARGIN Axis to iterate over (see \code{GetMargin})
    #' @param csize Block size; defaults to the dataset's chunk size on that
    #'   axis
    #' @return A two-column matrix (\code{start}, \code{end}), one row per block
    chunk.points = function(dataset,
                            MARGIN = getOption('scConvert.chunking.MARGIN', default = 'largest'),
                            csize = NULL) {
      if (!private$is.data(name = dataset, type = 'H5D')) {
        stop(sprintf("'%s' is not a dataset in this %s file", dataset, class(self)[1]),
             call. = FALSE)
      }
      ds <- self[[dataset]]
      axis <- GetMargin(dims = ds$dims, MARGIN = MARGIN)
      ChunkPoints(dsize = ds$dims[axis], csize = csize %||% ds$chunk_dims[axis])
    },
    #' @description Record the current time in an attribute
    #' @param name Object to stamp; \code{NULL} stamps the file itself
    #' @param attr Attribute name
    #' @param tz,format Time zone and \code{\link[base]{strptime}} format
    #' @return Invisibly, the object
    timestamp = function(name = NULL, attr = 'ts', tz = 'UTC',
                         format = TSFormats(type = 'R')) {
      if (isTRUE(Writeable(x = self, error = FALSE))) {
        target <- private$member(name)
        now <- format(Sys.time(), format = format, tz = tz)
        if (target$attr_exists(attr_name = attr)) target$attr_delete(attr_name = attr)
        target$create_attr(attr_name = attr, robj = now, dtype = GuessDType(x = now))
      }
      invisible(self)
    },
    #' @description Read a timestamp written by \code{timestamp()}
    #' @param name Object to read from; \code{NULL} reads the file itself
    #' @param attr Attribute name
    #' @param locale Show the time in the local time zone
    #' @param tz,format Time zone and format the stamp was written with
    #' @return The timestamp as a character string
    last.modified = function(name = NULL, attr = 'ts', locale = TRUE, tz = 'UTC',
                             format = TSFormats(type = 'R')) {
      target <- private$member(name)
      if (!AttrExists(x = target, name = attr)) {
        stop("No timestamp attribute '", attr, "' found", call. = FALSE)
      }
      FormatTime(time = h5attr(x = target, which = attr), locale = locale,
                 tz = tz, format = format)
    }
  ),
  private = list(
    # The file itself (name = NULL) or one of its members
    member = function(name) {
      if (is.null(name)) return(self)
      if (!Exists(x = self, name = name)) {
        stop("No member named '", name, "' in this file", call. = FALSE)
      }
      self[[name]]
    },
    # Canned error messages shared by the format classes
    errors = function(type = c('mode', 'ambiguous')) {
      switch(
        match.arg(type),
        mode = sprintf("This %s file is open read-only and cannot be modified",
                       class(self)[1]),
        ambiguous = "More than one dataset matches; please give a more specific name"
      )
    },
    # Does `name` exist and is it a dataset (H5D) / group (H5Group)?
    is.data = function(name, type = c('H5D', 'H5Group')) {
      type <- match.arg(type)
      isTRUE(self$exists(name = name)) && inherits(self[[name]], type)
    },
    # Format classes override this; it receives at least `validate`
    validate = function(...) {
      if (identical(class(self)[1], 'scdisk')) {
        stop("scdisk is abstract; open a file through a format class", call. = FALSE)
      }
      warning(class(self)[1], " defines no validation routine",
              call. = FALSE, immediate. = TRUE)
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Registry of format classes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Register and look up \code{scdisk} format classes
#'
#' Format classes (subclasses of \code{\link{scdisk}}) are registered under
#' their class name so that code which only knows a file type, such as
#' \code{\link{scConnect}}, can find the right generator. Packages that add a
#' format register it from their \code{.onLoad}.
#'
#' @param r6class For \code{GetSCDisk}: a class name (case-insensitive), or
#'   \code{NULL} to list registered names. For \code{RegisterSCDisk}: an
#'   \code{\link[R6]{R6Class}} generator inheriting from \code{scdisk}
#'
#' @return \code{GetSCDisk}: the generator for \code{r6class}, or the names
#'   of all registered classes; \code{RegisterSCDisk}: invisibly \code{NULL}
#'
#' @name RegisterSCDisk
#' @rdname RegisterSCDisk
#' @keywords internal
#' @export
#'
GetSCDisk <- function(r6class = NULL) {
  known <- names(scdisk.types)
  if (is.null(r6class)) return(known)
  hit <- known[tolower(known) == tolower(r6class)]
  if (!length(hit)) stop("Unknown file type: ", r6class, call. = FALSE)
  scdisk.types[[hit[1]]]
}

#' @rdname RegisterSCDisk
#' @export
#'
RegisterSCDisk <- function(r6class) {
  if (!isTRUE(IsSCDisk(r6class = r6class))) {
    warning(r6class$classname, " is not an scdisk subclass; not registered",
            call. = FALSE, immediate. = TRUE)
    return(invisible(NULL))
  }
  owner <- environmentName(r6class$parent_env)
  if (r6class$classname %in% names(scdisk.types) &&
      owner != environmentName(scdisk$parent_env)) {
    stop("A format class named ", r6class$classname,
         " is already registered (by ", owner, ")", call. = FALSE)
  }
  assign(r6class$classname, r6class, envir = scdisk.types)
  invisible(NULL)
}

# TRUE when `r6class` (a generator or its name) is scdisk or inherits from it
IsSCDisk <- function(r6class) {
  if (is.character(r6class)) r6class <- get(r6class, mode = "environment")
  while (is.R6Class(r6class)) {
    if (identical(r6class, scdisk)) return(TRUE)
    r6class <- r6class$get_inherit()
  }
  FALSE
}

# strptime() formats for timestamps: the default one, and loom's (which
# carries fractional seconds)
TSFormats <- function(type = c('R', 'loom')) {
  if (match.arg(type) == 'loom') '%Y%m%dT%H%M%OS6Z' else '%Y%m%dT%H%M%SZ'
}

# Parse a stored timestamp and render it, in the local time zone by default
FormatTime <- function(time, locale = TRUE, tz = 'UTC', format = TSFormats(type = 'R')) {
  parsed <- as.POSIXct(strptime(time, format = format, tz = tz))
  format(parsed, tz = if (isTRUE(locale)) Sys.timezone() else tz, usetz = TRUE)
}
