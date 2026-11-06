#' @include h5Seurat.R
#' @include loom.R
#' @importFrom hdf5r H5File h5attr
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Opening files
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Open a connection to an h5Seurat or loom file
#'
#' @param filename Path to the file
#' @param type File type (\code{"h5seurat"} or \code{"loom"}); guessed from
#'   the extension when \code{NULL}
#' @param mode \code{"r"} for read-only, \code{"r+"} for read-write
#' @param force When the file cannot be opened as \code{type} (it fails
#'   validation, or the type has no connection class, as for h5ad and h5mu),
#'   return a plain \code{\link[hdf5r]{H5File}} instead of an error
#'
#' @return An \code{\link{h5Seurat}} or \code{\link{loom}} connection (or an
#'   \code{H5File} with \code{force = TRUE})
#'
#' @export
#'
scConnect <- function(filename, type = NULL, mode = c('r', 'r+'), force = FALSE) {
  type <- type %||% FileType(file = filename)
  mode <- match.arg(mode)
  if (!file.exists(filename)) {
    stop("No ", type, " file at ", filename, call. = FALSE)
  }
  tryCatch(
    GetSCDisk(r6class = type)$new(filename = filename, mode = mode),
    error = function(err) {
      if (!isTRUE(force)) stop(conditionMessage(err), call. = FALSE)
      H5File$new(filename = filename, mode = mode)
    }
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SeuratObject accessors for h5Seurat connections
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' SeuratObject accessors for open h5Seurat files
#'
#' Read (and for the default assay and project name, change) object-level
#' information directly from an h5Seurat connection without loading it.
#'
#' @param x,object An \code{\link{h5Seurat}} connection, or for \code{Key} and
#'   \code{IsGlobal} an HDF5 group inside one
#' @param value New default assay or project name
#' @param reduction Reduction to read standard deviations from
#' @param ... Ignored
#'
#' @return \code{Cells}: cell names; \code{DefaultAssay}, \code{Project},
#'   \code{Key}: a string; \code{Idents}: a factor (or \code{NULL});
#'   \code{IsGlobal}: \code{TRUE}/\code{FALSE}; \code{Stdev}: a numeric
#'   vector; replacement functions return the connection invisibly
#'
#' @name h5Seurat-bindings
#' @rdname h5Seurat-bindings
#'
NULL

# Replace a scalar string attribute on the file root
.set_root_attr <- function(object, name, value) {
  if (object$attr_exists(attr_name = name)) object$attr_delete(attr_name = name)
  object$create_attr(attr_name = name, robj = value, dtype = GuessDType(x = value))
  invisible(object)
}

#' @importFrom SeuratObject Cells
#' @rdname h5Seurat-bindings
#' @method Cells h5Seurat
#' @export
#'
Cells.h5Seurat <- function(x, ...) {
  if (!x$exists(name = 'cell.names')) {
    stop("This h5Seurat file has no cell.names dataset", call. = FALSE)
  }
  ds <- x[['cell.names']]
  # some v5 writers stored cell names as an n x 1 table
  as.character(if (length(ds$dims) > 1L) ds[, 1] else ds[])
}

#' @importFrom SeuratObject DefaultAssay
#' @rdname h5Seurat-bindings
#' @method DefaultAssay h5Seurat
#' @export
#'
DefaultAssay.h5Seurat <- function(object, ...) {
  h5attr(x = object, which = 'active.assay')
}

#' @importFrom SeuratObject DefaultAssay<-
#' @rdname h5Seurat-bindings
#' @method DefaultAssay<- h5Seurat
#' @export
#'
"DefaultAssay<-.h5Seurat" <- function(object, ..., value) {
  if (!value %in% names(object[['assays']])) {
    stop("No assay named '", value, "' in this h5Seurat file", call. = FALSE)
  }
  .set_root_attr(object, 'active.assay', value)
}

#' @importFrom SeuratObject Idents
#' @rdname h5Seurat-bindings
#' @method Idents h5Seurat
#' @export
#'
Idents.h5Seurat <- function(object, ...) {
  if (!object$exists('active.ident')) return(NULL)
  ident <- object[['active.ident']]
  if (inherits(ident, "H5Group") && all(c("levels", "values") %in% names(ident))) {
    lv <- ident[['levels']][]
    return(factor(lv[as.integer(ident[['values']][])], levels = lv))
  }
  if (inherits(ident, "H5D")) return(factor(ident[]))
  NULL
}

#' @importFrom SeuratObject IsGlobal
#' @rdname h5Seurat-bindings
#' @method IsGlobal H5Group
#' @export
#'
IsGlobal.H5Group <- function(object, ...) {
  isTRUE(object$attr_exists(attr_name = 'global')) &&
    isTRUE(as.integer(h5attr(x = object, which = 'global')) == 1L)
}

#' @importFrom SeuratObject Key
#' @rdname h5Seurat-bindings
#' @method Key H5Group
#' @export
#'
Key.H5Group <- function(object, ...) {
  h5attr(x = object, which = 'key')
}

#' @importFrom SeuratObject Project
#' @rdname h5Seurat-bindings
#' @method Project h5Seurat
#' @export
#'
Project.h5Seurat <- function(object, ...) {
  h5attr(x = object, which = 'project')
}

#' @importFrom SeuratObject Project<-
#' @rdname h5Seurat-bindings
#' @method Project<- h5Seurat
#' @export
#'
"Project<-.h5Seurat" <- function(object, ..., value) {
  .set_root_attr(object, 'project', value)
}

#' @importFrom SeuratObject Stdev
#' @rdname h5Seurat-bindings
#' @method Stdev h5Seurat
#' @export
#'
Stdev.h5Seurat <- function(object, reduction = 'pca', ...) {
  path <- paste0('reductions/', reduction, '/stdev')
  if (Exists(x = object, name = path)) object[[path]][] else numeric(0)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Accessors for loom connections
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Accessors for open loom files
#'
#' @param object,x A \code{\link{loom}} connection
#' @param ... Ignored
#'
#' @return \code{DefaultAssay}: the Seurat assay recorded by
#'   \code{\link{writeLoom}} (or \code{NULL}); \code{dim}: genes and cells
#'
#' @name loom-bindings
#' @rdname loom-bindings
#'
NULL

#' @rdname loom-bindings
#' @method DefaultAssay loom
#' @export
#'
DefaultAssay.loom <- function(object, ...) {
  if (Exists(x = object, name = 'attrs/SEURAT_ASSAY')) {
    return(object[['attrs/SEURAT_ASSAY']][])
  }
  if (AttrExists(x = object, name = 'SEURAT_ASSAY')) {
    return(h5attr(x = object, which = 'SEURAT_ASSAY'))
  }
  warning("This loom file does not record a default assay", call. = FALSE, immediate. = TRUE)
  NULL
}

#' @rdname loom-bindings
#' @method dim loom
#' @export
#'
dim.loom <- function(x) {
  # hdf5r reports /matrix as cells x genes; loom's view is genes x cells
  rev(x[['matrix']]$dims)
}
