#' @include scdisk.R
#' @importFrom R6 R6Class
#' @importFrom hdf5r h5attr
#' @importFrom utils packageVersion
#' @importFrom cli symbol
#' @importFrom crayon red green yellow col_align
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# h5Seurat file connection
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Connection to an h5Seurat file
#'
#' An \code{\link{scdisk}} subclass for h5Seurat files, the HDF5 layout used
#' to store \code{Seurat} objects on disk. Opening an existing file checks
#' its structure and builds an index of what it contains (see
#' \code{$index()}); creating a file lays out the empty groups.
#'
#' @name h5Seurat-class
#' @aliases h5Seurat
#' @rdname h5Seurat-class
#' @docType class
#' @format An \code{\link[R6]{R6Class}} generator
#' @seealso \code{\link{scdisk}}, \code{\link{readH5Seurat}},
#'   \code{\link{writeH5Seurat}}
#'
#' @export
#'
h5Seurat <- R6Class(
  'h5Seurat',
  inherit = scdisk,
  lock_class = TRUE,
  cloneable = FALSE,
  public = list(
    #' @description Summary of the assays, layers, reductions, graphs,
    #'   images, commands and neighbors stored in the file
    #' @return An \code{h5SI} object (a nested list)
    index = function() {
      if (!length(private$index.internal)) {
        private$build.index(version = private$file.version())
      }
      private$index.internal
    },
    #' @description Set the h5Seurat version recorded in the file
    #' @param version A version string such as \code{"5.0.0"}
    #' @return Invisibly, the object
    set.version = function(version) {
      version <- as.character(version)
      if (!grepl(version.regex, version)) {
        stop("Not a valid version: ", version, call. = FALSE)
      }
      ClosestVersion(query = version, targets = private$versions)  # errors when unsupported
      if (self$attr_exists(attr_name = 'version')) self$attr_delete(attr_name = 'version')
      self$create_attr(attr_name = 'version', robj = version, dtype = GuessDType(x = version))
      if (!self$mode %in% modes$new) private$validate()
      invisible(self)
    },
    #' @description Version recorded in the file
    #' @return A character string
    version = function() {
      h5attr(x = self, which = 'version')
    },
    #' @description Seurat object generation the file was written from,
    #'   judged from its layout and, failing that, its version attribute
    #' @return \code{"V3"}, \code{"V4"} or \code{"V5"}
    detect.version = function() {
      first <- private$first.assay()
      if (!is.null(first)) {
        if (private$has.sparse.layers(first)) return("V5")
        if (first$exists("counts") && !inherits(first[["counts"]], "H5Group")) return("V4")
      }
      recorded <- tryCatch(self$version(), error = function(e) NULL)
      if (is.null(recorded)) return("V4")
      v <- numeric_version(recorded)
      if (v >= "5.0.0") "V5" else if (v >= "4.0.0") "V4" else "V3"
    },
    #' @description Was the file written from a Seurat v5 object?
    #' @return \code{TRUE} or \code{FALSE}
    is.v5 = function() {
      identical(self$detect.version(), "V5")
    }
  ),
  private = list(
    index.internal = list(),
    # Layout versions this class understands (each adds to the previous)
    versions = c('3.1.2', spatial.version, v5.version, '5.2.1'),

    file.version = function() {
      ClosestVersion(query = self$version(), targets = private$versions)
    },
    first.assay = function() {
      if (!self$exists('assays')) return(NULL)
      nm <- names(self[['assays']])
      if (length(nm)) self[['assays']][[nm[1]]] else NULL
    },
    # v5 assays keep each layer as a CSC group under <assay>/layers
    has.sparse.layers = function(assay) {
      if (!assay$exists("layers") || !inherits(assay[["layers"]], "H5Group")) return(FALSE)
      lyr <- assay[["layers"]]$names
      if (!length(lyr)) return(FALSE)
      first <- assay[["layers"]][[lyr[1]]]
      inherits(first, "H5Group") && all(c("data", "indices", "indptr") %in% names(first))
    },
    # Optional attribute as a string, NA when absent
    attr.or.na = function(obj, name) {
      if (obj$attr_exists(attr_name = name)) h5attr(x = obj, which = name) else NA_character_
    },
    missing.assay = function(what, assay) {
      warning(sprintf("%s refers to assay '%s', which is not in this h5Seurat file",
                      what, assay), call. = FALSE, immediate. = TRUE)
    },

    # Which slots/layers an assay group holds
    assay.entry = function(name, version) {
      grp <- self[['assays']][[name]]
      std <- c('counts', 'data', 'scale.data')
      if (version >= numeric_version(v5.version) && grp$attr_exists(attr_name = 's4class')) {
        layers <- if (grp$exists("layers") && inherits(grp[["layers"]], "H5Group")) {
          grp[["layers"]]$names
        } else {
          setdiff(grp$names, c("features", "scaled.features", "variable.features",
                               "meta.features", "misc", "layers"))
        }
        present <- stats::setNames(rep(TRUE, length(layers)), layers)
        present[setdiff(std, layers)] <- FALSE
        # scale.data is only usable together with the list of scaled features
        present[['scale.data']] <- present[['scale.data']] &&
          (grp$exists(name = 'scaled.features') ||
             (grp$exists(name = 'scale.data') &&
                grp[['scale.data']]$attr_exists(attr_name = 'scaled.features')))
        return(list(slots = present,
                    v5class = grp$attr_open(attr_name = 's4class')$read(),
                    layers = layers))
      }
      present <- stats::setNames(std %in% names(grp), std)
      present[['scale.data']] <- present[['scale.data']] && grp$exists(name = 'scaled.features')
      list(slots = present)
    },

    # v4-era files keep Visium images as SliceImage groups
    has.slice.images = function() {
      if (!self$exists(name = 'images')) return(FALSE)
      isTRUE(tryCatch(any(vapply(names(self[['images']]), function(img) {
        g <- self[['images']][[img]]
        g$attr_exists(attr_name = 's4class') && identical(h5attr(x = g, which = 's4class'), 'SliceImage')
      }, logical(1))), error = function(e) FALSE))
    },

    build.index = function(version) {
      version <- numeric_version(match.arg(version, choices = private$versions))
      assays <- names(self[['assays']])
      idx <- stats::setNames(lapply(assays, private$assay.entry, version = version), assays)
      known <- function(a) !is.na(a) && a %in% assays

      # Dimensional reductions: which components each one stores
      parts <- c('cell.embeddings', 'feature.loadings', 'feature.loadings.projected', 'jackstraw')
      for (red in names(self[['reductions']])) {
        grp <- self[['reductions']][[red]]
        assay <- h5attr(x = grp, which = 'active.assay')
        if (!known(assay)) {
          private$missing.assay(paste("Reduction", red), assay)
          next
        }
        has <- stats::setNames(parts %in% names(grp), parts)
        # loadings are only usable with the matching feature names
        has[['feature.loadings']] <- has[['feature.loadings']] && grp$exists(name = 'features')
        has[['feature.loadings.projected']] <- has[['feature.loadings.projected']] &&
          grp$exists(name = 'projected.features')
        idx[[assay]][['reductions']][[red]] <- has
        if (IsGlobal(object = grp)) idx$global$reductions <- c(idx$global$reductions, red)
      }

      # Graphs, grouped by the assay they were built from
      for (graph in names(self[['graphs']])) {
        assay <- private$attr.or.na(self[['graphs']][[graph]], 'assay.used')
        if (is.na(assay)) {
          idx$no.assay$graphs <- c(idx$no.assay$graphs, graph)
        } else if (known(assay)) {
          idx[[assay]]$graphs <- c(idx[[assay]]$graphs, graph)
        } else {
          private$missing.assay(paste("Graph", graph), assay)
        }
      }

      # Spatial images
      if (version >= numeric_version(spatial.version) || private$has.slice.images()) {
        for (img in names(self[['images']])) {
          grp <- self[['images']][[img]]
          if (!grp$attr_exists(attr_name = 'assay') || !grp$attr_exists(attr_name = 's4class')) next
          assay <- h5attr(x = grp, which = 'assay')
          if (known(assay)) {
            idx[[assay]]$images <- c(idx[[assay]]$images, img)
          } else {
            private$missing.assay(paste("Image", img), assay)
            idx$no.assay$images <- c(idx$no.assay$images, img)
          }
          if (IsGlobal(object = grp)) idx$global$images <- c(idx$global$images, img)
        }
      }

      # Command logs
      for (cmd in names(self[['commands']])) {
        assay <- private$attr.or.na(self[['commands']][[cmd]], 'assay.used')
        if (known(assay)) {
          idx[[assay]]$commands <- c(idx[[assay]]$commands, cmd)
        } else if (is.na(assay)) {
          idx$no.assay$commands <- c(idx$no.assay$commands, cmd)
        } else {
          private$missing.assay(paste("Command", cmd), assay)
        }
      }

      if (Exists(x = self, name = 'neighbors')) {
        idx$global$neighbors <- names(self[['neighbors']])
      }
      private$index.internal <- structure(idx, class = c('h5SI', 'list'),
                                          active.assay = DefaultAssay(object = self))
      invisible(NULL)
    },

    # Empty layout for a new file
    create = function(version, verbose = TRUE) {
      if (identical(self$mode, 'r')) stop(private$errors(type = 'mode'), call. = FALSE)
      version <- ClosestVersion(query = version, targets = private$versions)
      self$set.version(version = version)
      v <- numeric_version(version)
      for (g in c('assays', 'commands', 'neighbors', 'graphs', 'misc', 'reductions', 'tools')) {
        if (!private$is.data(name = g, type = 'H5Group')) self$create_group(name = g)
      }
      defaults <- c(active.assay = '', project = 'SeuratProject')
      for (a in names(defaults)) {
        if (!self$attr_exists(attr_name = a)) {
          self$create_attr(attr_name = a, robj = defaults[[a]], dtype = GuessDType(x = defaults[[a]]))
        }
      }
      if (v >= numeric_version(spatial.version)) self$create_group(name = 'images')
      invisible(self)
    },

    validate = function(verbose = TRUE, ...) {
      if (self$mode %in% modes$new) {
        version <- if (self$attr_exists(attr_name = 'version')) {
          h5attr(x = self, which = 'version')
        } else {
          packageVersion(pkg = 'Seurat')
        }
        private$create(version = version, verbose = verbose)
        return(invisible(NULL))
      }
      if (verbose) message("Validating h5Seurat file")
      if (!self$attr_exists(attr_name = 'version')) {
        stop("Not an h5Seurat file: the 'version' attribute is missing", call. = FALSE)
      }
      version <- private$file.version()
      private$check.core()
      if (numeric_version(version) >= numeric_version(v5.version)) private$check.v5()
      private$build.index(version = version)
      invisible(NULL)
    },

    # Requirements shared by every supported layout
    check.core = function() {
      for (a in c('project', 'active.assay', 'version')) {
        if (!self$attr_exists(attr_name = a)) {
          stop("h5Seurat file is missing the '", a, "' attribute", call. = FALSE)
        }
      }
      if (!private$is.data(name = 'cell.names')) {
        stop("h5Seurat file has no cell.names dataset", call. = FALSE)
      }
      n_cells <- self[['cell.names']]$dims
      if (length(n_cells) != 1L) stop("cell.names must be one-dimensional", call. = FALSE)
      if (private$is.data(name = 'meta.data')) {
        md <- self[['meta.data']]
        if (length(md$dims) != 1L) stop("meta.data must be one-dimensional", call. = FALSE)
        if (md$dims != n_cells) {
          stop("meta.data and cell.names have different numbers of cells", call. = FALSE)
        }
        if (!inherits(md$get_type(), 'H5T_COMPOUND')) {
          stop("A meta.data dataset must be a compound (data frame) type", call. = FALSE)
        }
      } else if (!private$is.data(name = 'meta.data', type = 'H5Group')) {
        stop("h5Seurat file has no meta.data", call. = FALSE)
      }
      if (!private$is.data(name = 'assays', type = 'H5Group')) {
        stop("h5Seurat file has no assays group", call. = FALSE)
      }
      if (!DefaultAssay(object = self) %in% names(self[['assays']])) {
        stop("The default assay is not stored in this file", call. = FALSE)
      }
      for (a in names(self[['assays']])) {
        if (!private$is.data(name = paste0('assays/', a), type = 'H5Group')) {
          stop("Assay '", a, "' is not stored as an HDF5 group", call. = FALSE)
        }
      }
      invisible(NULL)
    },
    check.v5 = function() {
      for (g in c('assays', 'reductions', 'graphs', 'misc', 'tools')) {
        if (!private$is.data(name = g, type = 'H5Group')) {
          stop("h5Seurat v5 file is missing the '", g, "' group", call. = FALSE)
        }
      }
      invisible(NULL)
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# h5SI: index of an h5Seurat file
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Index of an h5Seurat file
#'
#' \code{h5Seurat$index()} returns an \code{h5SI} object: one entry per assay
#' (stored slots/layers, reductions, graphs, images, commands) plus
#' \code{global} and \code{no.assay} entries.
#'
#' @param object,x An \code{h5SI} object
#' @param ... Ignored
#'
#' @return \code{DefaultAssay}: the default assay of the file;
#'   \code{print}: invisibly \code{x}
#'
#' @name h5SI
#' @rdname h5SI
#' @keywords internal
#'
NULL

#' @importFrom SeuratObject DefaultAssay
#' @rdname h5SI
#' @method DefaultAssay h5SI
#' @export
#'
DefaultAssay.h5SI <- function(object, ...) {
  attr(object, 'active.assay')
}

#' @rdname h5SI
#' @method print h5SI
#' @export
#'
print.h5SI <- function(x, ...) {
  mark <- function(flags) ifelse(flags, green(symbol$tick), red(symbol$cross))
  default <- DefaultAssay(object = x)
  assays <- setdiff(names(x), c('global', 'no.assay'))
  assays <- c(intersect(default, assays), setdiff(assays, default))
  part_labels <- c(cell.embeddings = 'Embeddings', feature.loadings = 'Loadings',
                   feature.loadings.projected = 'Projected', jackstraw = 'JackStraw')
  for (a in assays) {
    entry <- x[[a]]
    cat("Data for assay ", a,
        if (identical(a, default)) paste0(yellow(symbol$star), " (default assay)"), "\n", sep = "")
    slots <- c('counts', 'data', 'scale.data')
    width <- max(nchar(slots)) + 1L
    cat(col_align(slots, width = width, align = 'center'), "\n", sep = "")
    cat(col_align(mark(entry$slots[slots]), width = width, align = 'center'), "\n", sep = "")
    if (length(entry$reductions)) {
      cat("Dimensional reductions:\n")
      labels <- paste0(" ", names(entry$reductions), ": ")
      lw <- max(nchar(labels))
      cw <- max(nchar(part_labels)) + 1L
      cat(strrep(" ", lw), col_align(unname(part_labels), width = cw, align = 'center'),
          "\n", sep = "")
      for (i in seq_along(labels)) {
        has <- entry$reductions[[i]][names(part_labels)]
        cat(col_align(labels[i], width = lw), col_align(mark(has), width = cw, align = 'center'),
            "\n", sep = "")
      }
    }
    for (what in c('graphs', 'images')) {
      if (length(entry[[what]])) {
        cat(tools::toTitleCase(what), ":\n", sep = "")
        cat(paste0(" ", symbol$line, " ", entry[[what]]), sep = "\n")
      }
    }
  }
  invisible(x)
}
