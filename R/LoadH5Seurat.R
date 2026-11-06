#' @include zzz.R
#' @include h5Seurat.R
#' @include GetObject.R
#' @include AssembleObject.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal helpers for multi-library detection
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Detect if h5Seurat file contains multi-library spatial data
#' @param file h5Seurat file object
#' @return List with is_multilibrary and library_col
#' @keywords internal
DetectMultiLibrary <- function(file) {
  result <- list(is_multilibrary = FALSE, library_col = NULL)

  if (!file$exists(name = 'images')) {
    return(result)
  }

  # Count images with centroids
  # Use safe nested path check to avoid HDF5 errors on missing intermediate groups
  img_names <- names(file[['images']])
  images_with_centroids <- vapply(img_names, function(img_name) {
    img_group <- file[['images']][[img_name]]
    # Check each path component separately to avoid HDF5 nested path errors
    img_group$exists('boundaries') &&
      img_group[['boundaries']]$exists('centroids') &&
      img_group[['boundaries/centroids']]$exists('cells')
  }, logical(1))

  if (sum(images_with_centroids) < 2) {
    return(result)
  }

  # Look for library ID column
  if (file$exists('meta.data')) {
    meta_group <- file[['meta.data']]
    lib_col_names <- c('sangerID', 'library_id', 'sample', 'batch')
    for (col_name in lib_col_names) {
      if (meta_group$exists(col_name)) {
        result$library_col <- col_name
        result$is_multilibrary <- TRUE
        break
      }
    }
  }

  result
}

#' Read library IDs from h5Seurat metadata
#' @param file h5Seurat file object
#' @param library_col Name of the library ID column
#' @return Character vector of library IDs
#' @keywords internal
ReadLibraryIds <- function(file, library_col) {
  col_obj <- file[['meta.data']][[library_col]]

  # Handle factor structure
  if (inherits(col_obj, 'H5Group') &&
      col_obj$exists('values') && col_obj$exists('levels')) {
    return(.DecodeH5SeuratFactor(col_obj))
  }

  as.character(col_obj$read())
}

# Decode an h5Seurat factor group (`values` + `levels`). h5Seurat stores
# 1-based codes (as.integer(factor) on write; the h5ad converters add 1);
# out-of-range codes are missing values.
.DecodeH5SeuratFactor <- function(col_obj) {
  values <- as.integer(col_obj[['values']]$read())
  levels <- as.character(col_obj[['levels']]$read())
  values[!is.na(values) & (values < 1L | values > length(levels))] <- NA_integer_
  levels[values]
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load a Seurat object from an h5Seurat file
#'
#' Load a previously saved Seurat object from an h5Seurat file. This function supports
#' flexible loading options, allowing you to load only the components you need (e.g.,
#' specific assays, reductions) to minimize memory usage on large datasets.
#'
#' @param file,x Name of an h5Seurat file path (character) or connected h5Seurat file to load
#' @param assays One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all assays
#'  \item A character vector with names of assays to load (e.g., \code{c("RNA", "ADT")})
#'  \item A character vector specifying which data layers to load for all assays:
#'    \code{c("counts", "data")} loads only the counts and data layers, skipping scale.data
#'  \item A named list for fine-grained control, e.g., \code{list(RNA = "data", ADT = c("data", "scale.data"))}
#' }
#' @param reductions One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all reductions (PCA, UMAP, etc.)
#'  \item A character vector with names of specific reductions (e.g., \code{c("pca", "umap")})
#'  \item \code{NA}: Load only global (assay-independent) reductions
#'  \item \code{FALSE}: Skip loading all reductions
#' }
#' \strong{Note}: Only reductions associated with a loaded assay or marked as global will be loaded.
#' @param graphs One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all graphs
#'  \item A character vector with specific graph names (e.g., \code{c("RNA_snn", "ADT_snn")})
#'  \item \code{FALSE}: Skip loading graphs
#' }
#' \strong{Note}: Only graphs associated with loaded assays will be available.
#' @param neighbors One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all neighbor information
#'  \item A character vector with neighbor names
#'  \item \code{FALSE}: Skip neighbors
#' }
#' @param images One of:
#' \itemize{
#'  \item \code{NULL} (default): Load all images (for spatial experiments)
#'  \item A character vector with image names
#'  \item \code{NA}: Load only global images
#'  \item \code{FALSE}: Skip images
#' }
#' @param meta.data Logical; if \code{TRUE} (default), load cell metadata
#' @param commands Logical; if \code{TRUE} (default), load command history.
#'   Commands are only loaded if their associated assays are loaded.
#' @param misc Logical; if \code{TRUE} (default when all assays loaded), load miscellaneous data
#' @param tools Logical; if \code{TRUE} (default when all assays loaded), load tool-specific information
#' @param verbose Logical; if \code{TRUE} (default), show progress messages
#' @param ... Arguments passed to other methods
#'
#' @return A \code{Seurat} object containing the requested components
#'
#' @details
#' The h5Seurat format is highly flexible for selective loading. This is particularly useful when:
#' \itemize{
#'   \item Working with very large datasets where loading everything would exceed memory
#'   \item You only need specific assays or reductions for downstream analysis
#'   \item You want to quickly inspect object structure without full data loading
#' }
#'
#' @section Seurat V5 Layer Support:
#' For Seurat V5 objects with multiple layers, you can selectively load layers per assay.
#' For example, use \code{assays = list(RNA = "data")} to load only the normalized expression layer,
#' skipping raw counts and scaled data.
#'
#' @seealso
#' \code{\link{writeH5Seurat}} to save a Seurat object to h5Seurat format
#' \code{\link{scConvert}} to convert to other formats
#'
#' @examples
#' \dontrun{
#' library(scConvert)
#'
#' # Load entire h5Seurat file
#' seurat_obj <- readH5Seurat("data.h5seurat")
#'
#' # Load only specific assays
#' seurat_obj <- readH5Seurat("data.h5seurat", assays = c("RNA", "ADT"))
#'
#' # Load only specific data layers (memory-efficient for large files)
#' seurat_obj <- readH5Seurat("data.h5seurat", assays = c("data"))  # Only normalized expression
#'
#' # Load specific assays with different layers
#' seurat_obj <- readH5Seurat(
#'   "data.h5seurat",
#'   assays = list(RNA = c("data", "scale.data"), ADT = "data")
#' )
#'
#' # Load without reductions (faster)
#' seurat_obj <- readH5Seurat("data.h5seurat", reductions = FALSE)
#'
#' # Load UMAP and PCA reductions only
#' seurat_obj <- readH5Seurat("data.h5seurat", reductions = c("umap", "pca"))
#'
#' # Load spatial data without graphs (for Visium experiments)
#' seurat_obj <- readH5Seurat("visium.h5seurat", images = TRUE, graphs = FALSE)
#' }
#'
#' @export
#'
readH5Seurat <- function(file, ...) {
  UseMethod(generic = 'readH5Seurat', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Loading options shared by every readH5Seurat() method, collected from the
# calling method's frame
.h5s_load_options <- function(env = parent.frame()) {
  mget(c('assays', 'reductions', 'graphs', 'neighbors', 'images', 'meta.data',
         'commands', 'misc', 'tools', 'verbose'), envir = env)
}

#' @rdname readH5Seurat
#' @method readH5Seurat character
#' @export
#'
readH5Seurat.character <- function(file, assays = NULL, reductions = NULL, graphs = NULL,
                                   neighbors = NULL, images = NULL, meta.data = TRUE,
                                   commands = TRUE, misc = is.null(x = assays),
                                   tools = is.null(x = assays), verbose = TRUE, ...) {
  hfile <- h5Seurat$new(filename = file, mode = 'r')
  # HDF5 1.12.1 can error while closing handles that are still referenced;
  # the data has been read by then
  on.exit(tryCatch(hfile$close_all(), error = function(e) NULL))
  do.call(readH5Seurat, c(list(file = hfile), .h5s_load_options(), list(...)))
}

#' @rdname readH5Seurat
#' @method readH5Seurat H5File
#' @export
#'
readH5Seurat.H5File <- function(file, assays = NULL, reductions = NULL, graphs = NULL,
                                neighbors = NULL, images = NULL, meta.data = TRUE,
                                commands = TRUE, misc = is.null(x = assays),
                                tools = is.null(x = assays), verbose = TRUE, ...) {
  do.call(readH5Seurat, c(list(file = as.h5Seurat(x = file)), .h5s_load_options(), list(...)))
}

#' @importFrom Seurat as.Seurat
#'
#' @rdname readH5Seurat
#' @method readH5Seurat h5Seurat
#' @export
#'
readH5Seurat.h5Seurat <- function(file, assays = NULL, reductions = NULL, graphs = NULL,
                                  neighbors = NULL, images = NULL, meta.data = TRUE,
                                  commands = TRUE, misc = is.null(x = assays),
                                  tools = is.null(x = assays), verbose = TRUE, ...) {
  opts <- .h5s_load_options()
  libs <- DetectMultiLibrary(file)
  if (!libs$is_multilibrary || !is.null(images)) {
    return(do.call(as.Seurat, c(list(x = file), opts, list(...))))
  }
  # Several spatial libraries with their own images: one object per library
  if (verbose) message("Detected multi-library spatial data; returning one Seurat object per library")
  lib_of_cell <- ReadLibraryIds(file, libs$library_col)
  cells <- Cells(x = file)
  ids <- unique(lib_of_cell)
  if (verbose) message("Found ", length(ids), " libraries: ", paste(ids, collapse = ", "))
  # images are added per library after subsetting, to keep cells consistent
  opts$images <- character(0)
  out <- list()
  for (lib in ids) {
    keep <- cells[!is.na(lib_of_cell) & lib_of_cell == lib]
    if (!length(keep)) next
    if (verbose) message("\nLoading library: ", lib)
    obj <- subset(do.call(as.Seurat, c(list(x = file), opts, list(...))), cells = keep)
    img <- tryCatch(AssembleImage(image = lib, file = file, verbose = FALSE), error = function(e) {
      if (verbose) message("  Could not load image ", lib, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(img)) {
      obj@images[[lib]] <- img
      if (verbose) message("  Added image ", lib)
    }
    if (verbose) message("  Loaded ", ncol(obj), " cells from library ", lib)
    out[[lib]] <- obj
  }
  if (verbose) message("\nReturning ", length(out), " Seurat objects")
  out
}

# Assay objects for the requested assays/slots. Assays with nothing on
# disk become empty placeholders.
.h5s_load_assays <- function(x, index, requested, verbose) {
  out <- list()
  for (assay in names(requested)) {
    slots <- requested[[assay]]
    if (!is.character(slots) || !length(slots)) {
      # everything the index lists (counts/data/scale.data, split and extra
      # v5 layers), else the three standard slots
      slots <- tryCatch(names(Filter(isTRUE, index[[assay]]$slots)), error = function(e) NULL)
      if (!length(slots)) slots <- c('counts', 'data', 'scale.data')
    }
    on_disk <- tryCatch({
      grp <- x[['assays']][[assay]]
      found <- names(grp)
      if ('layers' %in% found && inherits(grp[['layers']], 'H5Group')) {
        found <- unique(c(found, names(grp[['layers']])))
      }
      found
    }, error = function(e) character(0))
    if (!length(on_disk)) {
      if (verbose) message(sprintf("Assay '%s' has no data on disk; adding an empty assay", assay))
      out[[assay]] <- SeuratObject::CreateAssayObject(counts = matrix(numeric(0), nrow = 0, ncol = 0))
      next
    }
    use <- intersect(slots, on_disk)
    if (!length(use)) {
      use <- on_disk[1]
      if (verbose) {
        message(sprintf("None of the requested slots exist for assay '%s'; using '%s'", assay, use))
      }
    }
    out[[assay]] <- AssembleAssay(assay = assay, file = x, slots = use, verbose = verbose)
  }
  out
}

# Should images be loaded? From the spatial layout on, or for v4 SliceImage
# objects when an older Seurat is installed
.h5s_images_supported <- function(x) {
  if (packageVersion(pkg = 'Seurat') >= numeric_version(spatial.version)) return(TRUE)
  if (!x$exists(name = 'images')) return(FALSE)
  isTRUE(tryCatch(any(vapply(names(x[['images']]), function(img) {
    g <- x[['images']][[img]]
    g$attr_exists(attr_name = 's4class') && identical(h5attr(x = g, which = 's4class'), 'SliceImage')
  }, logical(1))), error = function(e) FALSE))
}

#' @importClassesFrom SeuratObject Seurat
#' @importFrom methods slot<-
#' @importFrom SeuratObject DefaultAssay Cells Idents<- Idents Project<- Project AddMetaData
#'
#' @aliases as.Seurat
#'
#' @rdname readH5Seurat
#' @method as.Seurat h5Seurat
#' @export
#'
as.Seurat.h5Seurat <- function(x, assays = NULL, reductions = NULL, graphs = NULL,
                               neighbors = NULL, images = NULL, meta.data = TRUE,
                               commands = TRUE, misc = TRUE, tools = TRUE,
                               verbose = TRUE, ...) {
  index <- x$index()
  load_everything <- is.null(assays) && is.null(reductions) && is.null(graphs)
  requested <- GetAssays(assays = assays, index = index)
  active <- DefaultAssay(object = index)
  if (!active %in% names(requested)) {
    active <- names(requested)[1]
    warning("The default assay was not requested; using ", active, " instead",
            call. = FALSE, immediate. = TRUE)
  }
  assay_objs <- .h5s_load_assays(x, index, requested, verbose)
  object <- new(
    Class = 'Seurat',
    assays = stats::setNames(list(assay_objs[[active]]), active),
    active.assay = active,
    meta.data = data.frame(row.names = Cells(x = x)),
    version = package_version(x$version())
  )
  for (a in setdiff(names(assay_objs), active)) object[[a]] <- assay_objs[[a]]

  all_global <- isTRUE(getOption('scConvert.dimreducs.allglobal', default = FALSE))
  for (red in GetDimReducs(reductions = reductions, index = index, assays = requested)) {
    if (verbose) message("Adding reduction ", red)
    dr <- AssembleDimReduc(reduction = red, file = x, verbose = verbose)
    if (all_global) slot(dr, 'global') <- TRUE
    object[[red]] <- dr
  }
  for (g in GetGraphs(graphs = graphs, index = index, assays = requested)) {
    if (verbose) message("Adding graph ", g)
    object[[g]] <- AssembleGraph(graph = g, file = x, verbose = verbose)
  }
  for (nb in GetNeighbors(neighbors = neighbors, index = index)) {
    if (verbose) message("Adding neighbors ", nb)
    object[[nb]] <- AssembleNeighbor(neighbor = nb, file = x, verbose = verbose)
  }

  # images that cannot be attached to the object are kept in misc
  orphan_images <- list()
  if (.h5s_images_supported(x)) {
    for (img in GetImages(images = images, index = index, assays = requested)) {
      if (verbose) message("Adding image ", img)
      obj_img <- tryCatch(AssembleImage(image = img, file = x, verbose = verbose), error = function(e) {
        if (verbose) message("Failed to load image ", img, ": ", conditionMessage(e))
        NULL
      })
      if (is.null(obj_img)) next
      attached <- tryCatch({ object[[img]] <- obj_img; TRUE }, error = function(e) FALSE)
      if (!attached) {
        if (verbose) message("Storing image ", img, " in misc slot")
        orphan_images[[paste0('image_', img)]] <- obj_img
      }
    }
  }

  if (commands) {
    if (verbose) message("Adding command information")
    cmds <- GetCommands(index = index, assays = requested)
    slot(object, 'commands') <- stats::setNames(
      lapply(cmds, AssembleSeuratCommand, file = x, verbose = verbose), cmds)
  }

  if (meta.data && x$exists('meta.data')) {
    if (verbose) message("Adding cell-level metadata")
    tryCatch({
      md <- as.data.frame(x = x[['meta.data']], row.names = Cells(x = x))
      md <- md[, !vapply(md, function(col) is.null(col) || is.environment(col), logical(1)),
               drop = FALSE]
      if (ncol(md)) object <- AddMetaData(object = object, metadata = md)
    }, error = function(e) {
      if (verbose) {
        warning("Could not load metadata: ", conditionMessage(e), call. = FALSE, immediate. = TRUE)
      }
    })
  }
  Idents(object = object) <- Idents(object = x)
  Project(object = object) <- Project(object = x)

  if (misc) {
    if (verbose) message("Adding miscellaneous information")
    # SafeH5GroupToList copes with 3-D+ arrays
    slot(object, 'misc') <- SafeH5GroupToList(h5obj = x[['misc']], recursive = TRUE)
  }
  if (length(orphan_images)) slot(object, 'misc')[names(orphan_images)] <- orphan_images
  if (tools) {
    if (verbose) message("Adding tool-specific results")
    slot(object, 'tools') <- SafeH5GroupToList(h5obj = x[['tools']], recursive = TRUE)
  }

  # graphs and command logs not tied to any assay
  if (load_everything && !is.null(index$no.assay)) {
    if (verbose) message("Adding data that was not associated with an assay")
    for (g in index$no.assay$graphs) {
      object[[g]] <- AssembleGraph(graph = g, file = x, verbose = verbose)
    }
    for (cmd in index$no.assay$commands) {
      object[[cmd]] <- AssembleSeuratCommand(cmd = cmd, file = x, verbose = verbose)
    }
  }
  object
}
