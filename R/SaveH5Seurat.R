#' @include WriteH5Group.R
#'
NULL


#' Save a Seurat object to an h5Seurat file
#'
#' Save a Seurat object to the efficient HDF5-based h5Seurat format. This format
#' is optimized for large, complex single-cell datasets including multi-modal and
#' spatial data. h5Seurat files can be rapidly converted to other formats like h5ad
#' (AnnData) or h5mu (MuData) for interoperability with Python tools.
#'
#' @param object,x An object (typically a Seurat object)
#' @param filename Name of file to save the object to. If not provided, defaults to
#'   \code{<ProjectName>.h5seurat}. The .h5seurat extension is added automatically
#'   if not present.
#' @param overwrite Logical; if \code{TRUE}, overwrite an existing file with the
#'   same name. Default is \code{FALSE}.
#' @param verbose Show progress updates during save. Default is \code{TRUE}.
#' @param gzip Integer gzip compression level (0-9), or NULL to use the package default.
#' @param ... Arguments passed to other methods
#'
#' @return \code{writeH5Seurat}: Invisibly returns the filename of the saved file
#'
#' @details
#' The h5Seurat format stores:
#' \itemize{
#'   \item All assays with their layers (counts, data, scale.data, etc.) for Seurat V5
#'   \item Dimensional reductions (PCA, UMAP, etc.)
#'   \item Nearest-neighbor graphs and similarity graphs
#'   \item Spatial images and coordinates (for spatial experiments)
#'   \item Cell metadata and feature annotations
#'   \item Cell identity classes
#'   \item Command history
#'   \item Miscellaneous and tool-specific data
#' }
#'
#' The h5Seurat format is particularly useful for:
#' \itemize{
#'   \item Storing large datasets efficiently with HDF5 compression
#'   \item Rapid conversion to Python formats (h5ad, h5mu)
#'   \item Multi-modal and spatial transcriptomics experiments
#'   \item Preserving all Seurat V5 layer information
#' }
#'
#' @section Seurat V5 Layer Support:
#' When saving Seurat V5 objects with multiple layers (e.g., counts, data, scale.data),
#' all layers are preserved and can be selectively loaded using \code{\link{readH5Seurat}}.
#'
#' @seealso
#' \code{\link{readH5Seurat}} to load a saved h5Seurat file
#' \code{\link{as.h5Seurat}} for direct conversion without object assignment
#' \code{\link{scConvert}} for converting to other formats (h5ad, h5mu)
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scConvert)
#'
#' # Create a simple example Seurat object
#' seurat_obj <- CreateSeuratObject(counts = GetAssayData(pbmc_small))
#'
#' # Save to h5Seurat format
#' writeH5Seurat(seurat_obj, filename = "my_data.h5seurat")
#'
#' # Save with overwrite if file already exists
#' writeH5Seurat(seurat_obj, filename = "my_data.h5seurat", overwrite = TRUE)
#'
#' # Load the saved file back
#' seurat_obj <- readH5Seurat("my_data.h5seurat")
#'
#' # For multimodal data (e.g., CITE-seq)
#' # writeH5Seurat automatically saves all assays
#' writeH5Seurat(citeseq_obj, filename = "multimodal_data.h5seurat")
#'
#' # For spatial data (e.g., Visium)
#' writeH5Seurat(visium_obj, filename = "spatial_visium.h5seurat")
#' }
#'
#' @name writeH5Seurat
#' @rdname writeH5Seurat
#'
#' @export
#'
writeH5Seurat <- function(object, filename, overwrite = FALSE, verbose = TRUE,
                          gzip = NULL, ...) {
  if (!is.null(gzip)) {
    previous <- options(scConvert.compression.level = as.integer(gzip))
    on.exit(options(previous), add = TRUE)
  }
  UseMethod(generic = 'writeH5Seurat', object = object)
}

#' @return \code{as.h5Seurat}: an open \code{\link{h5Seurat}} connection
#'
#' @rdname writeH5Seurat
#' @export
#'
as.h5Seurat <- function(x, ...) {
  UseMethod(generic = 'as.h5Seurat', object = x)
}

# Coerce anything Seurat knows how to convert
.as_seurat_or_stop <- function(x, verbose, ...) {
  tryCatch(
    as.Seurat(x, verbose = verbose, ...),
    error = function(e) {
      stop("Cannot convert an object of class ", paste(class(x), collapse = ", "),
           " to a Seurat object", call. = FALSE)
    }
  )
}

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname writeH5Seurat
#' @method writeH5Seurat default
#' @export
#'
writeH5Seurat.default <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  object <- .as_seurat_or_stop(object, verbose = verbose, ...)
  if (missing(filename)) filename <- paste0(Project(object = object), '.h5Seurat')
  invisible(writeH5Seurat(object = object, filename = filename, overwrite = overwrite,
                          verbose = verbose, ...))
}

#' @rdname writeH5Seurat
#' @method writeH5Seurat Seurat
#' @export
#'
writeH5Seurat.Seurat <- function(object, filename = paste0(Project(object = object), '.h5Seurat'),
                                 overwrite = FALSE, verbose = TRUE, gzip = NULL, ...) {
  if (!is.null(gzip)) {
    previous <- options(scConvert.compression.level = as.integer(gzip))
    on.exit(options(previous), add = TRUE)
  }
  # opt-in compiled writer: options(scConvert.use_c_writer = TRUE)
  if (isTRUE(getOption("scConvert.use_c_writer")) &&
      is.loaded("C_write_h5seurat", PACKAGE = "scConvert")) {
    if (isTRUE(.writeH5Seurat_c(object, filename = filename, overwrite = overwrite,
                                verbose = verbose))) {
      return(invisible(filename))
    }
    if (verbose) message("C writer failed, falling back to R writer")
  }
  hfile <- as.h5Seurat(x = object, filename = filename, overwrite = overwrite,
                       verbose = verbose, ...)
  hfile$close_all()
  invisible(hfile$filename)
}

#' @rdname writeH5Seurat
#' @method as.h5Seurat default
#' @export
#'
as.h5Seurat.default <- function(x, filename, overwrite = FALSE, verbose = TRUE, ...) {
  x <- .as_seurat_or_stop(x, verbose = verbose, ...)
  if (missing(filename)) filename <- paste0(Project(object = x), '.h5Seurat')
  writeH5Seurat(object = x, filename = filename, overwrite = overwrite, verbose = verbose, ...)
}

#' @rdname writeH5Seurat
#' @method as.h5Seurat H5File
#' @export
#'
as.h5Seurat.H5File <- function(x, ...) {
  # reopen with the same intent as the existing handle
  h5Seurat$new(filename = x$filename,
               mode = if (Writeable(x = x, error = FALSE)) 'r+' else 'r')
}

#' @importFrom tools file_ext
#' @importFrom SeuratObject Assays Reductions Images Command Misc Tool Idents
#'
#' @rdname writeH5Seurat
#' @method as.h5Seurat Seurat
#' @export
#'
as.h5Seurat.Seurat <- function(x, filename = paste0(Project(object = x), '.h5seurat'),
                               overwrite = FALSE, verbose = TRUE, ...) {
  if (!grepl('^h5seurat$', file_ext(filename), ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5seurat')
  }
  if (file.exists(filename)) {
    if (!overwrite) stop("An h5Seurat file already exists at ", filename, call. = FALSE)
    warning("Overwriting ", filename, call. = FALSE, immediate. = TRUE)
    file.remove(filename)
  }
  obj_version <- as.character(slot(x, 'version'))
  if (verbose) message("Object version ", obj_version)
  if (numeric_version(obj_version) < numeric_version('3.1.2')) {
    stop("This Seurat object is too old to save; run UpdateSeuratObject() first ",
         "(at least v3.1.2 is needed)", call. = FALSE)
  }
  hfile <- h5Seurat$new(filename = filename, mode = 'w')
  tryCatch(hfile$set.version(version = obj_version), error = function(err) {
    hfile$close_all()
    file.remove(filename)
    stop(conditionMessage(err), call. = FALSE)
  })
  put <- function(value, name, where) {
    WriteH5Group(x = value, name = name, hgroup = where, verbose = verbose)
  }
  for (a in Assays(object = x)) put(x[[a]], a, hfile[['assays']])
  for (r in Reductions(object = x)) put(x[[r]], r, hfile[['reductions']])
  for (nb in names(slot(x, 'neighbors'))) put(x[[nb]], nb, hfile[['neighbors']])
  for (g in names(x)) {
    if (inherits(x[[g]], 'Graph')) put(x[[g]], g, hfile[['graphs']])
  }
  Project(object = hfile) <- Project(object = x)
  DefaultAssay(object = hfile) <- DefaultAssay(object = x)

  # Images: from the spatial-aware layout on, or v4 SliceImage objects
  images <- tryCatch(Images(object = x), error = function(e) character(0))
  spatial_layout <- numeric_version(obj_version) >= numeric_version(spatial.version)
  if (length(images) &&
      (spatial_layout || any(vapply(images, function(i) inherits(x[[i]], 'SliceImage'), logical(1))))) {
    for (img in images) {
      if (verbose) message("Adding image ", img)
      put(x[[img]], img, hfile[['images']])
    }
  }

  put(x[[]], 'meta.data', hfile)
  # cell names must be a plain 1-D character dataset
  cells <- as.character(colnames(x))
  hfile$create_dataset(name = 'cell.names', robj = cells, dtype = GuessDType(x = cells))
  put(Idents(object = x), 'active.ident', hfile)
  for (cmd in Command(object = x)) put(x[[cmd]], cmd, hfile[['commands']])
  for (m in names(Misc(object = x))) put(Misc(object = x, slot = m), m, hfile[['misc']])
  for (t in Tool(object = x)) put(Tool(object = x, slot = t), t, hfile[['tools']])
  hfile
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# C-accelerated h5Seurat writer
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Fast C-based h5Seurat writer
#'
#' Writes a Seurat object to h5Seurat format using native C HDF5 routines,
#' bypassing hdf5r's R6 method overhead (~10x faster).
#'
#' @param object Seurat object
#' @param filename Output file path
#' @param overwrite Allow overwriting existing file
#' @param verbose Show progress
#'
#' @return TRUE on success, FALSE on failure
#' @keywords internal
#'
.writeH5Seurat_c <- function(object, filename, overwrite = FALSE, verbose = TRUE) {
  if (!grepl('\\.h5seurat$', filename, ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5seurat')
  }
  if (file.exists(filename)) {
    if (overwrite) {
      file.remove(filename)
    } else {
      stop("H5Seurat file at ", filename, " already exists", call. = FALSE)
    }
  }

  assay_name <- DefaultAssay(object)

  # C writer only supports single-assay objects
  if (length(Assays(object)) > 1) {
    if (verbose) message("Multi-assay object, using R writer")
    return(FALSE)
  }

  assay_obj <- object[[assay_name]]

  # Extract sparse matrix (counts or data layer)
  mat <- NULL
  layer_name <- "counts"
  for (ln in c("counts", "data")) {
    tryCatch({
      m <- GetAssayData(assay_obj, layer = ln)
      if (inherits(m, "dgCMatrix") && length(m@x) > 0) {
        mat <- m
        layer_name <- ln
        break
      }
    }, error = function(e) NULL)
  }
  if (is.null(mat)) {
    if (verbose) message("C writer: no sparse matrix found, falling back")
    return(FALSE)
  }

  # Build mat list for C
  mat_list <- list(
    i = mat@i,
    p = mat@p,
    x = mat@x,
    dim = dim(mat),
    rownames = rownames(mat),
    colnames = colnames(mat)
  )
  # Add variable features if available
  vf <- tryCatch(VariableFeatures(object), error = function(e) character(0))
  if (length(vf) > 0) mat_list$variable.features <- vf

  # Extract metadata as named list
  meta <- as.list(object@meta.data)

  # Extract reductions as named list of matrices
  reductions <- list()
  for (rname in names(object@reductions)) {
    reduc <- object@reductions[[rname]]
    emb <- Embeddings(reduc)
    attr(emb, "key") <- Key(reduc)
    reductions[[rname]] <- emb
  }

  # Extract graphs as named list of sparse components
  graphs <- list()
  for (gname in names(object@graphs)) {
    g <- object@graphs[[gname]]
    if (inherits(g, "dgCMatrix") || inherits(g, "Graph")) {
      gmat <- as(g, "dgCMatrix")
      graphs[[gname]] <- list(
        i = gmat@i,
        p = gmat@p,
        x = gmat@x,
        dim = dim(gmat)
      )
    }
  }

  gzip_level <- GetCompressionLevel()
  if (verbose) message("Writing h5Seurat (C writer): ", filename)

  result <- .Call(C_write_h5seurat,
    filename, mat_list, meta, reductions, graphs,
    assay_name, as.integer(gzip_level)
  )

  if (isTRUE(result) && verbose) {
    message("  Written: ", ncol(mat), " cells, ", nrow(mat), " features")
  }
  return(isTRUE(result))
}
