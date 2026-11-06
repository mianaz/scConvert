#' @importFrom Seurat Project
#'
NULL

#' Save a Seurat object to a Loom file
#'
#' Export a Seurat object to Loom format (HDF5-based file format optimized for storing
#' annotated matrices). Loom files are compatible with the loompy Python package and
#' other tools in the bioinformatics community. This format is useful for sharing data
#' with Python-based analysis workflows or archiving analysis results.
#'
#' @param object A Seurat object to save
#' @param filename Name of file to save the object to. If not provided, defaults to
#'   \code{<ProjectName>.loom}. The .loom extension is added automatically if not present.
#' @param overwrite Logical; if \code{TRUE}, overwrite an existing file. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} (default), show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return \code{writeLoom}: Invisibly returns the filename of the saved file
#'
#' @details
#' The Loom format organizes data as follows:
#' \itemize{
#'   \item \code{/matrix}: Main expression matrix (features x cells)
#'   \item \code{/row_attrs}: Feature/gene-level annotations
#'   \item \code{/col_attrs}: Cell/sample-level metadata (cell names, cluster assignments, etc.)
#'   \item \code{/layers}: Additional expression layers if present
#' }
#'
#' When saving a Seurat object:
#' \itemize{
#'   \item Default assay data becomes the main matrix
#'   \item Cell metadata and feature annotations are preserved
#'   \item Dimensional reductions are stored under \code{/reductions}
#'   \item The SEURAT_ASSAY attribute stores the assay name for roundtrip loading
#' }
#'
#' @seealso
#' \code{\link{readLoom}} to load Loom files back as Seurat objects
#' \code{\link{writeH5Seurat}} to save in h5Seurat format
#' \code{\link{scConvert}} for converting between formats
#' \href{http://linnarssonlab.org/loompy/}{Loom documentation}
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scConvert)
#'
#' # Create a Seurat object (or use an existing one)
#' seurat_obj <- CreateSeuratObject(counts = pbmc_small$RNA@counts)
#'
#' # Save to Loom format
#' writeLoom(seurat_obj, filename = "my_data.loom")
#'
#' # Save with overwrite if needed
#' writeLoom(seurat_obj, filename = "my_data.loom", overwrite = TRUE)
#'
#' # Load it back
#' loaded_obj <- readLoom("my_data.loom")
#'
#' # For sharing with Python tools
#' writeLoom(seurat_obj, filename = "data_for_python.loom")
#' # Use in Python with: adata = loompy.connect("data_for_python.loom")
#' }
#'
#' @name writeLoom
#' @rdname writeLoom
#'
#' @export
#'
writeLoom <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  UseMethod(generic = 'writeLoom', object = object)
}

#' @param x A Seurat object to convert
#'
#' @return \code{as.loom}: an open \code{\link{loom}} connection
#'
#' @rdname writeLoom
#' @export
#'
as.loom <- function(x, ...) {
  UseMethod(generic = 'as.loom', object = x)
}

#' @rdname writeLoom
#' @method writeLoom default
#' @export
#'
writeLoom.default <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  object <- .as_seurat_or_stop(object, verbose = verbose, ...)
  if (missing(filename)) filename <- paste0(Project(object = object), '.loom')
  invisible(writeLoom(object = object, filename = filename, overwrite = overwrite,
                      verbose = verbose, ...))
}

#' @rdname writeLoom
#' @method writeLoom Seurat
#' @export
#'
writeLoom.Seurat <- function(object, filename = paste0(Project(object = object), '.loom'),
                             overwrite = FALSE, verbose = TRUE, ...) {
  con <- as.loom(x = object, filename = filename, overwrite = overwrite, verbose = verbose, ...)
  con$close_all()
  invisible(con$filename)
}

#' @rdname writeLoom
#' @method as.loom default
#' @export
#'
as.loom.default <- function(x, filename, overwrite = FALSE, verbose = TRUE, ...) {
  x <- .as_seurat_or_stop(x, verbose = verbose, ...)
  if (missing(filename)) filename <- paste0(Project(object = x), '.loom')
  as.loom(x = x, filename = filename, overwrite = overwrite, verbose = verbose, ...)
}

#' @rdname writeLoom
#' @method as.loom H5File
#' @export
#'
as.loom.H5File <- function(x, ...) {
  for (part in c('matrix', 'layers', 'row_attrs', 'col_attrs')) {
    if (!x$exists(name = part)) stop("Not a loom file: '", part, "' is missing", call. = FALSE)
  }
  loom$new(filename = x$filename,
           mode = if (x$get_file_id()$is_writeable()) 'r+' else 'r')
}

# /reductions/<name>: embeddings (components x cells), loadings, stdev
.loom_write_reduction <- function(lfile, red, dr) {
  base <- paste0('reductions/', red)
  lfile$create_group(name = base)
  emb <- Embeddings(object = dr)
  lfile$create_dataset(name = paste0(base, '/embeddings'), robj = t(emb),
                       dtype = GuessDType(x = emb[1, 1]))
  ld <- tryCatch(Loadings(object = dr), error = function(e) NULL)
  if (!is.null(ld) && !IsMatrixEmpty(x = ld)) {
    lfile$create_dataset(name = paste0(base, '/loadings'), robj = ld,
                         dtype = GuessDType(x = ld[1, 1]))
  }
  sdev <- tryCatch(Stdev(object = dr), error = function(e) NULL)
  if (length(sdev)) {
    lfile$create_dataset(name = paste0(base, '/stdev'), robj = sdev, dtype = GuessDType(x = sdev))
  }
  invisible(NULL)
}

# Factor levels (and ordering) of cell metadata columns, kept outside
# col_attrs: the loom spec only allows 1-D datasets there and loompy /
# scanpy reject sub-groups
.loom_write_factor_levels <- function(lfile, factors, verbose) {
  if (!length(factors)) return(invisible(NULL))
  tryCatch({
    root_path <- 'scConvert_extensions/col_factor_levels'
    for (p in c('scConvert_extensions', root_path)) {
      if (!lfile$exists(p)) lfile$create_group(p)
    }
    root <- lfile[[root_path]]
    for (col in names(factors)) {
      if (root$exists(col)) root$link_delete(col)
      g <- root$create_group(col)
      g$create_dataset(name = 'levels', robj = levels(factors[[col]]))
      g$create_dataset(name = 'ordered', robj = as.integer(is.ordered(factors[[col]])))
    }
    if (verbose) message("Preserved factor levels for ", length(factors), " meta.data column(s)")
  }, error = function(e) {
    warning("Could not write factor levels to loom file: ", conditionMessage(e), call. = FALSE)
  })
  invisible(NULL)
}

#' @importFrom tools file_ext
#' @importFrom SeuratObject Assays Reductions Embeddings Loadings Stdev Graphs GetAssayData
#'
#' @rdname writeLoom
#' @method as.loom Seurat
#' @export
#'
as.loom.Seurat <- function(x, filename = paste0(Project(object = x), '.loom'),
                           overwrite = FALSE, verbose = TRUE, ...) {
  if (!grepl('^loom$', file_ext(filename))) filename <- paste0(filename, '.loom')
  if (file.exists(filename)) {
    if (!isTRUE(overwrite)) stop("A loom file already exists at ", filename, call. = FALSE)
    warning("Overwriting ", filename, call. = FALSE, immediate. = TRUE)
    file.remove(filename)
  }
  lfile <- loom$new(filename = filename, mode = 'w')
  default <- DefaultAssay(object = x)
  layer_of <- function(layer, assay) {
    tryCatch(GetAssayData(object = x, layer = layer, assay = assay), error = function(e) NULL)
  }
  # counts and scale.data of an assay as loom layers, when they are full size
  add_extra_layers <- function(assay) {
    for (layer in c('counts', 'scale.data')) {
      m <- layer_of(layer, assay)
      if (is.null(m) || !identical(dim(m), dim(x))) next
      if (verbose) message("Adding layer ", layer, " for assay ", assay)
      lfile$add_layer(x = m, name = if (assay == default) layer else paste(assay, layer, sep = '_'),
                      verbose = verbose)
    }
  }

  # /matrix: normalized data of the default assay, else its counts
  main <- layer_of('data', default)
  if (is.null(main) || IsMatrixEmpty(x = main)) {
    main <- layer_of('counts', default)
    if (verbose && !is.null(main)) message("Data layer empty, using counts from ", default, " as /matrix")
  } else if (verbose) {
    message("Saving data from ", default, " as /matrix")
  }
  if (is.null(main) || IsMatrixEmpty(x = main)) {
    stop("No data or counts layer found in default assay", call. = FALSE)
  }
  .loom_write_matrix(lfile = lfile, x = main, name = 'matrix', verbose = verbose)
  add_extra_layers(default)
  # other assays with the same shape go to layers named after the assay
  others <- Filter(function(a) identical(dim(x[[a]]), dim(x)), setdiff(Assays(object = x), default))
  for (a in others) {
    if (verbose) message("Adding data for ", a)
    lfile$add_layer(x = layer_of('data', a), name = a, verbose = verbose)
    add_extra_layers(a)
  }

  reds <- Reductions(object = x)
  if (length(reds) && !lfile$exists('reductions')) lfile$create_group(name = 'reductions')
  for (red in reds) {
    if (verbose) message("Adding dimensional reduction ", red)
    tryCatch(.loom_write_reduction(lfile, red, x[[red]]), error = function(e) {
      warning("Failed to save dimensional reduction ", red, ": ", conditionMessage(e),
              call. = FALSE, immediate. = TRUE)
    })
  }

  for (g in tryCatch(Graphs(object = x), error = function(e) character(0))) {
    tryCatch(lfile$add_graph(x = x[[g]], name = g, type = 'col'), error = function(e) {
      warning("Could not add graph ", g, call. = FALSE, immediate. = TRUE)
    })
  }

  # cell attributes; read from @meta.data directly because x[[col]] drops
  # the 'ordered' class
  lfile$add_attribute(x = colnames(x), name = 'CellID', type = 'col')
  md <- x@meta.data
  for (col in colnames(md)) lfile$add_attribute(x = md[[col]], name = col, type = 'col')
  .loom_write_factor_levels(lfile, Filter(is.factor, as.list(md)), verbose)

  # gene attributes: names plus the default assay's feature metadata
  lfile$add_attribute(x = rownames(x), name = 'Gene', type = 'row')
  if (default %in% Assays(object = x)) {
    fm <- x[[default]][[]]
    for (col in colnames(fm)) {
      tryCatch(lfile$add_attribute(x = fm[[col]], name = col, type = 'row'), error = function(e) {
        if (verbose) {
          warning("Failed to save feature metadata '", col, "': ", conditionMessage(e),
                  call. = FALSE, immediate. = TRUE)
        }
      })
    }
  }
  lfile
}
