#' @include WriteH5Group.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
writeH5Seurat <- function(
  object,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  gzip = NULL,
  ...
) {
  # Optimization 3: configurable gzip level (0 = no compression, ~50-70% faster writes)
  if (!is.null(gzip)) {
    old_gzip <- getOption("scConvert.compression.level")
    on.exit(options("scConvert.compression.level" = old_gzip), add = TRUE)
    options("scConvert.compression.level" = as.integer(gzip))
  }
  UseMethod(generic = 'writeH5Seurat', object = object)
}

#' @return \code{as.h5Seurat}: An \code{\link{h5Seurat}} object
#'
#' @rdname writeH5Seurat
#'
#' @export
#'
as.h5Seurat <- function(x, ...) {
  UseMethod(generic = 'as.h5Seurat', object = x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname writeH5Seurat
#' @method writeH5Seurat default
#' @export
#'
writeH5Seurat.default <- function(
  object,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  tryCatch(
    expr = object <- as.Seurat(object = object, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = object), collapse = ', '),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = object), '.h5Seurat')
  }
  return(invisible(x = writeH5Seurat(
    object = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )))
}

#' @importFrom Seurat Project
#'
#' @rdname writeH5Seurat
#' @method writeH5Seurat Seurat
#' @export
#'
writeH5Seurat.Seurat <- function(
  object,
  filename = paste0(Project(object = object), '.h5Seurat'),
  overwrite = FALSE,
  verbose = TRUE,
  gzip = NULL,
  ...
) {
  # Optimization 3: configurable gzip
  if (!is.null(gzip)) {
    old_gzip <- getOption("scConvert.compression.level")
    on.exit(options("scConvert.compression.level" = old_gzip), add = TRUE)
    options("scConvert.compression.level" = as.integer(gzip))
  }
  # Try C writer (Optimization: ~10x faster with gzip=0)
  # Enable with options(scConvert.use_c_writer = TRUE)
  c_available <- isTRUE(getOption("scConvert.use_c_writer")) &&
                 is.loaded("C_write_h5seurat", PACKAGE = "scConvert")
  if (c_available) {
    result <- .writeH5Seurat_c(object, filename = filename, overwrite = overwrite,
                                verbose = verbose)
    if (isTRUE(result)) return(invisible(filename))
    if (verbose) message("C writer failed, falling back to R writer")
  }
  h5seurat <- as.h5Seurat(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  h5seurat$close_all()
  return(invisible(x = h5seurat$filename))
}

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname writeH5Seurat
#' @method as.h5Seurat default
#' @export
#'
as.h5Seurat.default <- function(
  x,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  tryCatch(
    expr = x <- as.Seurat(x = x, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = x), collapse = ', '),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = x), '.h5Seurat')
  }
  return(writeH5Seurat(
    object = x,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname writeH5Seurat
#' @method as.h5Seurat H5File
#' @export
#'
as.h5Seurat.H5File <- function(x, ...) {
  return(h5Seurat$new(
    filename = x$filename,
    mode = ifelse(test = Writeable(x = x, error = FALSE), yes = 'r', no = 'r+')
  ))
}

#' @importFrom tools file_ext
#' @importFrom Seurat Project Assays Reductions DefaultAssay<- DefaultAssay
#' Idents Command Misc Tool
#'
#' @rdname writeH5Seurat
#' @method as.h5Seurat Seurat
#' @export
#'
as.h5Seurat.Seurat <- function(
  x,
  filename = paste0(Project(object = x), '.h5seurat'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (!grepl(pattern = '^h5seurat$', x = file_ext(x = filename), ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5seurat')
  }
  if (file.exists(filename)) {
    if (overwrite) {
      warning(
        "Overwriting previous file ",
        filename,
        call. = FALSE,
        immediate. = TRUE
      )
      file.remove(filename)
    } else {
      stop("H5Seurat file at ", filename, " already exists", call. = FALSE)
    }
  }
  # Get the object version
  object.version <- as.character(x = slot(object = x, name = 'version'))

  # Display version message before creating file
  if (verbose) {
    message("Object version ", object.version)
  }

  # Create h5Seurat file (validation happens automatically but messages are suppressed in create method)
  hfile <- h5Seurat$new(filename = filename, mode = 'w')

  # Set the version to the object version
  tryCatch(
    expr = hfile$set.version(version = object.version),
    error = function(err) {
      file.remove(hfile$filename)
      test.message <- paste0(
        "All target versions greater than query version (",
        object.version,
        ")"
      )
      if (err$message == test.message) {
        stop(
          "Object too old to save, please update your Seurat object to at least v3.1.2 using UpdateSeuratObject",
          call. = FALSE
        )
      }
      stop(err$message, call. = FALSE)
    }
  )
  # Add Assays
  for (assay in Assays(object = x)) {
    WriteH5Group(
      x = x[[assay]],
      name = assay,
      hgroup = hfile[['assays']],
      verbose = verbose
    )
  }
  # Add DimReducs
  for (reduc in Reductions(object = x)) {
    WriteH5Group(
      x = x[[reduc]],
      name = reduc,
      hgroup = hfile[['reductions']],
      verbose = verbose
    )
  }
  # Add Neighbors
  neighbors <- names(x = slot(object = x, name = "neighbors"))
  for (neighbor in neighbors) {
    WriteH5Group(
      x = x[[neighbor]],
      name = neighbor,
      hgroup = hfile[['neighbors']],
      verbose = verbose
    )
  }
  # Add Graphs
  graphs <- Filter(
    f = function(g) {
      return(inherits(x = x[[g]], what = 'Graph'))
    },
    x = names(x = x)
  )
  for (graph in graphs) {
    WriteH5Group(
      x = x[[graph]],
      name = graph,
      hgroup = hfile[['graphs']],
      verbose = verbose
    )
  }
  # Add attributes for project and default assay
  Project(object = hfile) <- Project(object = x)
  DefaultAssay(object = hfile) <- DefaultAssay(object = x)
  # Add Images (if object has spatial data or SliceImage objects)
  has_slice_images <- FALSE
  if (package_version(x = object.version) < package_version(x = spatial.version)) {
    has_slice_images <- tryCatch({
      imgs <- Seurat::Images(object = x)
      if (length(imgs) > 0) {
        any(sapply(imgs, function(img) inherits(x[[img]], 'SliceImage')))
      } else {
        FALSE
      }
    }, error = function(e) {
      FALSE
    })
  }

  if (package_version(x = object.version) >= package_version(x = spatial.version) || has_slice_images) {
    has_images <- tryCatch({
      length(Seurat::Images(object = x)) > 0
    }, error = function(e) {
      FALSE
    })

    if (has_images) {
      for (image in Seurat::Images(object = x)) {
        if (verbose) {
          message("Adding image ", image)
        }
        WriteH5Group(
          x = x[[image]],
          name = image,
          hgroup = hfile[['images']],
          verbose = verbose
        )
      }
    }
  }
  # Add metadata, cell names, and identity classes
  WriteH5Group(x = x[[]], name = 'meta.data', hgroup = hfile, verbose = verbose)
  # Ensure cell names are written as a 1D vector, not a 2D matrix
  # The issue causing the dimension mismatch when reading h5Seurat files is here
  cell_names <- colnames(x = x)
  
  # Force cell names to be a standard character vector (not a matrix or data.frame)
  if (!is.null(dim(cell_names)) || !is.vector(cell_names)) {
    if (verbose) {
      message("Converting cell names to a 1D vector for compatibility")
    }
    cell_names <- as.character(cell_names)
  }
  
  # Create the dataset directly instead of using WriteH5Group to ensure proper dimensionality
  hfile$create_dataset(
    name = 'cell.names',
    robj = cell_names,
    dtype = GuessDType(x = cell_names)
  )
  WriteH5Group(
    x = Idents(object = x),
    name = 'active.ident',
    hgroup = hfile,
    verbose = verbose
  )
  # Add SeuratCommands
  for (cmd in Command(object = x)) {
    WriteH5Group(
      x = x[[cmd]],
      name = cmd,
      hgroup = hfile[['commands']],
      verbose = verbose
    )
  }
  # Add miscellaneous data
  for (misc in names(x = Misc(object = x))) {
    WriteH5Group(
      x = Misc(object = x, slot = misc),
      name = misc,
      hgroup = hfile[['misc']],
      verbose = verbose
    )
  }
  # Add tool data
  for (tool in Tool(object = x)) {
    WriteH5Group(
      x = Tool(object = x, slot = tool),
      name = tool,
      hgroup = hfile[['tools']],
      verbose = verbose
    )
  }
  return(hfile)
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

  result <- .Call("C_write_h5seurat",
    filename, mat_list, meta, reductions, graphs,
    assay_name, as.integer(gzip_level),
    PACKAGE = "scConvert"
  )

  if (isTRUE(result) && verbose) {
    message("  Written: ", ncol(mat), " cells, ", nrow(mat), " features")
  }
  return(isTRUE(result))
}
