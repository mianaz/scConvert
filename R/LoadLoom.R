#' @include zzz.R
#' @include loom.R
#' @include Bindings.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reading loom files
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load a loom file as a Seurat object
#'
#' Reads the main matrix, optional normalized/scaled layers, gene and cell
#' attributes, cell graphs, and the dimensional reductions and factor levels
#' that \code{\link{writeLoom}} stores, into a Seurat object.
#'
#' @param file,x Path to a loom file, or an open \code{\link{loom}} (or
#'   \code{H5File}) connection
#' @param assay Assay to store the data in; when \code{NULL}, the assay
#'   recorded by \code{writeLoom} (\code{/attrs/SEURAT_ASSAY}) or
#'   \code{"RNA"}
#' @param cells,features Datasets in \code{/col_attrs} and \code{/row_attrs}
#'   holding cell and gene names; both must exist
#' @param normalized Layer in \code{/layers} to load as normalized data, or
#'   \code{"matrix"} to store the main matrix as data instead of counts
#' @param scaled Layer in \code{/layers} to load as scaled data
#' @param filter Apply the \code{Valid} attributes: to \code{"cells"}
#'   (default), \code{"features"}, \code{"all"} or \code{"none"}
#' @param verbose Show progress messages
#' @param ... Passed on to other methods
#'
#' @return A \code{Seurat} object
#'
#' @details
#' A loom file stores the expression matrix in \code{/matrix} (genes x cells
#' in loompy's view), extra matrices in \code{/layers}, gene annotations in
#' \code{/row_attrs}, cell annotations in \code{/col_attrs} and cell graphs
#' (coordinate triplets) in \code{/col_graphs}. Loom 0.1 and 3.0 files are
#' read with the same rules.
#'
#' @seealso \code{\link{writeLoom}},
#'   \href{http://linnarssonlab.org/loompy/format/}{the loom file format}
#'
#' @examples
#' \dontrun{
#' obj <- readLoom("data.loom")
#' obj <- readLoom("data.loom", assay = "RNA", filter = "all")
#' obj <- readLoom("data.loom", normalized = "normalized", scaled = "scaled")
#' }
#'
#' @name readLoom
#' @rdname readLoom
#' @export
#'
readLoom <- function(file, assay = NULL, cells = 'CellID', features = 'Gene',
                     normalized = NULL, scaled = NULL,
                     filter = c('cells', 'features', 'all', 'none'), verbose = TRUE, ...) {
  UseMethod(generic = 'readLoom', object = file)
}

#' @rdname readLoom
#' @method readLoom character
#' @export
#'
readLoom.character <- function(file, ...) {
  con <- loom$new(filename = file, mode = 'r')
  on.exit(con$close_all())
  readLoom(file = con, ...)
}

#' @rdname readLoom
#' @method readLoom H5File
#' @export
#'
readLoom.H5File <- function(file, ...) {
  readLoom(file = as.loom(x = file), ...)
}

#' @rdname readLoom
#' @method readLoom loom
#' @export
#'
readLoom.loom <- function(file, ...) {
  as.Seurat(x = file, ...)
}

# Resolve a user-supplied layer name to "layers/<name>"; NULL (with a
# warning) when it is not in the file
.loom_layer_path <- function(x, name, what) {
  if (is.null(name)) return(NULL)
  path <- sub('^/', '', name)
  if (!startsWith(path, 'layers')) path <- paste0('layers/', path)
  if (lengths(regmatches(path, gregexpr('/', path))) != 1L) {
    stop("'", what, "' must name a dataset directly inside /layers", call. = FALSE)
  }
  if (!Exists(x = x, name = path)) {
    warning(sprintf("No layer %s in this loom file; %s data not loaded", basename(path), what),
            call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  path
}

#' @importFrom SeuratObject as.Seurat
#' @rdname readLoom
#' @method as.Seurat loom
#' @export
#'
as.Seurat.loom <- function(x, assay = NULL, cells = 'CellID', features = 'Gene',
                           normalized = NULL, scaled = NULL,
                           filter = c('cells', 'features', 'all', 'none'), verbose = TRUE, ...) {
  cells <- paste0('col_attrs/', cells)
  features <- paste0('row_attrs/', features)
  for (p in c(features, cells)) {
    if (!Exists(x = x, name = p) || !inherits(x[[p]], 'H5D')) {
      stop("No names dataset at ", p, call. = FALSE)
    }
    if (length(Dims(x = x[[p]])) != 1L) stop(p, " must be one-dimensional", call. = FALSE)
  }
  matrix_as_data <- !is.null(normalized) && grepl('^/?matrix$', normalized)
  if (!matrix_as_data) normalized <- .loom_layer_path(x, normalized, 'normalized')
  scaled <- .loom_layer_path(x, scaled, 'scaled')
  ClosestVersion(query = x$version(), targets = c('0.1.0', '3.0.0'))  # errors when unsupported
  .loom_to_seurat(x, assay = assay, cells = cells, features = features,
                  normalized = if (matrix_as_data) 'matrix' else normalized,
                  scaled = scaled, filter = match.arg(filter[1], c('cells', 'features', 'all', 'none')),
                  verbose = verbose)
}

# Read a genes x cells matrix (the loom /matrix or a layer) as an R matrix
.loom_read_matrix <- function(x, path, keep_features, keep_cells) {
  obj <- x[[path]]
  m <- if (inherits(obj, 'H5Group')) {
    # sparse stand-in (CSC features x cells) used by streaming conversion
    Matrix::sparseMatrix(i = obj[['indices']][] + 1, p = obj[['indptr']][],
                         x = as.numeric(obj[['data']][]),
                         dims = as.integer(h5attr(obj, 'shape')))
  } else {
    t(as.matrix(x = obj))
  }
  if (!is.null(keep_features)) m <- m[keep_features, , drop = FALSE]
  if (!is.null(keep_cells)) m <- m[, keep_cells, drop = FALSE]
  m
}

# Attribute table (row_attrs / col_attrs) without the name and Valid columns
.loom_attr_table <- function(x, group, drop, keep, names) {
  df <- as.data.frame(x = x[[group]])
  df <- df[, !colnames(df) %in% c(drop, 'Valid'), drop = FALSE]
  if (!is.null(keep)) df <- df[keep, , drop = FALSE]
  rownames(df) <- names
  df
}

# Restore factor columns from the scConvert extension group written by
# writeLoom (kept outside col_attrs for compatibility with the spec)
.loom_restore_factors <- function(x, md) {
  if (!Exists(x = x, name = 'scConvert_extensions/col_factor_levels')) return(md)
  tryCatch({
    root <- x[['scConvert_extensions/col_factor_levels']]
    for (col in intersect(names(root), colnames(md))) {
      g <- root[[col]]
      if (!Exists(x = g, name = 'levels')) next
      ordered <- Exists(x = g, name = 'ordered') && isTRUE(as.integer(g[['ordered']]$read()) == 1L)
      md[[col]] <- factor(as.character(md[[col]]), levels = as.character(g[['levels']]$read()),
                          ordered = ordered)
    }
    md
  }, error = function(e) {
    warning("Could not restore factor levels from loom: ", conditionMessage(e), call. = FALSE)
    md
  })
}

# Dimensional reductions stored by writeLoom under the (non-standard)
# /reductions group
.loom_add_reductions <- function(x, object, assay, keep_cells, keep_features, verbose) {
  if (!Exists(x = x, name = 'reductions')) return(object)
  for (red in names(x[['reductions']])) {
    object <- tryCatch({
      grp <- x[['reductions']][[red]]
      if (!Exists(x = grp, name = 'embeddings')) next
      if (verbose) message("Loading dimensional reduction ", red)
      emb <- t(as.matrix(x = grp[['embeddings']]))
      if (!is.null(keep_cells)) emb <- emb[keep_cells, , drop = FALSE]
      rownames(emb) <- colnames(object)
      dr <- CreateDimReducObject(embeddings = emb, key = paste0(red, '_'), assay = assay)
      if (Exists(x = grp, name = 'loadings')) {
        ld <- as.matrix(x = grp[['loadings']])
        if (!is.null(keep_features) && nrow(ld) == length(keep_features)) {
          ld <- ld[keep_features, , drop = FALSE]
        }
        Loadings(object = dr) <- ld
      }
      if (Exists(x = grp, name = 'stdev')) {
        tryCatch(slot(dr, 'stdev') <- grp[['stdev']][], error = function(e) {
          warning("Could not set stdev for reduction ", red, call. = FALSE, immediate. = TRUE)
        })
      }
      object[[red]] <- dr
      object
    }, error = function(e) {
      warning("Failed to load dimensional reduction ", red, ": ", conditionMessage(e),
              call. = FALSE, immediate. = TRUE)
      object
    })
  }
  object
}

#' @importFrom SeuratObject CreateAssayObject CreateSeuratObject SetAssayData Key<-
#'   DefaultAssay<- CreateDimReducObject Loadings<- as.Graph AddMetaData
#' @importFrom methods slot<-
#' @keywords internal
#' @noRd
.loom_to_seurat <- function(x, assay = NULL, cells = 'col_attrs/CellID',
                            features = 'row_attrs/Gene', normalized = NULL, scaled = NULL,
                            filter = 'cells', verbose = TRUE) {
  assay <- assay %||% suppressWarnings(DefaultAssay(object = x)) %||% 'RNA'
  # loom "Valid" flags select cells / genes
  valid <- function(axis) {
    path <- paste0(axis, '_attrs/Valid')
    applies <- filter == 'all' || (axis == 'col' && filter == 'cells') ||
      (axis == 'row' && filter == 'features')
    if (!applies || !Exists(x = x, name = path)) return(NULL)
    keep <- as.logical(x[[path]][])
    if (verbose) {
      message(sprintf("Found %s filter, keeping %d of %d %s",
                      if (axis == 'col') 'cell' else 'feature', sum(keep), length(keep),
                      if (axis == 'col') 'cells' else 'features'))
    }
    keep
  }
  keep_cells <- valid('col')
  keep_features <- valid('row')

  if (verbose) message("Reading in /matrix")
  mat <- .loom_read_matrix(x, 'matrix', keep_features, keep_cells)
  gene_names <- x[[features]][]
  cell_names <- x[[cells]][]
  if (!is.null(keep_features)) gene_names <- gene_names[keep_features]
  if (!is.null(keep_cells)) cell_names <- cell_names[keep_cells]
  if (anyDuplicated(gene_names)) {
    warning("Duplicate feature names found, making unique", call. = FALSE, immediate. = TRUE)
    gene_names <- make.unique(gene_names)
  }
  if (anyDuplicated(cell_names)) {
    warning("Duplicate cell names found, making unique", call. = FALSE, immediate. = TRUE)
    cell_names <- make.unique(cell_names)
  }
  dn <- list(features = gene_names, cells = cell_names)
  dimnames(mat) <- dn

  matrix_as_data <- identical(normalized, 'matrix')
  if (verbose) message("Storing /matrix as ", if (matrix_as_data) "data" else "counts")
  assay_obj <- if (matrix_as_data) CreateAssayObject(data = mat) else CreateAssayObject(counts = mat)
  Key(object = assay_obj) <- UpdateKey(key = tolower(assay))
  if (verbose) message("Saving /matrix to assay '", assay, "'")
  object <- CreateSeuratObject(counts = assay_obj, assay = assay)

  for (part in list(list(path = if (!matrix_as_data) normalized, slot = 'data', what = 'data'),
                    list(path = scaled, slot = 'scale.data', what = 'scaled data'))) {
    if (is.null(part$path)) next
    if (verbose) message("Saving ", basename(part$path), " as ", part$what)
    m <- .loom_read_matrix(x, part$path, keep_features, keep_cells)
    dimnames(m) <- dn
    object <- SetAssayData(object = object, assay = assay, slot = part$slot, new.data = m)
  }

  gene_meta <- .loom_attr_table(x, 'row_attrs', basename(features), keep_features, gene_names)
  if (ncol(gene_meta)) object[[assay]] <- AddMetaData(object = object[[assay]], metadata = gene_meta)
  cell_meta <- .loom_attr_table(x, 'col_attrs', basename(cells), keep_cells, cell_names)
  cell_meta <- .loom_restore_factors(x, cell_meta)
  if (ncol(cell_meta)) object <- AddMetaData(object = object, metadata = cell_meta)

  object <- .loom_add_reductions(x, object, assay, keep_cells, keep_features, verbose)

  n_all <- length(x[[cells]][])
  for (g in names(x[['col_graphs']])) {
    if (verbose) message("Loading graph ", g)
    adj <- LoadGraph(graph = x[['col_graphs']][[g]], n = n_all)
    if (!is.null(keep_cells)) adj <- adj[keep_cells, keep_cells, drop = FALSE]
    dimnames(adj) <- list(colnames(object), colnames(object))
    adj <- as.Graph(x = adj)
    DefaultAssay(object = adj) <- assay
    object[[g]] <- adj
  }
  object
}

# Sparse matrix from a loom graph group (0-based coordinate triplets a, b
# and weights w). `n` is the number of cells; without it the matrix would
# be cut short when the last cells have no edges.
LoadGraph <- function(graph, n = NULL) {
  if (!inherits(graph, 'H5Group')) {
    stop(graph$get_obj_name(), " is not a loom graph", call. = FALSE)
  }
  for (part in c('a', 'b', 'w')) {
    ok <- Exists(x = graph, name = part) && inherits(graph[[part]], 'H5D') &&
      length(Dims(x = graph[[part]])) == 1L
    if (!ok) stop("Loom graph is missing the one-dimensional dataset '", part, "'", call. = FALSE)
  }
  i <- graph[['a']][] + 1
  j <- graph[['b']][] + 1
  if (is.null(n)) n <- max(c(i, j, 0))
  Matrix::sparseMatrix(i = i, j = j, x = graph[['w']][], dims = c(n, n))
}
