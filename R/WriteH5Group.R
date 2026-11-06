#' @include H5Utils.R
#' @include V5LayerSupport.R
#' @importFrom methods setGeneric setMethod slot slotNames is
#' @importFrom SeuratObject DefaultAssay Embeddings Loadings Key IsGlobal Stdev JS
#'   Misc VariableFeatures S4ToList
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Serialising R / Seurat objects into an h5Seurat file
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Layout conventions (shared with the readers in AssembleObject.R and
# H5Coerce.R):
#   atomic vector        dataset; a named vector also gets "<name>__names__"
#   logical              integer dataset (0/1, NA = 2) with s3class "logical"
#   factor               group: values (1-based codes), levels, [ordered]
#   list                 group of members, "names" attribute, s3class when
#                        the list carries a class
#   data.frame           group of columns, "colnames"/"logicals" attributes,
#                        row names in "_index"
#   dgCMatrix / Graph    CSC group: indices, indptr, data + "dims"
#   other S4 objects     group of slots + "s4class" attribute

# Datasets longer than this are gzip-compressed (when compression is on)
.H5_COMPRESS_MIN <- 64L

# Create a dataset, compressing long ones
.h5_dataset <- function(hgroup, name, value, dtype = GuessDType(x = value), ...) {
  level <- GetCompressionLevel()
  if (level > 0L && length(value) > .H5_COMPRESS_MIN) {
    hgroup$create_dataset(name = name, robj = value, dtype = dtype, gzip_level = level, ...)
  } else {
    hgroup$create_dataset(name = name, robj = value, dtype = dtype, ...)
  }
  invisible(NULL)
}

# Add an attribute to an HDF5 object
.h5_attr <- function(obj, name, value, dtype = GuessDType(x = value), ...) {
  obj$create_attr(attr_name = name, robj = value, dtype = dtype, ...)
  invisible(NULL)
}

#' Write an R object into an h5Seurat file
#'
#' Serialises vectors, factors, lists, data frames, sparse matrices and the
#' Seurat classes (assays, reductions, graphs, neighbors, command logs,
#' images) into an HDF5 group using the h5Seurat layout.
#'
#' @param x Object to write
#' @param name Name of the new dataset or group
#' @param hgroup \code{\link[hdf5r]{H5File}} or \code{\link[hdf5r]{H5Group}}
#'   to write into
#' @param verbose Report progress
#'
#' @return Invisibly \code{NULL}
#'
#' @name WriteH5Group
#' @rdname WriteH5Group
#' @keywords internal
#'
setGeneric(
  name = 'WriteH5Group',
  def = function(x, name, hgroup, verbose = TRUE) {
    if (!inherits(hgroup, c('H5File', 'H5Group'))) {
      stop("'hgroup' must be an hdf5r file or group", call. = FALSE)
    }
    standardGeneric('WriteH5Group')
  },
  signature = c('x')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Basic R objects
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Vectors, arrays and lists
BasicWrite <- function(x, name, hgroup, verbose = TRUE) {
  if (is.null(x)) return(invisible(NULL))
  if (is.data.frame(x)) {
    return(WriteH5Group(x = x, name = name, hgroup = hgroup, verbose = verbose))
  }
  if (is.list(x)) {
    x <- PadNames(x = x)
    grp <- hgroup$create_group(name = name)
    for (i in seq_along(x)) {
      WriteH5Group(x = x[[i]], name = names(x)[i], hgroup = grp, verbose = verbose)
    }
    if (length(names(x))) {
      # members that could not be written (e.g. NULL) are left out
      .h5_attr(grp, 'names', intersect(names(x), names(grp)), dtype = GuessDType(x = names(x)[1]))
    }
    if (!all(class(x) == 'list')) .h5_attr(grp, 's3class', class(x), dtype = GuessDType(x = class(x)[1]))
    return(invisible(NULL))
  }
  if (is.vector(x) && !is.null(names(x))) {
    .h5_dataset(hgroup, paste0(name, '__names__'), names(x))
  }
  .h5_dataset(hgroup, name, x)
}

#' @rdname WriteH5Group
setMethod('WriteH5Group', signature = c(x = 'array'), definition = BasicWrite)

#' @rdname WriteH5Group
setMethod('WriteH5Group', signature = c(x = 'list'), definition = BasicWrite)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'logical'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    WriteH5Group(x = BoolToInt(x = x), name = name, hgroup = hgroup, verbose = verbose)
    .h5_attr(hgroup[[name]], 's3class', 'logical', dtype = CachedGuessDType(x = 'logical'))
  }
)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'factor'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    grp <- hgroup$create_group(name = name)
    WriteH5Group(x = as.integer(x), name = 'values', hgroup = grp, verbose = verbose)
    WriteH5Group(x = levels(x), name = 'levels', hgroup = grp, verbose = verbose)
    if (is.ordered(x)) {
      .h5_attr(grp, 'ordered', 1L, space = ScalarSpace())
    }
    invisible(NULL)
  }
)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'data.frame'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    cols <- colnames(x)
    is_lgl <- vapply(cols, function(cl) is.logical(x[[cl]]), logical(1), USE.NAMES = FALSE)
    is_fct <- vapply(cols, function(cl) is.factor(x[[cl]]), logical(1), USE.NAMES = FALSE)
    as_group <- any(is_fct) ||
      isTRUE(getOption('scConvert.dtypes.dataframe_as_group', default = TRUE))
    if (!as_group) {
      # a single compound dataset; logical columns become integer codes
      for (cl in cols[is_lgl]) x[[cl]] <- BoolToInt(x = x[[cl]])
      hgroup$create_dataset(name = name, robj = x, dtype = GuessDType(x = x))
      if (any(is_lgl)) {
        stored <- hgroup[[name]]$get_type()$get_cpd_labels()
        .h5_attr(hgroup[[name]], 'logicals', intersect(cols[is_lgl], stored),
                 dtype = GuessDType(x = cols))
      }
      return(invisible(NULL))
    }
    grp <- hgroup$create_group(name = name)
    for (cl in cols) {
      WriteH5Group(x = x[[cl]], name = cl, hgroup = grp, verbose = verbose)
    }
    .h5_attr(grp, 'colnames', intersect(cols, names(grp)), dtype = GuessDType(x = cols))
    if (any(is_lgl)) {
      .h5_attr(grp, 'logicals', intersect(cols[is_lgl], names(grp)), dtype = GuessDType(x = cols))
    }
    if (length(rownames(x))) {
      .h5_attr(grp, '_index', '_index')
      WriteH5Group(x = rownames(x), name = '_index', hgroup = grp, verbose = verbose)
    }
    invisible(NULL)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Matrices and graphs
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CSC sparse matrix as a group of indices / indptr / data plus "dims"
SparseWrite <- function(x, name, hgroup, verbose = TRUE) {
  grp <- hgroup$create_group(name = name)
  level <- GetCompressionLevel()
  parts <- list(indices = x@i, indptr = x@p, data = x@x)
  for (part in names(parts)) {
    v <- parts[[part]]
    if (level > 0L && length(v) > 0L) {
      grp$create_dataset(name = part, robj = v, dtype = GuessDType(x = v),
                         gzip_level = level, chunk_dims = min(length(v), 65536L))
    } else {
      grp$create_dataset(name = part, robj = v, dtype = GuessDType(x = v))
    }
  }
  .h5_attr(grp, 'dims', dim(x))
}

#' @rdname WriteH5Group
setMethod('WriteH5Group', signature = c(x = 'dgCMatrix'), definition = SparseWrite)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'Graph'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    SparseWrite(x = x, name = name, hgroup = hgroup, verbose = verbose)
    assay <- DefaultAssay(object = x)
    if (!is.null(assay)) .h5_attr(hgroup[[name]], 'assay.used', assay)
    invisible(NULL)
  }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Seurat classes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generic S4 object: one member per slot
.write_s4_slots <- function(x, name, hgroup, verbose = TRUE) {
  grp <- hgroup$create_group(name = name)
  cls <- attr(S4ToList(object = x), 'classDef')
  .h5_attr(grp, 's4class', cls)
  for (s in slotNames(x)) {
    WriteH5Group(x = slot(x, s), name = s, hgroup = grp, verbose = verbose)
  }
  invisible(NULL)
}

# Spatial images: assay / global / class as attributes, other slots as members
ImageWrite <- function(x, name, hgroup, verbose = TRUE) {
  if (!inherits(x, 'SpatialImage')) {
    stop("ImageWrite() needs a SpatialImage-derived object", call. = FALSE)
  }
  grp <- hgroup$create_group(name = name)
  assay <- DefaultAssay(object = x)
  global <- IsGlobal(object = x)
  cls <- GetClass(class = class(x)[1])
  .h5_attr(grp, 'assay', assay)
  .h5_attr(grp, 'global', BoolToInt(x = global), dtype = GuessDType(x = global))
  .h5_attr(grp, 's4class', cls)
  for (s in setdiff(slotNames(x), c('assay', 'global'))) {
    value <- tryCatch(slot(x, s), error = function(e) {
      # images restored from old files may lack a key slot
      if (identical(s, 'key')) paste0(name, '_') else NULL
    })
    if (is.null(value)) {
      if (verbose) message("Skipping slot '", s, "' of image '", name, "' (not accessible)")
      next
    }
    WriteH5Group(x = value, name = s, hgroup = grp, verbose = verbose)
  }
  invisible(NULL)
}

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'ANY'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    # Assay5 / SpatialImage subclasses can reach here when S4 dispatch does
    # not see the specific methods (e.g. Seurat not attached)
    if (inherits(x, 'Assay5')) return(WriteH5GroupAssay5(x, name, hgroup, verbose))
    if (inherits(x, 'SpatialImage')) return(ImageWrite(x, name, hgroup, verbose))
    if (isS4(x)) return(.write_s4_slots(x, name, hgroup, verbose))
    if (!is.null(x)) BasicWrite(x = x, name = name, hgroup = hgroup)
    invisible(NULL)
  }
)

# Key, variable features and (optionally) feature metadata shared by both
# assay writers
.write_assay_common <- function(x, name, grp, features, verbose) {
  if (!is.null(features)) {
    WriteH5Group(x = features, name = 'features', hgroup = grp, verbose = verbose)
  }
  .h5_attr(grp, 'key', Key(object = x))
  hvg <- VariableFeatures(object = x)
  if (length(hvg)) {
    if (verbose) message("Adding variable features for ", name)
    WriteH5Group(x = hvg, name = 'variable.features', hgroup = grp, verbose = verbose)
  } else if (verbose) {
    message("No variable features found for ", name)
  }
  invisible(NULL)
}

# Seurat v5 assays: every layer (including split and feature-subset layers)
# under layers/, with the dimnames of partial layers recorded alongside
WriteH5GroupAssay5 <- function(x, name, hgroup, verbose = TRUE) {
  grp <- hgroup$create_group(name = name)
  layers <- if (inherits(x, 'Assay5')) {
    tryCatch(SeuratObject::Layers(object = x), error = function(e) .STANDARD_LAYERS)
  } else {
    .STANDARD_LAYERS
  }
  cells <- tryCatch(colnames(x), error = function(e) NULL)
  features <- tryCatch(SeuratObject::Features(x), error = function(e) rownames(x))
  for (lyr in layers) {
    mat <- tryCatch(GetAssayData(object = x, layer = lyr), error = function(e) NULL)
    mat <- ConvertBPCellsMatrix(mat, verbose = verbose)
    if (is.null(mat) || IsMatrixEmpty(x = mat)) next
    if (verbose) message("Adding ", lyr, " for ", name)
    layer_features <- features %||% tryCatch(rownames(mat), error = function(e) rownames(x))
    WriteV5Layer(h5_group = grp, layer_name = lyr, matrix_data = mat,
                 features = layer_features, verbose = verbose)
    .WriteLayerDimnames(grp, lyr, mat, cells, layer_features)
  }
  .WriteDefaultLayer(grp, x)
  .write_assay_common(x, name, grp, features, verbose)
  .h5_attr(grp, 's4class', 'SeuratObject::Assay5',
           dtype = CachedGuessDType(x = 'SeuratObject::Assay5'))
  invisible(NULL)
}

# Record which cells / features a layer covers when it does not span the
# whole assay (split layers such as "counts.sample1", feature-subset
# layers), plus the default layer. Stored under `layer.cells/<layer>`,
# `layer.features/<layer>` and the `default.layer` attribute so that
# AssembleAssay() can rebuild the Assay5 exactly; readers that predate these
# fields simply ignore them.
.WriteLayerDimnames <- function(xgroup, layer, dat, all_cells, all_features) {
  lc <- colnames(x = dat); lf <- rownames(x = dat)
  if (identical(layer, 'scale.data') && !is.null(lf) && !xgroup$exists('scaled.features')) {
    # same convention as the v3 Assay writer; the index requires it before
    # scale.data is considered loadable
    xgroup$create_dataset(name = 'scaled.features', robj = lf, dtype = GuessDType(x = lf))
  }
  if (!is.null(lc) && !is.null(all_cells) && !identical(lc, all_cells)) {
    if (!xgroup$exists('layer.cells')) xgroup$create_group('layer.cells')
    xgroup[['layer.cells']]$create_dataset(name = layer, robj = lc, dtype = GuessDType(x = lc))
  }
  if (!is.null(lf) && !is.null(all_features) && !identical(lf, all_features) &&
      !identical(layer, 'scale.data')) {
    if (!xgroup$exists('layer.features')) xgroup$create_group('layer.features')
    xgroup[['layer.features']]$create_dataset(name = layer, robj = lf, dtype = GuessDType(x = lf))
  }
  invisible(NULL)
}

.WriteDefaultLayer <- function(xgroup, x) {
  dl <- tryCatch(SeuratObject::DefaultLayer(object = x), error = function(e) NULL)
  if (length(dl)) {
    xgroup$create_attr(attr_name = 'default.layer', robj = dl, dtype = GuessDType(x = dl))
  }
  lo <- tryCatch(SeuratObject::Layers(object = x), error = function(e) NULL)
  if (length(lo)) {
    xgroup$create_attr(attr_name = 'layer.order', robj = lo, dtype = GuessDType(x = lo))
  }
  invisible(NULL)
}

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'Assay5'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    WriteH5GroupAssay5(x = x, name = name, hgroup = hgroup, verbose = verbose)
  }
)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'Assay'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    if (inherits(x, 'Assay5')) return(WriteH5GroupAssay5(x, name, hgroup, verbose))
    grp <- hgroup$create_group(name = name)
    # Seurat v3/v4 assays: counts, data and scale.data directly in the group
    for (lyr in .STANDARD_LAYERS) {
      mat <- ConvertBPCellsMatrix(GetAssayDataCompat(object = x, layer_or_slot = lyr),
                                  verbose = verbose)
      if (is.null(mat) || IsMatrixEmpty(x = mat)) next
      if (verbose) message("Adding ", lyr, " for ", name)
      WriteH5Group(x = mat, name = lyr, hgroup = grp, verbose = verbose)
      if (lyr == 'scale.data') {
        WriteH5Group(x = rownames(mat), name = 'scaled.features', hgroup = grp, verbose = verbose)
      }
    }
    .write_assay_common(x, name, grp, rownames(x), verbose)
    feature_meta <- x[[]]
    if (ncol(feature_meta)) {
      if (verbose) message("Adding feature-level metadata for ", name)
      WriteH5Group(x = feature_meta, name = 'meta.features', hgroup = grp, verbose = verbose)
    } else if (verbose) {
      message("No feature-level metadata found for ", name)
    }
    WriteH5Group(x = Misc(object = x), name = 'misc', hgroup = grp, verbose = verbose)
    # subclasses (e.g. SCTAssay) keep their class and their extra slots
    if (!identical(class(x)[1], 'Assay')) {
      .h5_attr(grp, 's4class', attr(S4ToList(object = x), 'classDef'))
      extra <- setdiff(slotNames(x), slotNames(methods::getClassDef(Class = 'Assay')))
      for (s in extra) {
        if (verbose) message("Writing out ", s, " for ", name)
        WriteH5Group(x = slot(x, s), name = s, hgroup = grp, verbose = verbose)
      }
    }
    invisible(NULL)
  }
)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'DimReduc'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    grp <- hgroup$create_group(name = name)
    if (verbose) message("Adding cell embeddings for ", name)
    WriteH5Group(x = Embeddings(object = x), name = 'cell.embeddings', hgroup = grp,
                 verbose = verbose)
    for (projected in c(FALSE, TRUE)) {
      what <- if (projected) 'projected loadings' else 'loadings'
      ld <- Loadings(object = x, projected = projected)
      if (IsMatrixEmpty(x = ld)) {
        if (verbose) message("No ", what, " for ", name)
        next
      }
      if (verbose) message("Adding ", what, " for ", name)
      WriteH5Group(x = ld, hgroup = grp, verbose = verbose,
                   name = if (projected) 'feature.loadings.projected' else 'feature.loadings')
      WriteH5Group(x = rownames(ld), hgroup = grp, verbose = verbose,
                   name = if (projected) 'projected.features' else 'features')
    }
    global <- IsGlobal(object = x)
    .h5_attr(grp, 'active.assay', DefaultAssay(object = x))
    .h5_attr(grp, 'key', Key(object = x))
    .h5_attr(grp, 'global', BoolToInt(x = global), dtype = GuessDType(x = global))
    sdev <- Stdev(object = x)
    if (length(sdev)) {
      if (verbose) message("Adding standard deviations for ", name)
      WriteH5Group(x = sdev, name = 'stdev', hgroup = grp, verbose = verbose)
    } else if (verbose) {
      message("No standard deviations for ", name)
    }
    js <- JS(object = x)
    if (as.logical(x = js)) {
      if (verbose) message("Adding JackStraw information for ", name)
      WriteH5Group(x = js, name = 'jackstraw', hgroup = grp, verbose = verbose)
    } else if (verbose) {
      message("No JackStraw data for ", name)
    }
    WriteH5Group(x = slot(x, 'misc'), name = 'misc', hgroup = grp, verbose = verbose)
    invisible(NULL)
  }
)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'Neighbor'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    grp <- hgroup$create_group(name = name)
    for (s in slotNames(x)) {
      value <- slot(x, s)
      if (identical(s, 'alg.idx') && !is.null(value)) {
        # the Annoy/RcppHNSW index is an external pointer
        warning("Neighbor search indexes cannot be stored; save the index separately",
                call. = FALSE, immediate. = TRUE)
        next
      }
      WriteH5Group(x = value, name = s, hgroup = grp, verbose = verbose)
    }
    invisible(NULL)
  }
)

#' @rdname WriteH5Group
setMethod(
  'WriteH5Group',
  signature = c(x = 'SeuratCommand'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    # parameters become the group's members, everything else an attribute
    params <- Filter(Negate(is.function), as.list(x))
    WriteH5Group(x = params, name = name, hgroup = hgroup, verbose = verbose)
    for (s in setdiff(slotNames(x), 'params')) {
      value <- slot(x, s)
      if (!is.null(value)) .h5_attr(hgroup[[name]], s, as.character(value))
    }
    invisible(NULL)
  }
)
