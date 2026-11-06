#' @include zzz.R
#' @include H5Utils.R
#' @importFrom hdf5r H5File h5attr h5attr_names H5S h5types
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Streaming (out-of-memory) conversion between HDF5 single-cell formats
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Expression matrices are the only parts of a single-cell file that do not
# fit in memory; cell/feature metadata, embeddings, graphs, images and
# unstructured data are small. A streaming conversion (h5ad, h5mu,
# h5Seurat, loom, any pair):
#
#   1. writes a *skeleton* of the source: an HDF5-level copy in which every
#      expression matrix is replaced by a placeholder sparse matrix of the
#      same shape holding a single entry whose value is a unique tag;
#   2. reads the skeleton with the regular reader, checks that each
#      placeholder kept the source's cell and feature order, computes
#      nCount_/nFeature_ from the real counts, and writes the object with
#      the regular writer, so all metadata is handled exactly as in an
#      in-memory conversion;
#   3. finds every tag in the output and writes the real matrix there:
#      a plain HDF5 copy when the layout is unchanged (CSR of cells x genes
#      has the same arrays as CSC of genes x cells), otherwise in blocks of
#      cells sized by options(scConvert.stream_chunk_bytes).
#
# If a reader or writer rearranges a matrix in a way the engine cannot
# reproduce, the conversion falls back to the in-memory path.

.STREAM_TAG_BASE <- 9.5e6          # exactly representable in float32 too
.stream_tag <- function(k) .STREAM_TAG_BASE + k
.is_tag <- function(v) v >= .STREAM_TAG_BASE & v < .STREAM_TAG_BASE + 1e5 & v == round(v)

# Resolve the `stream` argument: TRUE/FALSE, or "auto" (stream files larger
# than options(scConvert.stream_threshold), 2 GB by default)
.use_streaming <- function(stream, path) {
  if (is.null(stream)) stream <- getOption('scConvert.stream', 'auto')
  if (isTRUE(stream) || isFALSE(stream)) return(stream)
  if (!identical(stream, 'auto')) stop("'stream' must be TRUE, FALSE or \"auto\"", call. = FALSE)
  size <- suppressWarnings(file.size(path))
  isTRUE(size > getOption('scConvert.stream_threshold', 2e9))
}

# Elements per streamed block (~12 bytes per non-zero: value + index)
.stream_block <- function() {
  max(1e5, floor(getOption('scConvert.stream_chunk_bytes', 64 * 1024^2) / 12))
}

.stream_fallback <- function(msg) {
  structure(class = c('scConvert_stream_fallback', 'error', 'condition'),
            list(message = msg, call = NULL))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# HDF5 helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.is_sparse_group <- function(obj) {
  inherits(obj, 'H5Group') && all(c('data', 'indices', 'indptr') %in% names(obj))
}

.attr_or <- function(obj, name, default = NULL) {
  if (obj$attr_exists(attr_name = name)) h5attr(obj, name) else default
}

# Copy every attribute of `from` onto `to` (same type and dataspace)
.copy_attrs <- function(from, to, skip = character(0)) {
  for (nm in setdiff(h5attr_names(from), skip)) {
    a <- from$attr_open(attr_name = nm)
    if (to$attr_exists(attr_name = nm)) to$attr_delete(attr_name = nm)
    to$create_attr(attr_name = nm, robj = a$read(), dtype = a$get_type(), space = a$get_space())
    a$close()
  }
  invisible(NULL)
}

# Copy the tree under `from` into `to`, writing placeholders at `replace`
# (named list: in-file path -> function(parent, name) creating it)
.copy_tree <- function(from, to, replace, prefix = '') {
  .copy_attrs(from, to)
  for (nm in names(from)) {
    path <- if (nzchar(prefix)) paste0(prefix, '/', nm) else nm
    if (!is.null(replace[[path]])) {
      replace[[path]](to, nm)
    } else if (any(startsWith(names(replace), paste0(path, '/')))) {
      .copy_tree(from[[nm]], to$create_group(nm), replace, path)
    } else {
      to$obj_copy_from(src_loc = from, src_name = nm, dst_name = nm)
    }
  }
  invisible(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Source matrices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# A source matrix is described by
#   path       location in the source file
#   kind       "cells"     compressed sparse, one compressed slice per cell
#                          (h5ad/h5mu CSR, h5Seurat CSC)
#              "features"  compressed sparse, one slice per feature (h5ad CSC)
#              "dense"     dense dataset; `cell_axis` is the cell axis in
#                          hdf5r's (R) dimension order
#   features, cells        names in source order
#   row_map    0-based output row for each source feature (NA = dropped);
#              NULL keeps all features in order
#   col_keep   logical per source cell (NULL = all)

.h5ad_df_index <- function(df) {
  key <- .attr_or(df, '_index', '_index')
  if (!df$exists(key)) key <- intersect(c('_index', 'index'), names(df))[1]
  as.character(df[[key]][])
}

.anndata_kind <- function(obj) {
  if (inherits(obj, 'H5D')) return('dense')
  if (identical(.attr_or(obj, 'encoding-type', 'csr_matrix'), 'csc_matrix')) 'features' else 'cells'
}

# Matrices of one AnnData group (an h5ad file, or one h5mu modality)
.anndata_matrices <- function(grp, prefix) {
  genes <- .h5ad_df_index(grp[['var']])
  cells <- .h5ad_df_index(grp[['obs']])
  p <- function(x) if (nzchar(prefix)) paste0(prefix, '/', x) else x
  out <- list()
  add <- function(path, obj, feats, row_map = NULL) {
    out[[length(out) + 1]] <<- list(path = p(path), kind = .anndata_kind(obj), cell_axis = 2L,
                                    features = feats, cells = cells, row_map = row_map,
                                    fmt = 'anndata')
  }
  if (grp$exists('X')) add('X', grp[['X']], genes)
  if (grp$exists('layers')) {
    for (nm in names(grp[['layers']])) add(paste0('layers/', nm), grp[['layers']][[nm]], genes)
  }
  if (grp$exists('raw') && grp[['raw']]$exists('X')) {
    raw_genes <- if (grp[['raw']]$exists('var')) .h5ad_df_index(grp[['raw/var']]) else genes
    # the readers keep raw features shared with X, in X order (matched
    # after their "_" -> "-" renaming)
    rm <- if (identical(raw_genes, genes)) NULL else
      match(gsub('_', '-', raw_genes, fixed = TRUE), gsub('_', '-', genes, fixed = TRUE)) - 1L
    add('raw/X', grp[['raw/X']], raw_genes, rm)
  }
  out
}

.source_matrices <- function(src, fmt) {
  switch(
    fmt,
    h5ad = .anndata_matrices(src, ''),
    h5mu = unlist(lapply(names(src[['mod']]), function(m) {
      .anndata_matrices(src[['mod']][[m]], paste0('mod/', m))
    }), recursive = FALSE),
    h5seurat = .h5seurat_matrices(src),
    loom = .loom_matrices(src),
    stop("Streaming is not supported for ", fmt, " files", call. = FALSE)
  )
}

.h5seurat_matrices <- function(src) {
  cells <- as.character(Cells.h5Seurat(src))
  out <- list()
  for (a in names(src[['assays']])) {
    ag <- src[['assays']][[a]]
    feats <- if (ag$exists('features')) as.character(ag[['features']][]) else character(0)
    sub_names <- function(group, layer, default) {
      if (ag$exists(group) && ag[[group]]$exists(layer)) as.character(ag[[group]][[layer]][]) else default
    }
    add <- function(path, obj, layer) {
      f <- if (identical(layer, 'scale.data') && ag$exists('scaled.features')) {
        as.character(ag[['scaled.features']][])
      } else sub_names('layer.features', layer, feats)
      out[[length(out) + 1]] <<- list(
        path = paste0('assays/', a, '/', path),
        kind = if (inherits(obj, 'H5D')) 'dense' else 'cells', cell_axis = 2L,
        features = f, cells = sub_names('layer.cells', layer, cells), row_map = NULL,
        fmt = 'h5seurat')
    }
    for (slot in intersect(c('counts', 'data', 'scale.data'), names(ag))) add(slot, ag[[slot]], slot)
    if (ag$exists('layers') && inherits(ag[['layers']], 'H5Group')) {
      for (lyr in names(ag[['layers']])) {
        obj <- ag[['layers']][[lyr]]
        if (inherits(obj, 'H5D') || .is_sparse_group(obj)) add(paste0('layers/', lyr), obj, lyr)
      }
    }
  }
  out
}

.loom_matrices <- function(src) {
  genes <- as.character(src[['row_attrs/Gene']][])
  cells <- as.character(src[['col_attrs/CellID']][])
  # readLoom()'s default filter: cells flagged in col_attrs/Valid
  keep <- if (Exists(src, 'col_attrs/Valid')) as.logical(src[['col_attrs/Valid']][]) else NULL
  paths <- c('matrix', if (src$exists('layers')) paste0('layers/', names(src[['layers']])))
  lapply(paths, function(p) {
    list(path = p, kind = 'dense', cell_axis = 1L, features = genes, cells = cells,
         row_map = NULL, col_keep = keep, fmt = 'loom')
  })
}

# Output (row, column) positions, 1-based, of source entry (feature f, cell c)
.out_position <- function(d, f, c) {
  row <- if (is.null(d$row_map)) f else d$row_map[f] + 1L
  col <- if (is.null(d$col_keep)) c else sum(d$col_keep[seq_len(c)])
  c(row, col)
}

.n_out <- function(d) {
  c(if (is.null(d$row_map)) length(d$features) else sum(!is.na(d$row_map)),
    if (is.null(d$col_keep)) length(d$cells) else sum(d$col_keep))
}

# Placeholder written in place of a source matrix, in the source format's
# own sparse layout: one entry (tag) at the first kept feature / cell
.placeholder_writer <- function(d, tag, template) {
  force(tag)
  force(template)
  f0 <- if (is.null(d$row_map)) 1L else which(!is.na(d$row_map))[1]
  c0 <- if (is.null(d$col_keep)) 1L else which(d$col_keep)[1]
  n_cells <- length(d$cells)
  n_feat <- length(d$features)
  function(parent, name) {
    g <- parent$create_group(name)
    g$create_dataset('data', robj = tag, dtype = h5types$H5T_NATIVE_DOUBLE)
    g$create_dataset('indices', robj = as.integer(f0 - 1L), dtype = h5types$H5T_NATIVE_INT32)
    g$create_dataset('indptr', robj = c(rep(0L, c0), rep(1L, n_cells - c0 + 1L)),
                     dtype = h5types$H5T_NATIVE_INT32)
    sc <- function(nm, v) g$create_attr(nm, v, dtype = CachedGuessDType(v), space = ScalarSpace())
    if (d$fmt == 'anndata') {
      sc('encoding-type', 'csr_matrix'); sc('encoding-version', '0.1.0')
      g$create_attr('shape', as.integer(c(n_cells, n_feat)), dtype = GuessDType(1L))
    } else {
      # h5Seurat (v5 layers: shape + csc encoding; v3 slots: dims) and the
      # loom reader's sparse stand-in: CSC features x cells
      if (inherits(template, 'H5Group')) .copy_attrs(template, g)
      for (a in c('shape', 'dims')) {
        if (g$attr_exists(attr_name = a)) g$attr_delete(attr_name = a)
        g$create_attr(a, as.integer(c(n_feat, n_cells)), dtype = GuessDType(1L))
      }
    }
    invisible(g)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reading source matrices in blocks of cells
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Consecutive index ranges (1-based, inclusive) of at most `block` entries
# according to the cumulative pointer `ptr` (a single larger item is its own
# range)
.major_blocks <- function(ptr, block) {
  n <- length(ptr) - 1
  if (n < 1) return(matrix(integer(0), ncol = 2))
  starts <- integer(0)
  ends <- integer(0)
  s <- 1
  while (s <= n) {
    e <- max(s, min(findInterval(ptr[s] + block, ptr) - 1, n))
    starts <- c(starts, s)
    ends <- c(ends, e)
    s <- e + 1
  }
  cbind(starts, ends)
}

# Call fun(block, out_cols) for consecutive blocks of the source matrix,
# where `block` is a features_out x cells dgCMatrix already reduced to the
# kept features (in output order) and kept cells, and `out_cols` are the
# output column numbers of those cells
.for_blocks <- function(src_file, d, fun) {
  obj <- src_file[[d$path]]
  n_feat <- length(d$features)
  n_cells <- length(d$cells)
  sel <- if (is.null(d$row_map)) NULL else match(seq_len(sum(!is.na(d$row_map))) - 1L, d$row_map)
  out_col <- if (is.null(d$col_keep)) seq_len(n_cells) else cumsum(d$col_keep)
  emit <- function(m, cells) {
    keep <- if (is.null(d$col_keep)) rep(TRUE, length(cells)) else d$col_keep[cells]
    if (!is.null(sel)) m <- m[sel, , drop = FALSE]
    if (!all(keep)) m <- m[, keep, drop = FALSE]
    if (ncol(m)) fun(methods::as(m, 'CsparseMatrix'), out_col[cells][keep])
  }
  if (d$kind == 'cells') {
    ptr <- as.numeric(obj[['indptr']][])
    spans <- .major_blocks(ptr, .stream_block())
    for (b in seq_len(nrow(spans))) {
      cells <- spans[b, 1]:spans[b, 2]
      a <- ptr[spans[b, 1]] + 1
      z <- ptr[spans[b, 2] + 1]
      rng <- if (z >= a) seq.int(a, z) else integer(0)
      m <- Matrix::sparseMatrix(
        i = if (length(rng)) obj[['indices']][rng] + 1 else integer(0),
        p = ptr[spans[b, 1]:(spans[b, 2] + 1)] - ptr[spans[b, 1]],
        x = if (length(rng)) as.numeric(obj[['data']][rng]) else numeric(0),
        dims = c(n_feat, length(cells)))
      emit(m, cells)
    }
  } else if (d$kind == 'dense') {
    per <- max(1L, floor(getOption('scConvert.stream_chunk_bytes', 64 * 1024^2) / (8 * n_feat)))
    for (s in seq.int(1L, n_cells, by = per)) {
      cells <- s:min(n_cells, s + per - 1L)
      v <- if (d$cell_axis == 1L) t(obj[cells, , drop = FALSE]) else obj[, cells, drop = FALSE]
      emit(methods::as(Matrix::Matrix(v, sparse = TRUE), 'CsparseMatrix'), cells)
    }
  } else {
    # one compressed slice per feature: a transpose needs the whole matrix
    ptr <- obj[['indptr']][]
    m <- Matrix::sparseMatrix(i = obj[['indices']][] + 1, p = ptr, x = as.numeric(obj[['data']][]),
                              dims = c(n_cells, length(ptr) - 1))
    emit(Matrix::t(m), seq_len(n_cells))
  }
  invisible(NULL)
}

# nCount / nFeature per output cell
.block_stats <- function(src_file, d) {
  n <- .n_out(d)[2]
  sums <- numeric(n)
  nfeat <- integer(n)
  .for_blocks(src_file, d, function(m, cols) {
    sums[cols] <<- Matrix::colSums(m)
    nfeat[cols] <<- as.integer(Matrix::colSums(m > 0))
  })
  list(sum = sums, n = nfeat)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Writing into the output
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Replace the arrays of placeholder group `dst` (one compressed slice per
# cell) with the source matrix
.fill_cells_major <- function(src_file, d, dst) {
  for (nm in c('data', 'indices', 'indptr')) dst$link_delete(nm)
  src <- src_file[[d$path]]
  if (d$kind == 'cells' && is.null(d$row_map) && is.null(d$col_keep)) {
    # same layout: let HDF5 copy the three datasets
    for (nm in c('data', 'indices', 'indptr')) {
      dst$obj_copy_from(src_loc = src, src_name = nm, dst_name = nm)
    }
    return(invisible(NULL))
  }
  level <- GetCompressionLevel()
  grow <- function(name, dtype) {
    dst$create_dataset(name = name, dtype = dtype, space = H5S$new(dims = 0, maxdims = Inf),
                       chunk_dims = 65536, gzip_level = level)
  }
  d_x <- grow('data', h5types$H5T_NATIVE_DOUBLE)
  d_i <- grow('indices', h5types$H5T_NATIVE_INT32)
  counts <- numeric(.n_out(d)[2])
  pos <- 0
  .for_blocks(src_file, d, function(m, cols) {
    n <- length(m@x)
    if (n) {
      d_x$set_extent(pos + n)
      d_i$set_extent(pos + n)
      d_x[(pos + 1):(pos + n)] <- m@x
      d_i[(pos + 1):(pos + n)] <- m@i
      pos <<- pos + n
    }
    counts[cols] <<- diff(m@p)
  })
  ptr <- c(0, cumsum(counts))
  dst$create_dataset('indptr', robj = ptr, gzip_level = level,
                     dtype = if (pos > .Machine$integer.max) h5types$H5T_STD_I64LE else h5types$H5T_NATIVE_INT32,
                     chunk_dims = max(1, min(length(ptr), 65536)))
  invisible(NULL)
}

# Placeholder group stored with one slice per feature: transpose in memory
.fill_features_major <- function(src_file, d, dst, verbose) {
  if (verbose) message("  ", d$path, ": the destination stores it by feature; transposing in memory")
  parts <- list()
  .for_blocks(src_file, d, function(m, cols) parts[[length(parts) + 1]] <<- list(m = m, cols = cols))
  full <- do.call(cbind, lapply(parts, `[[`, 'm'))
  full <- full[, order(unlist(lapply(parts, `[[`, 'cols'))), drop = FALSE]
  tm <- methods::as(Matrix::t(full), 'CsparseMatrix')
  for (nm in c('data', 'indices', 'indptr')) dst$link_delete(nm)
  level <- GetCompressionLevel()
  for (part in list(list('data', tm@x), list('indices', tm@i), list('indptr', tm@p))) {
    dst$create_dataset(part[[1]], robj = part[[2]], gzip_level = level,
                       chunk_dims = max(1, min(length(part[[2]]), 65536)))
  }
  invisible(NULL)
}

# Dense dataset in the output (loom): written block by block
.fill_dense <- function(src_file, d, ds) {
  dims <- ds$dims
  n_out <- .n_out(d)
  cells_first <- identical(as.numeric(dims), as.numeric(rev(n_out)))
  .for_blocks(src_file, d, function(m, cols) {
    v <- as.matrix(m)
    if (cells_first) ds[cols, ] <- t(v) else ds[, cols] <- v
  })
  invisible(NULL)
}

# Every tag in the output: data.frame(path, tag, type) where type is
# "sparse" (one-entry compressed group) or "dense" (dataset whose entry at
# the tag's position holds the tag)
.find_output_tags <- function(out, descs, prefix = '') {
  found <- data.frame(path = character(0), tag = numeric(0), type = character(0),
                      stringsAsFactors = FALSE)
  for (nm in names(out)) {
    path <- if (nzchar(prefix)) paste0(prefix, '/', nm) else nm
    obj <- tryCatch(out[[nm]], error = function(e) NULL)
    if (inherits(obj, 'H5Group')) {
      if (.is_sparse_group(obj) && obj[['data']]$dims <= 2) {
        v <- as.numeric(obj[['data']][])
        if (length(v) && .is_tag(v[1])) {
          found[nrow(found) + 1, ] <- list(path, v[1], 'sparse')
          next
        }
      }
      found <- rbind(found, .find_output_tags(obj, descs, path))
    } else if (inherits(obj, 'H5D') && length(obj$dims) == 2 &&
               inherits(obj$get_type(), c('H5T_FLOAT', 'H5T_INTEGER'))) {
      for (tg in names(descs)) {
        pos <- descs[[tg]]$out_pos
        n_out <- .n_out(descs[[tg]])
        at <- if (identical(as.numeric(obj$dims), as.numeric(rev(n_out)))) rev(pos) else
          if (identical(as.numeric(obj$dims), as.numeric(n_out))) pos else NULL
        if (!is.null(at) && isTRUE(obj[at[1], at[2]] == as.numeric(tg))) {
          found[nrow(found) + 1, ] <- list(path, as.numeric(tg), 'dense')
          break
        }
      }
    }
  }
  found
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Checks on the skeleton object
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Layers of `obj` that hold a placeholder: list of (assay, layer, tag, row, col)
.object_tags <- function(obj) {
  out <- list()
  for (a in SeuratObject::Assays(obj)) {
    layers <- tryCatch(SeuratObject::Layers(obj[[a]]), error = function(e) character(0))
    for (l in layers) {
      m <- tryCatch(suppressWarnings(SeuratObject::LayerData(obj, assay = a, layer = l)),
                    error = function(e) NULL)
      if (is.null(m) || !length(m)) next
      if (inherits(m, 'dgCMatrix')) {
        hit <- which(.is_tag(m@x))
        if (!length(hit) || length(m@x) > 2) next
        row <- m@i[hit[1]] + 1
        col <- findInterval(hit[1] - 1, m@p, rightmost.closed = FALSE)
      } else if (is.matrix(m)) {
        w <- which(.is_tag(m), arr.ind = TRUE)
        if (!nrow(w)) next
        row <- w[1, 1]
        col <- w[1, 2]
      } else next
      v <- if (inherits(m, 'dgCMatrix')) m@x[hit[1]] else m[row, col]
      out[[length(out) + 1]] <- list(assay = a, layer = l, tag = v, row = row, col = col,
                                     dim = dim(m), dimnames = dimnames(m))
    }
  }
  out
}

.same_names <- function(found, expected) {
  if (is.null(found)) return(TRUE)
  norm <- function(x) make.unique(gsub('_', '-', as.character(x), fixed = TRUE))
  identical(as.character(found), as.character(expected)) || identical(as.character(found), norm(expected))
}

# Each placeholder must still be the whole source matrix, with the
# source's features and cells in the same order
.check_object <- function(tags, descs) {
  for (t in tags) {
    d <- descs[[as.character(t$tag)]]
    if (is.null(d)) next
    feats <- if (is.null(d$row_map)) d$features else d$features[order(d$row_map, na.last = NA)]
    cells <- if (is.null(d$col_keep)) d$cells else d$cells[d$col_keep]
    ok <- identical(as.numeric(t$dim), as.numeric(.n_out(d))) &&
      identical(as.numeric(c(t$row, t$col)), as.numeric(d$out_pos)) &&
      .same_names(t$dimnames[[1]], feats) && .same_names(t$dimnames[[2]], cells)
    if (!ok) {
      stop(.stream_fallback(sprintf("the reader rearranged %s (%s/%s)", d$path, t$assay, t$layer)))
    }
  }
  invisible(TRUE)
}

# Readers compute nCount_/nFeature_ from the counts they loaded, which here
# are placeholders; recompute them from the real matrix
.fix_count_stats <- function(obj, tags, descs, src) {
  md <- obj@meta.data
  for (t in Filter(function(t) identical(t$layer, 'counts'), tags)) {
    d <- descs[[as.character(t$tag)]]
    cols <- paste0(c('nCount_', 'nFeature_'), t$assay)
    looks_placeholder <- vapply(cols, function(cl) {
      v <- md[[cl]]
      !is.null(v) && sum(v != 0, na.rm = TRUE) <= 1 && max(v, na.rm = TRUE) <= .STREAM_TAG_BASE + 1e5
    }, logical(1))
    if (!any(looks_placeholder) || is.null(d)) next
    stats <- .block_stats(src, d)
    if (looks_placeholder[1]) md[[cols[1]]] <- stats$sum
    if (looks_placeholder[2]) md[[cols[2]]] <- stats$n
  }
  obj@meta.data <- md
  obj
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Engine
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SeuratObject Assays Layers LayerData
#' @keywords internal
#' @noRd
.stream_convert <- function(source, dest, stype, dtype, assay = NULL, overwrite = FALSE,
                            verbose = TRUE, standardize = FALSE) {
  if (file.exists(dest) && !overwrite) {
    stop("Destination file '", dest, "' already exists; set overwrite = TRUE", call. = FALSE)
  }
  res <- tryCatch(
    .stream_convert_impl(source, dest, stype, dtype, assay, verbose, standardize),
    scConvert_stream_fallback = function(e) e
  )
  if (inherits(res, 'scConvert_stream_fallback')) {
    if (verbose) message("Streaming not possible (", conditionMessage(res), "); converting in memory")
    if (file.exists(dest)) unlink(dest, recursive = TRUE)
    extra <- if (dtype == 'h5ad') list(standardize = standardize) else list()
    return(do.call(.hub_from_connection, c(list(source, dest, stype = stype, dtype = dtype,
                                                 assay = assay %||% 'RNA', overwrite = TRUE,
                                                 verbose = verbose), extra)))
  }
  invisible(dest)
}

.stream_convert_impl <- function(source, dest, stype, dtype, assay, verbose, standardize) {
  src <- H5File$new(source, mode = 'r')
  on.exit(tryCatch(src$close_all(), error = function(e) NULL), add = TRUE)
  descs <- .source_matrices(src, stype)
  # h5ad: raw/X lacking some of X's features makes the reader subset X too
  raw <- Filter(function(d) endsWith(d$path, 'raw/X') && !is.null(d$row_map), descs)
  for (d in raw) {
    x_genes <- descs[[which(vapply(descs, function(e) identical(e$path, sub('raw/X$', 'X', d$path)), logical(1)))[1]]]$features
    if (sum(!is.na(d$row_map)) < length(x_genes)) {
      stop(.stream_fallback("raw/X lacks some of X's features"))
    }
  }
  tags <- vapply(seq_along(descs), .stream_tag, numeric(1))
  names(descs) <- as.character(tags)
  replace <- list()
  for (i in seq_along(descs)) {
    d <- descs[[i]]
    f0 <- if (is.null(d$row_map)) 1L else which(!is.na(d$row_map))[1]
    c0 <- if (is.null(d$col_keep)) 1L else which(d$col_keep)[1]
    descs[[i]]$out_pos <- .out_position(d, f0, c0)
    replace[[d$path]] <- .placeholder_writer(d, tags[i], src[[d$path]])
  }
  if (!length(descs)) stop(.stream_fallback("no expression matrices found"))

  skel_path <- tempfile(fileext = paste0('.', stype))
  on.exit(unlink(skel_path), add = TRUE)
  skel <- H5File$new(skel_path, mode = 'w')
  .copy_tree(src, skel, replace)
  skel$close_all()

  if (verbose) message("Converting metadata, embeddings and graphs (matrices are streamed)")
  obj <- suppressWarnings(GetLoader(ext = stype)(file = skel_path, assay = assay %||% 'RNA',
                                                 verbose = FALSE))
  if (!inherits(obj, 'Seurat')) stop(.stream_fallback("the source does not read as one Seurat object"))
  if (!is.null(obj@misc$scConvert_read)) obj@misc$scConvert_read$source_file <- normalizePath(source)
  if (stype == 'h5seurat' && !is.null(assay) && assay %in% Assays(obj)) {
    SeuratObject::DefaultAssay(obj) <- assay
  }
  obj_tags <- .object_tags(obj)
  .check_object(obj_tags, descs)
  obj <- .fix_count_stats(obj, obj_tags, descs, src)

  saver_args <- list(object = obj, filename = dest, overwrite = TRUE, verbose = FALSE)
  if (dtype == 'h5ad') saver_args$standardize <- standardize
  old <- options(scConvert.stream_placeholders = TRUE)
  on.exit(options(old), add = TRUE)
  do.call(GetSaver(ext = dtype), saver_args)
  options(old)
  rm(obj)

  out <- H5File$new(dest, mode = 'r+')
  on.exit(tryCatch(out$close_all(), error = function(e) NULL), add = TRUE)
  found <- .find_output_tags(out, descs)
  placed <- unique(found$tag)
  needed <- unique(vapply(obj_tags, function(t) t$tag, numeric(1)))
  if (!length(placed) && length(needed)) stop(.stream_fallback("the writer stored no matrix"))
  for (i in seq_len(nrow(found))) {
    d <- descs[[as.character(found$tag[i])]]
    if (verbose) message("Streaming ", d$path, " -> ", found$path[i])
    target <- out[[found$path[i]]]
    if (found$type[i] == 'dense') {
      .fill_dense(src, d, target)
      next
    }
    n_major <- target[['indptr']]$dims - 1
    n_out <- .n_out(d)
    enc <- .attr_or(target, 'encoding-type', '')
    by_feature <- identical(enc, 'csc_matrix') && dtype %in% c('h5ad', 'h5mu')
    if (by_feature && n_major == n_out[1]) {
      .fill_features_major(src, d, target, verbose)
    } else if (!by_feature && n_major == n_out[2]) {
      .fill_cells_major(src, d, target)
    } else {
      stop(.stream_fallback(sprintf("the writer stored %s with an unexpected shape", d$path)))
    }
  }
  out$close_all()
  if (verbose) message("Streaming ", stype, " -> ", dtype, " complete: ", dest)
  invisible(dest)
}

# Kept for callers of the first streaming release
.stream_h5ad_to_h5seurat <- function(source, dest, assay = 'RNA', overwrite = FALSE,
                                     verbose = TRUE, ...) {
  .stream_convert(source, dest, 'h5ad', 'h5seurat', assay = assay, overwrite = overwrite,
                  verbose = verbose)
}

.stream_h5seurat_to_h5ad <- function(source, dest, assay = NULL, overwrite = FALSE,
                                     verbose = TRUE, standardize = FALSE, ...) {
  .stream_convert(source, dest, 'h5seurat', 'h5ad', assay = assay, overwrite = overwrite,
                  verbose = verbose, standardize = standardize)
}
