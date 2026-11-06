#' @include zzz.R
#' @include H5Utils.R
#' @importFrom hdf5r H5File h5attr h5attr_names H5S h5types
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Streaming (out-of-memory) h5ad <-> h5Seurat conversion
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# The expression matrices are the only parts of a single-cell file that do
# not fit in memory; everything else (cell/feature metadata, embeddings,
# graphs, images, uns/misc) is small. Streaming conversion therefore:
#
#   1. writes a *skeleton* copy of the source in which every expression
#      matrix is replaced by a placeholder sparse matrix of the same shape
#      holding one entry whose value is a unique tag; all other members are
#      copied at the HDF5 level (no R round trip);
#   2. converts the skeleton with the regular reader and writer, so all
#      metadata is handled by exactly the same code as an in-memory
#      conversion;
#   3. finds every tagged placeholder in the output (wherever the writer put
#      it) and replaces it with the real matrix, copied from the source.
#
# CSR of a cells x genes matrix has the same arrays as CSC of its genes x
# cells transpose, so h5ad X <-> h5Seurat layers is a plain copy of the
# data/indices/indptr datasets, done by HDF5 itself. Only when rows have to
# be subset or reordered (raw/X with its own feature list) is the matrix
# processed in blocks of columns sized by options(scConvert.stream_chunk_bytes).
# Peak memory is a few blocks plus the small components, independent of the
# number of non-zeros.

# Tag values: large, distinctive and exactly representable in float32
.STREAM_TAG_BASE <- 9.5e6

.stream_tag <- function(k) .STREAM_TAG_BASE + k

# Resolve the `stream` argument: TRUE/FALSE, or "auto" (stream files larger
# than options(scConvert.stream_threshold), 2 GB by default)
.use_streaming <- function(stream, path) {
  if (is.null(stream)) stream <- getOption('scConvert.stream', 'auto')
  if (isTRUE(stream) || isFALSE(stream)) return(stream)
  if (!identical(stream, 'auto')) stop("'stream' must be TRUE, FALSE or \"auto\"", call. = FALSE)
  size <- suppressWarnings(file.size(path))
  isTRUE(size > getOption('scConvert.stream_threshold', 2e9))
}

# Elements per streamed block (data + index ~ 12 bytes per non-zero)
.stream_block <- function() {
  budget <- getOption('scConvert.stream_chunk_bytes', 64 * 1024^2)
  max(1e5, floor(budget / 12))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# HDF5 helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.is_sparse_group <- function(obj) {
  inherits(obj, 'H5Group') && all(c('data', 'indices', 'indptr') %in% names(obj))
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

# Copy members of group `from` into group `to`, except `skip`
.copy_members <- function(from, to, skip = character(0)) {
  for (nm in setdiff(names(from), skip)) {
    to$obj_copy_from(src_loc = from, src_name = nm, dst_name = nm)
  }
  invisible(NULL)
}

# Placeholder compressed-sparse matrix: `n_major` columns (CSC) / rows (CSR),
# one entry with value `tag` at (minor index `at`, major index 0)
.write_placeholder <- function(parent, name, n_major, tag, at = 0L, template = NULL) {
  g <- parent$create_group(name)
  g$create_dataset('data', robj = tag, dtype = h5types$H5T_NATIVE_DOUBLE)
  g$create_dataset('indices', robj = as.integer(at), dtype = h5types$H5T_NATIVE_INT32)
  g$create_dataset('indptr', robj = c(0L, rep(1L, n_major)), dtype = h5types$H5T_NATIVE_INT32)
  if (!is.null(template)) .copy_attrs(template, g)
  invisible(g)
}

# All tagged placeholders in an open file: data.frame(path, tag)
.find_tags <- function(h5, prefix = '') {
  out <- data.frame(path = character(0), tag = numeric(0), stringsAsFactors = FALSE)
  for (nm in names(h5)) {
    path <- if (nzchar(prefix)) paste0(prefix, '/', nm) else nm
    obj <- tryCatch(h5[[nm]], error = function(e) NULL)
    if (!inherits(obj, 'H5Group')) next
    if (.is_sparse_group(obj) && identical(as.numeric(obj[['data']]$dims), 1)) {
      v <- as.numeric(obj[['data']][])
      if (v >= .STREAM_TAG_BASE && v < .STREAM_TAG_BASE + 1e5 && v == round(v)) {
        out <- rbind(out, data.frame(path = path, tag = v, stringsAsFactors = FALSE))
        next
      }
    }
    out <- rbind(out, .find_tags(obj, path))
  }
  out
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Moving the real matrices into place
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Replace the arrays of placeholder group `dst` by those of sparse group
# `src`. With `row_map` (new 0-based minor index for every source minor
# index, NA = drop) entries are filtered/reindexed in column blocks;
# otherwise the three datasets are copied by HDF5 directly.
.fill_sparse <- function(src, dst, row_map = NULL, verbose = FALSE) {
  for (nm in c('data', 'indices', 'indptr')) dst$link_delete(nm)
  if (is.null(row_map)) {
    .copy_members(src, dst, skip = setdiff(names(src), c('data', 'indices', 'indptr')))
    return(invisible(NULL))
  }
  indptr <- as.numeric(src[['indptr']]$read())
  n_major <- length(indptr) - 1
  block <- .stream_block()
  spans <- .major_blocks(indptr, block)
  # pass 1: how many entries survive per column
  kept <- numeric(n_major)
  for (b in seq_len(nrow(spans))) {
    rng <- .entry_range(indptr, spans[b, ])
    if (!length(rng)) next
    idx <- src[['indices']][rng]
    major <- rep.int(spans[b, 1]:spans[b, 2], diff(indptr[spans[b, 1]:(spans[b, 2] + 1)]))
    kept[spans[b, 1]:spans[b, 2]] <- tabulate(major[!is.na(row_map[idx + 1])] - spans[b, 1] + 1,
                                              nbins = spans[b, 2] - spans[b, 1] + 1)
  }
  new_ptr <- c(0, cumsum(kept))
  nnz <- new_ptr[n_major + 1]
  level <- GetCompressionLevel()
  chunk <- max(1, min(nnz, 65536))
  ptr_type <- if (nnz > .Machine$integer.max) h5types$H5T_STD_I64LE else h5types$H5T_NATIVE_INT32
  mk <- function(name, dtype, n) {
    dst$create_dataset(name = name, dtype = dtype, space = H5S$new(dims = n, maxdims = n),
                       chunk_dims = max(1, min(n, 65536)), gzip_level = level)
  }
  d_data <- mk('data', h5types$H5T_NATIVE_DOUBLE, nnz)
  d_idx <- mk('indices', h5types$H5T_NATIVE_INT32, nnz)
  # pass 2: filter, reindex, keep rows sorted within each column
  pos <- 0
  for (b in seq_len(nrow(spans))) {
    rng <- .entry_range(indptr, spans[b, ])
    if (!length(rng)) next
    idx <- src[['indices']][rng]
    val <- src[['data']][rng]
    major <- rep.int(spans[b, 1]:spans[b, 2], diff(indptr[spans[b, 1]:(spans[b, 2] + 1)]))
    new_idx <- row_map[idx + 1]
    keep <- !is.na(new_idx)
    ord <- order(major[keep], new_idx[keep])
    n_keep <- sum(keep)
    if (n_keep) {
      at <- (pos + 1):(pos + n_keep)
      d_data[at] <- as.numeric(val[keep][ord])
      d_idx[at] <- as.integer(new_idx[keep][ord])
      pos <- pos + n_keep
    }
    if (verbose && nrow(spans) > 1) message("  block ", b, "/", nrow(spans))
  }
  dst$create_dataset(name = 'indptr', robj = new_ptr, dtype = ptr_type,
                     chunk_dims = max(1, min(length(new_ptr), 65536)), gzip_level = level)
  invisible(NULL)
}

# Consecutive major-index ranges (1-based, inclusive) holding at most
# `block` entries each (a single major larger than that is its own block)
.major_blocks <- function(indptr, block) {
  n <- length(indptr) - 1
  if (n < 1) return(matrix(integer(0), ncol = 2))
  starts <- integer(0)
  ends <- integer(0)
  s <- 1
  while (s <= n) {
    e <- findInterval(indptr[s] + block, indptr) - 1
    e <- max(s, min(e, n))
    starts <- c(starts, s)
    ends <- c(ends, e)
    s <- e + 1
  }
  cbind(starts, ends)
}

# 1-based positions of the entries of majors span[1]..span[2]
.entry_range <- function(indptr, span) {
  a <- indptr[span[1]] + 1
  b <- indptr[span[2] + 1]
  if (b < a) integer(0) else seq.int(a, b)
}

# Sparse matrix that needs a real transpose (CSC h5ad X, or a writer that
# chose the other orientation): this one matrix is read into memory
.fill_transposed <- function(src, dst, verbose = FALSE) {
  if (verbose) message("  transposing ", src$get_obj_name(), " in memory (different orientation)")
  ptr <- src[['indptr']][]
  # rows from the recorded shape (h5ad `shape`, h5Seurat `shape`/`dims`);
  # the largest index would miss trailing empty rows
  key <- intersect(c('shape', 'dims'), h5attr_names(src))[1]
  n_rows <- as.integer(h5attr(src, key))[1]
  m <- Matrix::sparseMatrix(i = src[['indices']][] + 1, p = ptr, x = as.numeric(src[['data']][]),
                            dims = c(n_rows, length(ptr) - 1))
  t_m <- methods::as(Matrix::t(m), 'CsparseMatrix')
  for (nm in c('data', 'indices', 'indptr')) dst$link_delete(nm)
  level <- GetCompressionLevel()
  for (part in list(list('data', t_m@x), list('indices', t_m@i), list('indptr', t_m@p))) {
    dst$create_dataset(part[[1]], robj = part[[2]], gzip_level = level,
                       chunk_dims = max(1, min(length(part[[2]]), 65536)))
  }
  invisible(NULL)
}

# Replace a placeholder group (at `path` in `out`) with the dense dataset
# `src`; both formats lay dense matrices out the same way
.fill_dense <- function(src, out, path) {
  parent <- if (grepl('/', path)) out[[dirname(path)]] else out
  parent$link_delete(basename(path))
  parent$obj_copy_from(src_loc = src$get_file_id(), src_name = src$get_obj_name(),
                       dst_name = basename(path))
  invisible(NULL)
}

# Column sums and non-zero counts per major index (nCount / nFeature),
# streamed. Entries of one major index are contiguous, so per-major totals
# are differences of a running sum taken at the segment ends.
.major_stats <- function(src) {
  indptr <- as.numeric(src[['indptr']]$read())
  n <- length(indptr) - 1
  sums <- numeric(n)
  nfeat <- integer(n)
  spans <- .major_blocks(indptr, .stream_block())
  for (b in seq_len(nrow(spans))) {
    rng <- .entry_range(indptr, spans[b, ])
    if (!length(rng)) next
    val <- as.numeric(src[['data']][rng])
    majors <- spans[b, 1]:spans[b, 2]
    ends <- indptr[majors + 1] - indptr[spans[b, 1]]   # segment ends within the block
    cs <- c(0, cumsum(val))
    cp <- c(0, cumsum(val > 0))
    starts <- c(0, ends[-length(ends)])
    sums[majors] <- cs[ends + 1] - cs[starts + 1]
    nfeat[majors] <- as.integer(cp[ends + 1] - cp[starts + 1])
  }
  list(sum = sums, n = nfeat)
}

# Fill every tagged placeholder of `out` from the source matrix it stands
# for. `sources`: list indexed by tag (as character) of list(obj, row_map)
.fill_placeholders <- function(out, sources, verbose) {
  tags <- .find_tags(out)
  for (i in seq_len(nrow(tags))) {
    s <- sources[[as.character(tags$tag[i])]]
    if (is.null(s)) next
    if (verbose) message("Streaming ", s$label, " -> ", tags$path[i])
    dst <- out[[tags$path[i]]]
    if (inherits(s$obj, 'H5D')) {
      .fill_dense(s$obj, out, tags$path[i])
    } else if (isTRUE(s$transpose)) {
      .fill_transposed(s$obj, dst, verbose)
    } else {
      .fill_sparse(s$obj, dst, row_map = s$row_map, verbose = verbose)
    }
  }
  invisible(tags)
}

# Seurat recomputes nCount_*/nFeature_* from the (placeholder) counts when
# the source does not provide them; restore the real values
.patch_counts_meta <- function(out_md, counts_group, assay, verbose) {
  stats <- NULL
  for (col in c(paste0('nCount_', assay), paste0('nFeature_', assay))) {
    if (!out_md$exists(col)) next
    v <- out_md[[col]][]
    placeholder <- length(v) && all(v[-1] == 0) && v[1] <= 1 + .STREAM_TAG_BASE + 1e5
    if (!placeholder) next
    if (verbose) message("Recomputing ", col, " from the streamed counts")
    if (is.null(stats)) stats <- .major_stats(counts_group)
    new <- if (startsWith(col, 'nCount')) stats$sum else as.integer(stats$n)
    dt <- out_md[[col]]$get_type()
    out_md$link_delete(col)
    out_md$create_dataset(col, robj = new, dtype = dt)
  }
  invisible(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# h5ad -> h5Seurat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.h5ad_index_names <- function(df_group) {
  key <- if (df_group$attr_exists(attr_name = '_index')) h5attr(df_group, '_index') else '_index'
  if (!df_group$exists(key)) key <- intersect(c('_index', 'index'), names(df_group))[1]
  as.character(df_group[[key]][])
}

.stream_h5ad_to_h5seurat <- function(source, dest, assay = 'RNA', overwrite = FALSE,
                                     verbose = TRUE, ...) {
  if (file.exists(dest)) {
    if (!overwrite) stop("Destination file '", dest, "' already exists; set overwrite = TRUE", call. = FALSE)
    file.remove(dest)
  }
  src <- H5File$new(source, mode = 'r')
  on.exit(tryCatch(src$close_all(), error = function(e) NULL), add = TRUE)
  # raw/X lacking some of X's features makes the reader subset X itself;
  # that case is converted in memory
  if (identical(.h5ad_resolve_counts_source(src), 'raw/X') && src$exists('raw/var')) {
    missing_raw <- setdiff(gsub('_', '-', .h5ad_index_names(src[['var']]), fixed = TRUE),
                           gsub('_', '-', .h5ad_index_names(src[['raw/var']]), fixed = TRUE))
    if (length(missing_raw)) {
      if (verbose) message("raw/X lacks ", length(missing_raw), " of X's features; converting in memory")
      src$close_all()
      return(.hub_from_connection(source, dest, stype = 'h5ad', dtype = 'h5seurat', assay = assay,
                                  overwrite = overwrite, verbose = verbose))
    }
  }
  skel_path <- tempfile(fileext = '.h5ad')
  on.exit(unlink(skel_path), add = TRUE)
  skel <- H5File$new(skel_path, mode = 'w')
  .copy_attrs(src, skel)
  .copy_members(src, skel, skip = c('X', 'layers', 'raw'))

  n_obs <- length(.h5ad_index_names(src[['obs']]))
  x_genes <- .h5ad_index_names(src[['var']])
  sources <- list()
  k <- 0L
  add_matrix <- function(obj, parent, name, label, genes) {
    k <<- k + 1L
    tag <- .stream_tag(k)
    is_csc <- !inherits(obj, 'H5D') &&
      identical(tryCatch(h5attr(obj, 'encoding-type'), error = function(e) ''), 'csc_matrix')
    row_map <- NULL
    at <- 0L
    if (!identical(genes, x_genes)) {
      # raw/X with its own features: keep those shared with X, in X order
      # (matched after the reader's "_" -> "-" feature renaming)
      row_map <- match(gsub('_', '-', genes, fixed = TRUE), gsub('_', '-', x_genes, fixed = TRUE)) - 1L
      at <- max(0L, which(!is.na(row_map))[1] - 1L, na.rm = TRUE)
    }
    tmpl <- if (inherits(obj, 'H5Group')) obj else NULL
    shape <- c(n_obs, length(genes))
    g <- .write_placeholder(parent, name, n_major = n_obs, tag = tag, at = at, template = tmpl)
    for (a in c('encoding-type', 'encoding-version', 'shape')) {
      if (g$attr_exists(attr_name = a)) g$attr_delete(attr_name = a)
    }
    g$create_attr('encoding-type', 'csr_matrix', dtype = CachedGuessDType('csr_matrix'), space = ScalarSpace())
    g$create_attr('encoding-version', '0.1.0', dtype = CachedGuessDType('0.1.0'), space = ScalarSpace())
    g$create_attr('shape', as.integer(shape), dtype = GuessDType(1L))
    sources[[as.character(tag)]] <<- list(obj = obj, row_map = row_map, label = label,
                                          transpose = is_csc)
    invisible(tag)
  }
  if (src$exists('X')) add_matrix(src[['X']], skel, 'X', 'X', x_genes)
  if (src$exists('layers')) {
    lg <- skel$create_group('layers')
    .copy_attrs(src[['layers']], lg)
    for (nm in names(src[['layers']])) {
      add_matrix(src[['layers']][[nm]], lg, nm, paste0('layers/', nm), x_genes)
    }
  }
  if (src$exists('raw')) {
    rg <- skel$create_group('raw')
    .copy_attrs(src[['raw']], rg)
    .copy_members(src[['raw']], rg, skip = 'X')
    if (src[['raw']]$exists('X')) {
      raw_genes <- if (src[['raw']]$exists('var')) .h5ad_index_names(src[['raw']][['var']]) else x_genes
      add_matrix(src[['raw']][['X']], rg, 'X', 'raw/X', raw_genes)
    }
  }
  skel$close_all()

  if (verbose) message("Converting metadata, embeddings and graphs (matrices are streamed)")
  obj <- suppressWarnings(readH5AD(skel_path, assay.name = assay, verbose = FALSE))
  if (!is.null(obj@misc$scConvert_read)) obj@misc$scConvert_read$source_file <- normalizePath(source)
  writeH5Seurat(obj, filename = dest, overwrite = TRUE, verbose = FALSE)
  rm(obj)

  out <- H5File$new(dest, mode = 'r+')
  on.exit(tryCatch(out$close_all(), error = function(e) NULL), add = TRUE)
  tags <- .fill_placeholders(out, sources, verbose)
  if (src$exists('X') && !.stream_tag(1L) %in% tags$tag) {
    stop("Streaming conversion could not place X in the output; ",
         "retry with stream = FALSE", call. = FALSE)
  }
  # nCount/nFeature: from whatever was placed in the counts layer
  counts_path <- grep(paste0('^assays/', assay, '/(layers/)?counts$'), tags$path, value = TRUE)
  if (length(counts_path) && out$exists('meta.data')) {
    s <- sources[[as.character(tags$tag[tags$path == counts_path[1]])]]
    if (!inherits(s$obj, 'H5D')) {
      .patch_counts_meta(out[['meta.data']], out[[counts_path[1]]], assay, verbose)
    }
  }
  out$close_all()
  if (verbose) message("Streaming h5ad -> h5Seurat complete: ", dest)
  invisible(dest)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# h5Seurat -> h5ad
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.stream_h5seurat_to_h5ad <- function(source, dest, assay = NULL, overwrite = FALSE,
                                     verbose = TRUE, standardize = FALSE, ...) {
  if (file.exists(dest)) {
    if (!overwrite) stop("Destination file '", dest, "' already exists; set overwrite = TRUE", call. = FALSE)
    file.remove(dest)
  }
  src <- H5File$new(source, mode = 'r')
  on.exit(tryCatch(src$close_all(), error = function(e) NULL), add = TRUE)
  skel_path <- tempfile(fileext = '.h5seurat')
  on.exit(unlink(skel_path), add = TRUE)
  skel <- H5File$new(skel_path, mode = 'w')
  .copy_attrs(src, skel)
  .copy_members(src, skel, skip = 'assays')
  sources <- list()
  k <- 0L
  n_cells <- src[['cell.names']]$dims[1]
  placeholder_for <- function(obj, parent, name, label) {
    k <<- k + 1L
    tag <- .stream_tag(k)
    if (inherits(obj, 'H5D')) {
      # dense (e.g. scale.data): sparse stand-in with the same dims;
      # hdf5r reports h5Seurat dense layers as features x cells
      g <- .write_placeholder(parent, name, n_major = obj$dims[2], tag = tag)
      for (a in c('dims', 'shape')) g$create_attr(a, as.integer(obj$dims), dtype = GuessDType(1L))
    } else {
      .write_placeholder(parent, name, n_major = obj[['indptr']]$dims - 1, tag = tag, template = obj)
    }
    sources[[as.character(tag)]] <<- list(obj = obj, row_map = NULL, label = label)
  }
  ag_src <- src[['assays']]
  ag <- skel$create_group('assays')
  matrix_slots <- c('counts', 'data', 'scale.data')
  for (a in names(ag_src)) {
    s_a <- ag_src[[a]]
    d_a <- ag$create_group(a)
    .copy_attrs(s_a, d_a)
    .copy_members(s_a, d_a, skip = c('layers', matrix_slots))
    for (slot in intersect(matrix_slots, names(s_a))) {
      placeholder_for(s_a[[slot]], d_a, slot, paste0('assays/', a, '/', slot))
    }
    if (s_a$exists('layers') && inherits(s_a[['layers']], 'H5Group')) {
      s_l <- s_a[['layers']]
      d_l <- d_a$create_group('layers')
      .copy_attrs(s_l, d_l)
      for (lyr in names(s_l)) {
        obj <- s_l[[lyr]]
        if (inherits(obj, 'H5D') || .is_sparse_group(obj)) {
          placeholder_for(obj, d_l, lyr, paste0('assays/', a, '/layers/', lyr))
        } else {
          d_l$obj_copy_from(src_loc = s_l, src_name = lyr, dst_name = lyr)
        }
      }
    }
  }
  skel$close_all()

  if (verbose) message("Converting metadata, embeddings and graphs (matrices are streamed)")
  obj <- suppressWarnings(readH5Seurat(skel_path, verbose = FALSE))
  if (!is.null(assay) && assay %in% SeuratObject::Assays(obj)) SeuratObject::DefaultAssay(obj) <- assay
  .h5ad_saver(object = obj, filename = dest, overwrite = TRUE, verbose = FALSE,
              standardize = standardize)
  rm(obj)

  out <- H5File$new(dest, mode = 'r+')
  on.exit(tryCatch(out$close_all(), error = function(e) NULL), add = TRUE)
  # h5ad sparse matrices are CSR (cells x genes): the CSC arrays of the
  # h5Seurat genes x cells layers are copied as is; a CSC slot would need a
  # real transpose
  tags <- .find_tags(out)
  for (i in seq_len(nrow(tags))) {
    enc <- tryCatch(h5attr(out[[tags$path[i]]], 'encoding-type'), error = function(e) '')
    s <- sources[[as.character(tags$tag[i])]]
    if (!is.null(s) && identical(enc, 'csc_matrix') && !inherits(s$obj, 'H5D')) {
      sources[[as.character(tags$tag[i])]]$transpose <- TRUE
    }
  }
  .fill_placeholders(out, sources, verbose)
  out$close_all()
  if (verbose) message("Streaming h5Seurat -> h5ad complete: ", dest)
  invisible(dest)
}
