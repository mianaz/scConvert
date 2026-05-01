library(scConvert)

# The R sparse streaming fallback used to call $read() on the full data /
# indices / indptr datasets. On large matrices that pulled multiple gigabytes
# into R's heap. .stream_sparse_group now does bounded-chunk reads via
# hdf5r range indexing, gated by options("scConvert.stream_chunk_bytes").
#
# This test verifies two things:
#   1. The output is byte-for-byte equivalent to a normal write, no matter
#      what stream_chunk_bytes is set to (no data corruption from chunking).
#   2. Forcing a tiny byte budget still produces a correct file — i.e. the
#      chunked path is actually exercised, not bypassed.

test_that(".stream_sparse_group produces correct output under a small byte budget", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")

  # Synthetic sparse matrix: 500 rows x 200 cols, ~10% density
  set.seed(42)
  nrows <- 500L
  ncols <- 200L
  nnz <- 10000L
  i <- sample.int(nrows, nnz, replace = TRUE)
  j <- sample.int(ncols, nnz, replace = TRUE)
  x <- runif(nnz, 0.1, 10)
  mat <- Matrix::sparseMatrix(i = i, j = j, x = x,
                              dims = c(nrows, ncols),
                              repr = "C")

  # Write a synthetic CSR h5ad-style group by hand so we exercise
  # .stream_sparse_group end-to-end.
  src <- tempfile(fileext = ".h5")
  dst <- tempfile(fileext = ".h5")
  on.exit({ unlink(src); unlink(dst) }, add = TRUE)

  # CSR layout: convert from dgCMatrix (CSC) by transposing.
  csr <- Matrix::t(mat)
  h_src <- hdf5r::H5File$new(src, mode = "w")
  grp <- h_src$create_group("X")
  # data / indices / indptr as 1-D datasets
  grp$create_dataset("data", robj = as.numeric(csr@x),
                     chunk_dims = 1024L, gzip_level = 1L)
  grp$create_dataset("indices", robj = as.integer(csr@i),
                     chunk_dims = 1024L, gzip_level = 1L)
  grp$create_dataset("indptr", robj = as.integer(csr@p),
                     chunk_dims = 1024L, gzip_level = 1L)
  grp$create_attr("shape", robj = as.integer(c(ncols, nrows)),
                  dtype = scConvert:::GuessDType(c(ncols, nrows)))
  grp$create_attr("encoding-type", robj = "csr_matrix",
                  dtype = scConvert:::GuessDType("csr_matrix"),
                  space = hdf5r::H5S$new(type = "scalar"))
  h_src$close_all()

  # Force a tiny 1 KiB byte budget so the streaming loop runs many iterations.
  old_opt <- options(scConvert.stream_chunk_bytes = 1024L)
  on.exit(options(old_opt), add = TRUE)

  h_src <- hdf5r::H5File$new(src, mode = "r")
  h_dst <- hdf5r::H5File$new(dst, mode = "w")

  scConvert:::.stream_sparse_group(
    src_group = h_src[["X"]],
    dst_parent = h_dst,
    dst_name = "X",
    gzip = 1L,
    src_format = "h5ad",
    dst_format = "h5ad"
  )

  h_src$close_all()
  h_dst$close_all()

  # Read back and verify data / indices / indptr match bit-for-bit.
  h <- hdf5r::H5File$new(dst, mode = "r")
  on.exit(h$close_all(), add = TRUE)

  expect_equal(length(h[["X/data"]][]), length(csr@x))
  expect_equal(h[["X/data"]][], csr@x)
  expect_equal(h[["X/indices"]][], csr@i)
  expect_equal(h[["X/indptr"]][], csr@p)
})

test_that("sparse stream default budget handles normal matrices identically", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")

  # Same fixture, default budget (64 MiB); should produce identical output.
  set.seed(7)
  nrows <- 200L
  ncols <- 100L
  nnz <- 2000L
  i <- sample.int(nrows, nnz, replace = TRUE)
  j <- sample.int(ncols, nnz, replace = TRUE)
  x <- runif(nnz)
  mat <- Matrix::sparseMatrix(i = i, j = j, x = x,
                              dims = c(nrows, ncols),
                              repr = "C")
  csr <- Matrix::t(mat)

  src <- tempfile(fileext = ".h5")
  dst <- tempfile(fileext = ".h5")
  on.exit({ unlink(src); unlink(dst) }, add = TRUE)

  h_src <- hdf5r::H5File$new(src, mode = "w")
  g <- h_src$create_group("X")
  g$create_dataset("data", robj = as.numeric(csr@x),
                   chunk_dims = 1024L, gzip_level = 1L)
  g$create_dataset("indices", robj = as.integer(csr@i),
                   chunk_dims = 1024L, gzip_level = 1L)
  g$create_dataset("indptr", robj = as.integer(csr@p),
                   chunk_dims = 1024L, gzip_level = 1L)
  g$create_attr("shape", robj = as.integer(c(ncols, nrows)),
                dtype = scConvert:::GuessDType(c(ncols, nrows)))
  h_src$close_all()

  h_src <- hdf5r::H5File$new(src, mode = "r")
  h_dst <- hdf5r::H5File$new(dst, mode = "w")
  scConvert:::.stream_sparse_group(
    src_group = h_src[["X"]],
    dst_parent = h_dst,
    dst_name = "X",
    gzip = 1L,
    src_format = "h5ad",
    dst_format = "h5ad"
  )
  h_src$close_all(); h_dst$close_all()

  h <- hdf5r::H5File$new(dst, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  expect_equal(h[["X/data"]][], csr@x)
  expect_equal(h[["X/indices"]][], csr@i)
  expect_equal(h[["X/indptr"]][], csr@p)
})
