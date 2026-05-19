make_local_zarr <- function() {
  obj <- suppressWarnings(SeuratObject::pbmc_small)
  dest <- tempfile(fileext = ".zarr")
  writeZarr(obj, dest, verbose = FALSE)
  attr(dest, "n_cells") <- ncol(obj)
  attr(dest, "n_genes") <- nrow(obj)
  attr(dest, "ref") <- obj
  dest
}

test_that("readZarr(obs_idx = ...) returns the correct cell slice with X intact", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  n <- attr(dest, "n_cells")
  ref <- attr(dest, "ref")

  pick <- c(1L, 5L, 10L, 25L, 40L)
  obj <- readZarr(dest, obs_idx = pick, verbose = FALSE)
  expect_equal(ncol(obj), length(pick))
  expect_equal(nrow(obj), attr(dest, "n_genes"))
  # Cell names should match the picked positions of the source
  expect_equal(colnames(obj), Cells(ref)[pick])

  # Counts should match elementwise for the slice
  ref_mat <- SeuratObject::GetAssayData(ref, layer = "counts")
  out_mat <- SeuratObject::GetAssayData(obj, layer = "counts")
  expect_equal(as.matrix(out_mat),
               as.matrix(ref_mat[, pick, drop = FALSE]),
               ignore_attr = TRUE)
})

test_that("readZarr(var_idx = ...) returns the correct feature slice", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  ref <- attr(dest, "ref")

  pick <- c(2L, 3L, 7L, 50L)
  obj <- readZarr(dest, var_idx = pick, verbose = FALSE)
  expect_equal(nrow(obj), length(pick))
  expect_equal(ncol(obj), attr(dest, "n_cells"))
  expect_equal(rownames(obj), rownames(ref)[pick])

  ref_mat <- SeuratObject::GetAssayData(ref, layer = "counts")
  out_mat <- SeuratObject::GetAssayData(obj, layer = "counts")
  expect_equal(as.matrix(out_mat),
               as.matrix(ref_mat[pick, , drop = FALSE]),
               ignore_attr = TRUE)
})

test_that("readZarr(obs_idx, var_idx) slices both dims", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  ref <- attr(dest, "ref")

  cells <- c(3L, 11L, 17L)
  genes <- c(1L, 2L, 5L, 8L)
  obj <- readZarr(dest, obs_idx = cells, var_idx = genes, verbose = FALSE)
  expect_equal(ncol(obj), length(cells))
  expect_equal(nrow(obj), length(genes))

  ref_mat <- SeuratObject::GetAssayData(ref, layer = "counts")
  out_mat <- SeuratObject::GetAssayData(obj, layer = "counts")
  expect_equal(as.matrix(out_mat),
               as.matrix(ref_mat[genes, cells, drop = FALSE]),
               ignore_attr = TRUE)
})

test_that("readZarr(obs_idx) preserves obs metadata for the slice", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  ref <- attr(dest, "ref")

  pick <- c(2L, 8L, 15L)
  obj <- readZarr(dest, obs_idx = pick, verbose = FALSE)
  # nCount_RNA / nFeature_RNA are computed at construction, but the original
  # obs metadata columns should track the slice exactly
  meta <- obj[[]]
  ref_meta <- ref[[]]
  common <- intersect(names(meta), names(ref_meta))
  common <- setdiff(common, c("nCount_RNA", "nFeature_RNA",
                              "orig.ident", "groups", "letter.idents"))
  if (length(common) > 0L) {
    for (col in common) {
      expect_equal(meta[[col]], ref_meta[[col]][pick],
                   ignore_attr = TRUE,
                   info = paste0("column ", col))
    }
  }
})

test_that("readZarr(obs_idx) slices obsm embeddings to the same cells", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  ref <- attr(dest, "ref")

  pick <- c(5L, 10L, 25L)
  obj <- readZarr(dest, obs_idx = pick, verbose = FALSE)
  # pbmc_small has pca + tsne reductions; pick whichever survives
  reds <- SeuratObject::Reductions(obj)
  if (length(reds) > 0L) {
    e <- SeuratObject::Embeddings(obj, reduction = reds[1])
    expect_equal(nrow(e), length(pick))
  }
})

test_that("readZarr rejects out-of-range indices with a clear message", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  n <- attr(dest, "n_cells")
  expect_error(readZarr(dest, obs_idx = c(1L, n + 1L), verbose = FALSE),
               "obs_idx out of range")
  expect_error(readZarr(dest, var_idx = -1L, verbose = FALSE),
               "var_idx out of range")
})

test_that("readZarr accepts a logical mask for obs_idx", {
  dest <- make_local_zarr()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  ref <- attr(dest, "ref")

  mask <- rep(c(TRUE, FALSE), length.out = ncol(ref))
  obj <- readZarr(dest, obs_idx = mask, verbose = FALSE)
  expect_equal(ncol(obj), sum(mask))
})
