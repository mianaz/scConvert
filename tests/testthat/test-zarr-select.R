test_that(".zarr_select_keep returns the right slice for each filter mode", {
  available <- c("X_pca", "X_umap", "X_harmony")
  expect_equal(scConvert:::.zarr_select_keep(NULL, available), available)
  expect_equal(scConvert:::.zarr_select_keep(character(0), available),
               character(0))
  expect_equal(scConvert:::.zarr_select_keep(c("X_pca"), available), "X_pca")
  expect_equal(scConvert:::.zarr_select_keep(c("X_pca", "missing"), available),
               "X_pca")
})

test_that(".zarr_select_keep errors on non-character non-NULL", {
  expect_error(scConvert:::.zarr_select_keep(1L, c("a")),
               "filter must be NULL or a character")
})

# Set up a small local zarr store for end-to-end selection tests.
local_zarr_fixture <- function() {
  obj <- suppressWarnings(SeuratObject::pbmc_small)
  dest <- tempfile(fileext = ".zarr")
  writeZarr(obj, dest, verbose = FALSE)
  dest
}

test_that("readZarr(layers = character(0)) skips all layer reads", {
  dest <- local_zarr_fixture()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  obj <- readZarr(dest, layers = character(0), verbose = FALSE)
  expect_s4_class(obj, "Seurat")
  # Layers other than counts/data should not be present beyond defaults
  expect_true(length(SeuratObject::Layers(obj)) <= 2L)
})

test_that("readZarr(obsm = character(0)) drops all reductions", {
  dest <- local_zarr_fixture()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  obj <- readZarr(dest, obsm = character(0), verbose = FALSE)
  expect_length(SeuratObject::Reductions(obj), 0L)
})

test_that("readZarr(obsm = 'X_pca') keeps only the named reduction", {
  obj <- suppressWarnings(SeuratObject::pbmc_small)
  dest <- tempfile(fileext = ".zarr")
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  writeZarr(obj, dest, verbose = FALSE)

  loaded <- readZarr(dest, obsm = "X_pca", verbose = FALSE)
  # pbmc_small writes pca + tsne; with the filter only pca should land
  expect_true("pca" %in% SeuratObject::Reductions(loaded))
  expect_false("tsne" %in% SeuratObject::Reductions(loaded))
})

test_that("readZarr(include_x = FALSE) returns a Seurat with empty X", {
  dest <- local_zarr_fixture()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  obj <- readZarr(dest, include_x = FALSE, verbose = FALSE)
  expect_s4_class(obj, "Seurat")
  # All zero entries in the counts layer
  expect_equal(sum(SeuratObject::GetAssayData(obj, layer = "counts")), 0)
  # Cell + gene names still populated
  expect_gt(ncol(obj), 0L)
  expect_gt(nrow(obj), 0L)
})

test_that("readZarr(uns = character(0)) skips uns reads", {
  dest <- local_zarr_fixture()
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))
  obj <- readZarr(dest, uns = character(0), verbose = FALSE)
  # misc should not contain typical uns items
  expect_false(any(grepl("_present$", names(obj@misc))))
})
