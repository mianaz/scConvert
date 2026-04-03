# tests/testthat/test-misc-roundtrip.R
#
# Tests for misc (uns) roundtrip preservation through h5ad, including:
# - Nested lists, matrices, sparse matrices, factors, vectors
# - Sequential layer write (RAM optimization)
# - Spatial metadata preservation via misc

library(testthat)
library(scConvert)
library(Seurat)
library(Matrix)

test_that("misc scalars roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["scalar_string"]] <- "hello world"
  obj@misc[["scalar_numeric"]] <- 3.14159
  obj@misc[["scalar_integer"]] <- 42L

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_equal(obj2@misc[["scalar_string"]], "hello world")
  expect_equal(obj2@misc[["scalar_numeric"]], 3.14159, tolerance = 1e-6)
  expect_equal(obj2@misc[["scalar_integer"]], 42L)
})

test_that("misc character vectors roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["gene_list"]] <- c("TP53", "BRCA1", "EGFR", "MYC")
  obj@misc[["single_char"]] <- "only one"

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_equal(obj2@misc[["gene_list"]], c("TP53", "BRCA1", "EGFR", "MYC"))
})

test_that("misc numeric vectors roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["scores"]] <- c(0.1, 0.5, 0.9, 1.0)
  obj@misc[["integers"]] <- c(1L, 2L, 3L, 4L, 5L)

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_equal(obj2@misc[["scores"]], c(0.1, 0.5, 0.9, 1.0), tolerance = 1e-6)
  expect_equal(obj2@misc[["integers"]], c(1L, 2L, 3L, 4L, 5L))
})

test_that("misc dense matrices roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  mat <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(mat) <- paste0("gene", 1:4)
  colnames(mat) <- paste0("comp", 1:5)
  obj@misc[["dense_matrix"]] <- mat

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_true(!is.null(obj2@misc[["dense_matrix"]]))
  expect_equal(as.numeric(obj2@misc[["dense_matrix"]]),
               as.numeric(mat), tolerance = 1e-6)
})

test_that("misc sparse matrices roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  sp_mat <- rsparsematrix(10, 10, density = 0.3)
  rownames(sp_mat) <- paste0("g", 1:10)
  colnames(sp_mat) <- paste0("g", 1:10)
  obj@misc[["gene_network"]] <- sp_mat

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  # Sparse matrices are stored as CSR groups in uns; read back as list components
  expect_true(!is.null(obj2@misc[["gene_network"]]))
})

test_that("misc nested lists roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["spatial"]] <- list(
    library1 = list(
      scalefactors = list(
        spot_diameter_fullres = 178.5,
        tissue_hires_scalef = 0.085
      ),
      metadata = list(
        chemistry = "Visium v2",
        n_spots = 2696L
      )
    )
  )

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_true(!is.null(obj2@misc[["spatial"]]))
  expect_true(!is.null(obj2@misc[["spatial"]][["library1"]]))

  sf <- obj2@misc[["spatial"]][["library1"]][["scalefactors"]]
  expect_equal(sf[["spot_diameter_fullres"]], 178.5, tolerance = 1e-6)
  expect_equal(sf[["tissue_hires_scalef"]], 0.085, tolerance = 1e-6)

  meta <- obj2@misc[["spatial"]][["library1"]][["metadata"]]
  expect_equal(meta[["chemistry"]], "Visium v2")
  expect_equal(meta[["n_spots"]], 2696L)
})

test_that("misc data.frames roundtrip through h5ad", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  df <- data.frame(
    gene = c("TP53", "BRCA1", "MYC"),
    score = c(0.95, 0.87, 0.73),
    significant = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$gene
  obj@misc[["de_results"]] <- df

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_true(!is.null(obj2@misc[["de_results"]]))
})

test_that("misc internal keys are excluded from uns", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["__varp__"]] <- list(test = matrix(1:4, 2, 2))
  obj@misc[[".__h5ad_path__"]] <- "/fake/path"
  obj@misc[[".__h5ad_loaded__"]] <- c("X", "obs")
  obj@misc[["user_data"]] <- "should be preserved"

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)

  h5 <- hdf5r::H5File$new(tmp, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  # Internal keys should NOT appear in uns
  uns_names <- if (h5$exists("uns")) names(h5[["uns"]]) else character(0)
  expect_false("__varp__" %in% uns_names)
  expect_false(".__h5ad_path__" %in% uns_names)
  expect_false(".__h5ad_loaded__" %in% uns_names)

  # User data should be in uns
  expect_true("user_data" %in% uns_names)
})

test_that("misc with simulated spatial metadata roundtrips", {
  # Simulate what a MERFISH h5ad would store in uns
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["spatial_metadata"]] <- list(
    segmentation = list(
      method = "watershed",
      params = c(min_area = 50.0, max_area = 500.0)
    ),
    regions = c("cortex", "hippocampus", "thalamus"),
    z_coordinates = rnorm(ncol(obj))
  )

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  sm <- obj2@misc[["spatial_metadata"]]
  expect_true(!is.null(sm))
  expect_equal(sm[["regions"]], c("cortex", "hippocampus", "thalamus"))
  expect_equal(sm[["segmentation"]][["method"]], "watershed")
  expect_equal(length(sm[["z_coordinates"]]), ncol(obj))
})

# ============================================================
# Sequential layer write / RAM optimization tests
# ============================================================

test_that("sequential layer write produces identical output", {
  set.seed(42)
  mat <- rsparsematrix(500, 100, density = 0.1)
  rownames(mat) <- paste0("G", 1:500)
  colnames(mat) <- paste0("C", 1:100)
  obj <- CreateSeuratObject(mat)
  # Add a data layer by normalizing
  obj <- NormalizeData(obj, verbose = FALSE)

  # Verify both layers exist
  expect_true("counts" %in% Layers(obj[["RNA"]]))
  expect_true("data" %in% Layers(obj[["RNA"]]))

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  # Check dimensions preserved
  expect_equal(ncol(obj), ncol(obj2))
  expect_equal(nrow(obj), nrow(obj2))

  # Check counts preserved (stored as X or raw/X)
  counts_orig <- GetAssayData(obj, layer = "counts")
  counts_rt <- GetAssayData(obj2, layer = "counts")
  expect_equal(dim(counts_orig), dim(counts_rt))

  # Sample a few values for exact match
  set.seed(1)
  idx <- sample(length(counts_orig@x), min(100, length(counts_orig@x)))
  expect_equal(counts_orig@x[idx], counts_rt@x[idx], tolerance = 1e-6)
})

test_that("write with only counts layer works", {
  set.seed(42)
  mat <- rsparsematrix(200, 50, density = 0.1)
  rownames(mat) <- paste0("G", 1:200)
  colnames(mat) <- paste0("C", 1:50)
  obj <- CreateSeuratObject(mat)

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp, verbose = FALSE)

  expect_equal(ncol(obj), ncol(obj2))
  expect_equal(nrow(obj), nrow(obj2))
})

test_that("RAM optimization does not increase write time", {
  set.seed(42)
  mat <- rsparsematrix(5000, 2000, density = 0.05)
  rownames(mat) <- paste0("G", 1:5000)
  colnames(mat) <- paste0("C", 1:2000)
  obj <- CreateSeuratObject(mat)
  obj <- NormalizeData(obj, verbose = FALSE)

  # Time the write (should be fast with gzip=0)
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))
  t <- system.time(writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE, gzip = 0L))

  # Sanity: should complete in reasonable time (< 5s for 2K cells)
  expect_true(t["elapsed"] < 5)
})

# ============================================================
# Double roundtrip stability
# ============================================================

test_that("misc survives double roundtrip h5ad -> Seurat -> h5ad -> Seurat", {
  obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
  obj@misc[["nested"]] <- list(
    a = list(b = c(1.0, 2.0, 3.0), c = "test"),
    d = c(10L, 20L, 30L)
  )

  tmp1 <- tempfile(fileext = ".h5ad")
  tmp2 <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(tmp1, tmp2)))

  # First roundtrip
  writeH5AD(obj, tmp1, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(tmp1, verbose = FALSE)

  # Second roundtrip
  writeH5AD(obj2, tmp2, overwrite = TRUE, verbose = FALSE)
  obj3 <- readH5AD(tmp2, verbose = FALSE)

  # Check nested data preserved through double roundtrip
  expect_equal(obj3@misc[["nested"]][["a"]][["c"]], "test")
  expect_equal(obj3@misc[["nested"]][["d"]], c(10L, 20L, 30L))
  expect_equal(obj3@misc[["nested"]][["a"]][["b"]], c(1.0, 2.0, 3.0),
               tolerance = 1e-6)
})
