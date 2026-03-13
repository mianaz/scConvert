# Test scConvert.Seurat method and exported writeH5AD

library(scConvert)

skip_if_not_installed("Seurat")
skip_if_not_installed("Matrix")

# Helper function to create a test Seurat object
create_test_seurat <- function(n_cells = 50, n_features = 100) {
  set.seed(42)
  counts <- matrix(
    data = rpois(n_features * n_cells, lambda = 5),
    nrow = n_features,
    ncol = n_cells,
    dimnames = list(
      paste0("Gene", seq_len(n_features)),
      paste0("Cell", seq_len(n_cells))
    )
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "TestProject")
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj$cluster <- factor(sample(c("A", "B", "C"), n_cells, replace = TRUE))
  return(seurat_obj)
}

# -----------------------------------------------------------------------------
# scConvert.Seurat: h5Seurat output
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat, dest='*.h5seurat') creates valid h5Seurat file", {
  srt <- create_test_seurat()
  dest <- tempfile(fileext = ".h5Seurat")
  result <- scConvert(srt, dest = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))
  expect_true(hdf5r::is_hdf5(dest))
  expect_equal(result, dest)
})

# -----------------------------------------------------------------------------
# scConvert.Seurat: h5ad output
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat, dest='*.h5ad') creates valid h5ad file", {
  srt <- create_test_seurat()
  dest <- tempfile(fileext = ".h5ad")
  result <- scConvert(srt, dest = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))
  expect_true(hdf5r::is_hdf5(dest))
  expect_equal(result, dest)
})

# -----------------------------------------------------------------------------
# scConvert.Seurat: loom output
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat, dest='*.loom') creates valid loom file", {
  srt <- create_test_seurat()
  dest <- tempfile(fileext = ".loom")
  result <- scConvert(srt, dest = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))
  expect_true(hdf5r::is_hdf5(dest))
  expect_equal(result, dest)
})

# -----------------------------------------------------------------------------
# Bare type dest uses Project name
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat, dest='h5ad') uses Project name for filename", {
  srt <- create_test_seurat()
  # Run in a temp directory to avoid polluting working dir
  old_wd <- getwd()
  tmp_dir <- tempdir()
  setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  result <- scConvert(srt, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
  expect_equal(result, "TestProject.h5ad")
  expect_true(file.exists("TestProject.h5ad"))
  file.remove("TestProject.h5ad")
})

# -----------------------------------------------------------------------------
# Error: missing dest
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat) errors when dest is missing", {
  srt <- create_test_seurat()
  expect_error(scConvert(srt), "dest.*must be provided")
})

# -----------------------------------------------------------------------------
# Error: overwrite = FALSE on existing file
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat, dest='*.h5ad') errors when file exists and overwrite=FALSE", {
  srt <- create_test_seurat()
  dest <- tempfile(fileext = ".h5ad")
  # Create the file first
  scConvert(srt, dest = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))
  expect_error(
    scConvert(srt, dest = dest, overwrite = FALSE, verbose = FALSE),
    "exists"
  )
})

# -----------------------------------------------------------------------------
# Error: unsupported format
# -----------------------------------------------------------------------------

test_that("scConvert(Seurat) errors on unsupported format", {
  srt <- create_test_seurat()
  expect_error(
    scConvert(srt, dest = "output.csv", verbose = FALSE),
    "Cannot convert Seurat objects"
  )
})

# -----------------------------------------------------------------------------
# writeH5AD exported convenience function
# -----------------------------------------------------------------------------

test_that("writeH5AD() creates valid h5ad file", {
  srt <- create_test_seurat()
  dest <- tempfile(fileext = ".h5ad")
  result <- writeH5AD(srt, filename = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(dest))
  expect_true(hdf5r::is_hdf5(dest))
  expect_equal(result, dest)
})

test_that("writeH5AD() adds .h5ad extension if missing", {
  srt <- create_test_seurat()
  dest <- tempfile()
  result <- writeH5AD(srt, filename = dest, overwrite = TRUE, verbose = FALSE)
  expect_true(grepl("\\.h5ad$", result))
  expect_true(file.exists(result))
})

test_that("writeH5AD() errors on non-Seurat input", {
  expect_error(writeH5AD("not_a_seurat"), "must be a Seurat object")
})
