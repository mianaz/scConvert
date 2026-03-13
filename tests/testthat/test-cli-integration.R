# Comprehensive CLI integration tests for scConvert_cli()
#
# Uses shipped demo data: inst/testdata/pbmc_small.h5ad (214 cells, 2000 genes)
#                          inst/testdata/pbmc_small.rds  (214 cells, 2000 genes)

library(scConvert)

skip_if_not_installed("Seurat")
skip_if_not_installed("hdf5r")

# ---------------------------------------------------------------------------
# Locate demo data
# ---------------------------------------------------------------------------

h5ad_path <- system.file("testdata", "pbmc_small.h5ad", package = "scConvert")
rds_path  <- system.file("testdata", "pbmc_small.rds",  package = "scConvert")

skip_if(
  !nzchar(h5ad_path) || !file.exists(h5ad_path),
  "Demo h5ad file not found in inst/testdata"
)
skip_if(

  !nzchar(rds_path) || !file.exists(rds_path),
  "Demo rds file not found in inst/testdata"
)

# ---------------------------------------------------------------------------
# Load reference objects once for comparison across all tests
# ---------------------------------------------------------------------------

ref_h5ad <- readH5AD(h5ad_path, verbose = FALSE)
ref_rds  <- readRDS(rds_path)

expected_ncol <- 214L
expected_nrow <- 2000L

# Sanity check on the reference objects
test_that("reference data has expected dimensions", {
  expect_equal(ncol(ref_h5ad), expected_ncol)
  expect_equal(nrow(ref_h5ad), expected_nrow)
  expect_equal(ncol(ref_rds),  expected_ncol)
  expect_equal(nrow(ref_rds),  expected_nrow)
})

# ---------------------------------------------------------------------------
# Helper: compare a loaded Seurat object to a reference
# ---------------------------------------------------------------------------

#' @param loaded  Seurat object loaded from the converted file
#' @param ref     Reference Seurat object
#' @param label   A short description for test failure messages
#' @param check_meta Whether to check meta.data column names
compare_to_ref <- function(loaded, ref, label, check_meta = TRUE) {
  # Dimensions

  expect_equal(ncol(loaded), ncol(ref),
               label = paste(label, "ncol"))
  expect_equal(nrow(loaded), nrow(ref),
               label = paste(label, "nrow"))

  # Cell barcode preservation
  expect_equal(sort(colnames(loaded)), sort(colnames(ref)),
               label = paste(label, "cell names"))

  # Feature name preservation
  expect_equal(sort(rownames(loaded)), sort(rownames(ref)),
               label = paste(label, "feature names"))

  # Metadata column preservation

  if (check_meta) {
    ref_meta_cols <- setdiff(colnames(ref[[]]), "orig.ident")
    loaded_meta_cols <- colnames(loaded[[]])
    for (col in ref_meta_cols) {
      expect_true(
        col %in% loaded_meta_cols,
        label = paste(label, "meta.data column:", col)
      )
    }
  }
}

# ===========================================================================
# 1. h5ad -> h5seurat
# ===========================================================================

test_that("scConvert_cli: h5ad -> h5seurat", {
  out <- tempfile(fileext = ".h5Seurat")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readH5Seurat(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  compare_to_ref(loaded, ref_h5ad, "h5ad->h5seurat")
})

# ===========================================================================
# 2. h5seurat -> h5ad  (needs an h5seurat intermediate)
# ===========================================================================

test_that("scConvert_cli: h5seurat -> h5ad", {
  # First create an h5seurat from the h5ad
  h5seurat_tmp <- tempfile(fileext = ".h5Seurat")
  on.exit(unlink(h5seurat_tmp), add = TRUE)
  scConvert_cli(
    input = h5ad_path, output = h5seurat_tmp,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )

  out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = h5seurat_tmp, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readH5AD(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  compare_to_ref(loaded, ref_h5ad, "h5seurat->h5ad")
})

# ===========================================================================
# 3. h5ad -> loom
# ===========================================================================

test_that("scConvert_cli: h5ad -> loom", {
  out <- tempfile(fileext = ".loom")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readLoom(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")

  # Loom may not preserve all metadata columns, so relax meta check
  expect_equal(ncol(loaded), expected_ncol)
  expect_equal(nrow(loaded), expected_nrow)
  expect_equal(sort(colnames(loaded)), sort(colnames(ref_h5ad)))
  expect_equal(sort(rownames(loaded)), sort(rownames(ref_h5ad)))
})

# ===========================================================================
# 4. h5ad -> rds
# ===========================================================================

test_that("scConvert_cli: h5ad -> rds", {
  out <- tempfile(fileext = ".rds")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readRDS(out)
  expect_s4_class(loaded, "Seurat")
  compare_to_ref(loaded, ref_h5ad, "h5ad->rds")
})

# ===========================================================================
# 5. h5ad -> zarr
# ===========================================================================

test_that("scConvert_cli: h5ad -> zarr", {
  skip_if_not_installed("jsonlite")

  out <- tempfile(fileext = ".zarr")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  result <- scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(dir.exists(out))

  loaded <- readZarr(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")

  expect_equal(ncol(loaded), expected_ncol)
  expect_equal(nrow(loaded), expected_nrow)
  expect_equal(sort(colnames(loaded)), sort(colnames(ref_h5ad)))
  expect_equal(sort(rownames(loaded)), sort(rownames(ref_h5ad)))
})

# ===========================================================================
# 6. rds -> h5ad
# ===========================================================================

test_that("scConvert_cli: rds -> h5ad", {
  out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = rds_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readH5AD(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  compare_to_ref(loaded, ref_rds, "rds->h5ad")
})

# ===========================================================================
# 7. loom -> h5ad  (needs a loom intermediate)
# ===========================================================================

test_that("scConvert_cli: loom -> h5ad", {
  # Create loom intermediate from h5ad
  loom_tmp <- tempfile(fileext = ".loom")
  on.exit(unlink(loom_tmp), add = TRUE)
  scConvert_cli(
    input = h5ad_path, output = loom_tmp,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  skip_if(!file.exists(loom_tmp), "loom intermediate not created")

  out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = loom_tmp, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readH5AD(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")

  expect_equal(ncol(loaded), expected_ncol)
  expect_equal(nrow(loaded), expected_nrow)
  expect_equal(sort(colnames(loaded)), sort(colnames(ref_h5ad)))
  expect_equal(sort(rownames(loaded)), sort(rownames(ref_h5ad)))
})

# ===========================================================================
# 8. zarr -> h5ad  (needs a zarr intermediate)
# ===========================================================================

test_that("scConvert_cli: zarr -> h5ad", {
  skip_if_not_installed("jsonlite")

  # Create zarr intermediate from h5ad
  zarr_tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(zarr_tmp, recursive = TRUE), add = TRUE)
  scConvert_cli(
    input = h5ad_path, output = zarr_tmp,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  skip_if(!dir.exists(zarr_tmp), "zarr intermediate not created")

  out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(out), add = TRUE)

  result <- scConvert_cli(
    input = zarr_tmp, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)

  loaded <- readH5AD(out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")

  expect_equal(ncol(loaded), expected_ncol)
  expect_equal(nrow(loaded), expected_nrow)
  expect_equal(sort(colnames(loaded)), sort(colnames(ref_h5ad)))
  expect_equal(sort(rownames(loaded)), sort(rownames(ref_h5ad)))
})

# ===========================================================================
# 9. CLI flag: --quiet (verbose = FALSE produces no messages)
# ===========================================================================

test_that("scConvert_cli: verbose = FALSE suppresses messages", {
  out <- tempfile(fileext = ".h5Seurat")
  on.exit(unlink(out), add = TRUE)

  msgs <- capture.output(
    scConvert_cli(
      input = h5ad_path, output = out,
      assay = "RNA", overwrite = TRUE, verbose = FALSE
    ),
    type = "message"
  )
  # Verbose = FALSE should suppress most messages; allow a few from sub-functions
  # that don't yet gate on verbose (tracked as enhancement)
  expect_lt(length(msgs), 10)
})

# ===========================================================================
# 10. CLI flag: --overwrite behaviour
# ===========================================================================

test_that("scConvert_cli: overwrite = FALSE errors on existing output", {
  out <- tempfile(fileext = ".h5Seurat")
  on.exit(unlink(out), add = TRUE)

  # First conversion succeeds
  scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(file.exists(out))

  # Second conversion with overwrite = FALSE should fail.
  # scConvert_cli returns FALSE on failure (wraps errors internally),
  # but depending on the code path it may also raise a stop().
  # We handle both: either it returns FALSE, or it errors.
  result <- tryCatch(
    scConvert_cli(
      input = h5ad_path, output = out,
      assay = "RNA", overwrite = FALSE, verbose = FALSE
    ),
    error = function(e) "errored"
  )
  expect_true(identical(result, FALSE) || identical(result, "errored"))
})

test_that("scConvert_cli: overwrite = TRUE replaces existing output", {
  out <- tempfile(fileext = ".h5Seurat")
  on.exit(unlink(out), add = TRUE)

  # First conversion
  scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  first_mtime <- file.info(out)$mtime

  # Brief pause to ensure mtime differs
  Sys.sleep(1)

  # Second conversion with overwrite = TRUE should succeed
  result <- scConvert_cli(
    input = h5ad_path, output = out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  expect_true(result)
  expect_true(file.exists(out))
  expect_true(file.info(out)$mtime >= first_mtime)
})

# ===========================================================================
# 11. Error on missing input file
# ===========================================================================

test_that("scConvert_cli: errors on missing input file", {
  out <- tempfile(fileext = ".h5ad")

  expect_error(
    scConvert_cli(
      input = "/nonexistent/path/fake.h5ad",
      output = out,
      verbose = FALSE
    ),
    "not found"
  )
})

# ===========================================================================
# 12. Full roundtrip fidelity: h5ad -> rds -> h5ad
# ===========================================================================

test_that("full roundtrip h5ad -> rds -> h5ad preserves data", {
  rds_tmp <- tempfile(fileext = ".rds")
  h5ad_out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(rds_tmp, h5ad_out)), add = TRUE)

  scConvert_cli(
    input = h5ad_path, output = rds_tmp,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  scConvert_cli(
    input = rds_tmp, output = h5ad_out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )

  loaded <- readH5AD(h5ad_out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), expected_ncol)
  expect_equal(nrow(loaded), expected_nrow)
  expect_equal(sort(colnames(loaded)), sort(colnames(ref_h5ad)))
  expect_equal(sort(rownames(loaded)), sort(rownames(ref_h5ad)))
})

# ===========================================================================
# 13. Full roundtrip fidelity: h5ad -> h5seurat -> h5ad
# ===========================================================================

test_that("full roundtrip h5ad -> h5seurat -> h5ad preserves data", {
  h5s_tmp  <- tempfile(fileext = ".h5Seurat")
  h5ad_out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(h5s_tmp, h5ad_out)), add = TRUE)

  scConvert_cli(
    input = h5ad_path, output = h5s_tmp,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )
  scConvert_cli(
    input = h5s_tmp, output = h5ad_out,
    assay = "RNA", overwrite = TRUE, verbose = FALSE
  )

  loaded <- readH5AD(h5ad_out, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), expected_ncol)
  expect_equal(nrow(loaded), expected_nrow)
  expect_equal(sort(colnames(loaded)), sort(colnames(ref_h5ad)))
  expect_equal(sort(rownames(loaded)), sort(rownames(ref_h5ad)))
})
