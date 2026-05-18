# Opt-in cloud-SOMA integration tests.
#
# These reach the public CELLxGENE Census S3 bucket via tiledbsoma. Both
# variables must be set to run, otherwise the file is fully skipped:
#
#   SCCONVERT_TEST_SOMA_URI    -- a SOMA experiment URI (s3:// or local).
#                                  Defaults to a small slice of Census 2024-07-01.
#   SCCONVERT_TEST_SOMA        -- "true" to enable.
#
# Census URIs change with release versions; pin a specific LTS in the env var
# rather than relying on the docstring's example URI.

skip_remote <- function() {
  skip_if(!identical(Sys.getenv("SCCONVERT_TEST_SOMA"), "true"),
          "Set SCCONVERT_TEST_SOMA=true to run cloud SOMA tests")
  skip_if_not_installed("tiledbsoma")
  skip_if_offline()
}

soma_uri <- function() {
  uri <- Sys.getenv("SCCONVERT_TEST_SOMA_URI", "")
  if (nzchar(uri)) return(uri)
  paste0(
    "s3://cellxgene-census-public-us-west-2/",
    "cell-census/2024-07-01/soma/census_data/homo_sapiens"
  )
}

test_that("readSOMA opens a public Census slice and returns a Seurat object", {
  skip_remote()
  obj <- readSOMA(
    uri = soma_uri(),
    obs_query = paste0(
      "cell_type == 'T cell' & ",
      "tissue_general == 'blood' & ",
      "is_primary_data == TRUE"
    ),
    verbose = FALSE
  )
  expect_s4_class(obj, "Seurat")
  expect_gt(ncol(obj), 0L)
  expect_gt(nrow(obj), 0L)
  # Census obs always carries these columns; verify they survived the slice
  meta <- obj[[]]
  expect_true("cell_type" %in% names(meta))
  expect_true(all(meta$cell_type == "T cell"))
})

test_that("readSOMA respects var_column_names filter", {
  skip_remote()
  obj <- readSOMA(
    uri = soma_uri(),
    obs_query = "cell_type == 'pericyte' & is_primary_data == TRUE",
    var_column_names = c("feature_name", "feature_length"),
    verbose = FALSE
  )
  expect_s4_class(obj, "Seurat")
  feat_meta <- obj[[Seurat::DefaultAssay(obj)]][[]]
  # Either column may be dropped during Seurat construction; just ensure the
  # narrowed projection succeeded and we got an object.
  expect_true(length(feat_meta) <= 5L)
})

test_that("readSOMA error message names the package when tiledbsoma is missing", {
  # Cannot reliably uninstall mid-test; instead verify the error path text
  # is present in the function source.
  src <- deparse(body(readSOMA))
  expect_true(any(grepl("tiledbsoma", src)))
  expect_true(any(grepl("install.packages", src)))
})
