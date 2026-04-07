library(scConvert)
library(Seurat)
library(Matrix)

test_that("varp roundtrip via h5ad preserves pairwise variable annotations", {

  # Create a small Seurat object
  set.seed(42)
  ngenes <- 50
  ncells <- 100
  mat <- sparseMatrix(
    i = sample.int(ngenes, 200, replace = TRUE),
    j = sample.int(ncells, 200, replace = TRUE),
    x = as.double(rpois(200, 5) + 1),
    dims = c(ngenes, ncells),
    dimnames = list(paste0("Gene", seq_len(ngenes)),
                    paste0("Cell", seq_len(ncells)))
  )
  obj <- CreateSeuratObject(counts = mat, project = "varp_test")

  # Add a varp matrix (simulating a gene-gene correlation from SCENIC/scvi)
  varp_corr <- matrix(rnorm(ngenes * ngenes), nrow = ngenes, ncol = ngenes)
  varp_corr <- (varp_corr + t(varp_corr)) / 2  # symmetrize
  diag(varp_corr) <- 1
  rownames(varp_corr) <- paste0("Gene", seq_len(ngenes))
  colnames(varp_corr) <- paste0("Gene", seq_len(ngenes))

  # Also test a sparse varp matrix
  varp_sparse <- sparseMatrix(
    i = sample.int(ngenes, 100, replace = TRUE),
    j = sample.int(ngenes, 100, replace = TRUE),
    x = runif(100),
    dims = c(ngenes, ngenes),
    dimnames = list(paste0("Gene", seq_len(ngenes)),
                    paste0("Gene", seq_len(ngenes)))
  )

  obj@misc[["__varp__"]] <- list(
    gene_correlation = varp_corr,
    regulatory_network = varp_sparse
  )

  # Write to h5ad
  h5ad_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)
  writeH5AD(obj, h5ad_path, overwrite = TRUE, verbose = FALSE)

  # Verify varp group exists in h5ad file
  h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
  expect_true(h5$exists("varp"))
  expect_true("gene_correlation" %in% names(h5[["varp"]]))
  h5$close_all()

  # Read back
  obj_rt <- readH5AD(h5ad_path, verbose = FALSE)

  # Check varp is preserved
  expect_false(is.null(obj_rt@misc[["__varp__"]]))
  expect_true("gene_correlation" %in% names(obj_rt@misc[["__varp__"]]))
  expect_true("regulatory_network" %in% names(obj_rt@misc[["__varp__"]]))

  # Check dense matrix values match
  rt_corr <- obj_rt@misc[["__varp__"]][["gene_correlation"]]
  expect_equal(dim(rt_corr), dim(varp_corr))
  expect_equal(as.numeric(rt_corr), as.numeric(varp_corr), tolerance = 1e-6)

  # Check sparse matrix dimensions
  rt_sparse <- obj_rt@misc[["__varp__"]][["regulatory_network"]]
  expect_equal(dim(rt_sparse), dim(varp_sparse))
})

test_that("varp is not written when not present", {

  set.seed(42)
  mat <- sparseMatrix(
    i = sample.int(20, 50, replace = TRUE),
    j = sample.int(30, 50, replace = TRUE),
    x = as.double(rpois(50, 5) + 1),
    dims = c(20, 30),
    dimnames = list(paste0("Gene", 1:20), paste0("Cell", 1:30))
  )
  obj <- CreateSeuratObject(counts = mat)

  h5ad_path <- tempfile(fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)
  writeH5AD(obj, h5ad_path, overwrite = TRUE, verbose = FALSE)

  # varp group should either not exist or be empty
  h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
  if (h5$exists("varp")) {
    expect_equal(length(names(h5[["varp"]])), 0)
  } else {
    expect_false(h5$exists("varp"))
  }
  h5$close_all()
})

# Round-trip varp through the C CLI binary: h5ad (with varp) -> h5seurat
# (misc/__varp__) -> h5ad. Closes the manuscript limitation
# "CLI does not preserve varp".
test_that("varp survives a CLI h5ad -> h5seurat -> h5ad round-trip", {
  cli_path <- file.path("..", "..", "..", "src", "scconvert")
  if (!file.exists(cli_path)) {
    # Also try relative to the package source tree during devtools::test.
    cli_path <- system.file("../src/scconvert", package = "scConvert")
  }
  if (!nzchar(cli_path) || !file.exists(cli_path)) {
    skip("scconvert CLI binary not built; run `cd src && make`")
  }

  set.seed(7)
  ngenes <- 30L
  ncells <- 40L
  mat <- sparseMatrix(
    i = sample.int(ngenes, 120, replace = TRUE),
    j = sample.int(ncells, 120, replace = TRUE),
    x = as.double(rpois(120, 5) + 1),
    dims = c(ngenes, ncells),
    dimnames = list(paste0("Gene", seq_len(ngenes)),
                    paste0("Cell", seq_len(ncells)))
  )
  obj <- CreateSeuratObject(counts = mat, project = "varp_cli")

  # Build a symmetric dense varp (gene x gene) with deterministic values.
  varp_mat <- matrix(seq_len(ngenes * ngenes), nrow = ngenes, ncol = ngenes)
  varp_mat <- (varp_mat + t(varp_mat)) / 2
  dimnames(varp_mat) <- list(rownames(mat), rownames(mat))
  obj@misc[["__varp__"]] <- list(gene_corr = varp_mat)

  h5ad_in  <- tempfile(fileext = ".h5ad")
  h5s_mid  <- tempfile(fileext = ".h5seurat")
  h5ad_out <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(h5ad_in, h5s_mid, h5ad_out)), add = TRUE)

  writeH5AD(obj, h5ad_in, overwrite = TRUE, verbose = FALSE)

  # Step 1: h5ad -> h5seurat via CLI.
  rc1 <- system2(cli_path, c(h5ad_in, h5s_mid, "--overwrite", "--quiet"),
                  stdout = FALSE, stderr = FALSE)
  expect_equal(rc1, 0, info = "CLI h5ad -> h5seurat should succeed")

  # Verify misc/__varp__ group was written by the CLI.
  h5 <- hdf5r::H5File$new(h5s_mid, mode = "r")
  expect_true(h5$exists("misc/__varp__"),
               info = "CLI must write misc/__varp__ in h5seurat")
  expect_true("gene_corr" %in% names(h5[["misc/__varp__"]]))
  h5$close_all()

  # Step 2: h5seurat -> h5ad via CLI.
  rc2 <- system2(cli_path, c(h5s_mid, h5ad_out, "--overwrite", "--quiet"),
                  stdout = FALSE, stderr = FALSE)
  expect_equal(rc2, 0, info = "CLI h5seurat -> h5ad should succeed")

  # Verify varp group is populated in the output h5ad.
  h5 <- hdf5r::H5File$new(h5ad_out, mode = "r")
  expect_true(h5$exists("varp"))
  expect_true("gene_corr" %in% names(h5[["varp"]]))
  h5$close_all()

  # Step 3: load the final h5ad via R and check the values round-trip.
  obj_rt <- readH5AD(h5ad_out, verbose = FALSE)
  expect_false(is.null(obj_rt@misc[["__varp__"]]))
  rt_mat <- obj_rt@misc[["__varp__"]][["gene_corr"]]
  expect_equal(dim(rt_mat), dim(varp_mat))
  expect_equal(as.numeric(rt_mat), as.numeric(varp_mat), tolerance = 1e-9)
})
