# Tests for H5MU multimodal support

library(scConvert)

test_that("readH5MU function exists and is exported", {
  expect_true(exists("readH5MU"))
  expect_true(is.function(readH5MU))
})

test_that("writeH5MU function exists and is exported", {
  expect_true(exists("writeH5MU"))
  expect_true(is.function(writeH5MU))
})

test_that("readH5MU validates file existence", {
  expect_error(
    readH5MU("nonexistent_file.h5mu"),
    "File not found"
  )
})

test_that("writeH5MU validates input object", {
  expect_error(
    writeH5MU("not_a_seurat", "output.h5mu"),
    "must be a Seurat object"
  )
})

test_that("writeH5MU works natively without MuDataSeurat", {
  counts <- matrix(rpois(12, 5), nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error(writeH5MU(obj, tmp, verbose = FALSE))
  expect_true(file.exists(tmp))
})

test_that("readH5MU works natively without MuDataSeurat", {
  counts <- matrix(rpois(12, 5), nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  writeH5MU(obj, tmp, verbose = FALSE)
  loaded <- readH5MU(tmp, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), 4)
  expect_equal(nrow(loaded), 3)
})

test_that("GetDefaultModalityMapping returns correct mappings", {
  mapping <- scConvert:::GetDefaultModalityMapping(c("rna", "prot", "atac"))
  expect_equal(mapping[["rna"]], "RNA")
  expect_equal(mapping[["prot"]], "ADT")
  expect_equal(mapping[["atac"]], "ATAC")
})

test_that("GetDefaultModalityMapping preserves unknown modalities", {
  mapping <- scConvert:::GetDefaultModalityMapping(c("rna", "custom_assay"))
  expect_equal(mapping[["rna"]], "RNA")
  expect_equal(mapping[["custom_assay"]], "custom_assay")
})

test_that("GetDefaultAssayToModalityMapping returns correct reverse mappings", {
  mapping <- scConvert:::GetDefaultAssayToModalityMapping(c("RNA", "ADT", "ATAC"))
  expect_equal(mapping[["RNA"]], "rna")
  expect_equal(mapping[["ADT"]], "prot")
  expect_equal(mapping[["ATAC"]], "atac")
})

test_that("GetDefaultAssayToModalityMapping handles conflicts", {
  # RNA and SCT both map to "rna" by default - should resolve
  mapping <- scConvert:::GetDefaultAssayToModalityMapping(c("RNA", "SCT"))
  expect_false(anyDuplicated(mapping) > 0)
  expect_equal(mapping[["RNA"]], "rna")
  expect_equal(mapping[["SCT"]], "sct")
})

test_that("ValidateMultimodalObject accepts valid objects", {
  skip_if_not_installed("Seurat")

  rna <- matrix(rpois(100, 1), 10, 10)
  rownames(rna) <- paste0("Gene", 1:10)
  colnames(rna) <- paste0("Cell", 1:10)
  adt <- matrix(rpois(50, 5), 5, 10)
  rownames(adt) <- paste0("ADT", 1:5)
  colnames(adt) <- colnames(rna)

  obj <- Seurat::CreateSeuratObject(counts = rna)
  obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)

  result <- scConvert:::ValidateMultimodalObject(obj, c("RNA", "ADT"))
  expect_true(result$valid)
})

test_that("writeH5MU rejects overwrite when file exists", {
  skip_if_not_installed("Seurat")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5mu")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    writeH5MU(obj, tmp, overwrite = FALSE),
    "already exists"
  )
})

test_that("Single assay works in writeH5MU", {
  skip_if_not_installed("Seurat")

  counts <- matrix(rpois(12, 5), nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  # Should work fine with a single assay
  expect_no_error(writeH5MU(obj, tmp, verbose = FALSE))
  expect_true(file.exists(tmp))
})

test_that("Multimodal round-trip via writeH5MU/readH5MU preserves structure", {
  skip_if_not_installed("Seurat")

  # Create multimodal object (RNA + ADT)
  set.seed(42)
  rna <- matrix(rpois(200, 5), 20, 10)
  rownames(rna) <- paste0("Gene", 1:20)
  colnames(rna) <- paste0("Cell", 1:10)
  adt <- matrix(rpois(50, 10), 5, 10)
  rownames(adt) <- paste0("ADT", 1:5)
  colnames(adt) <- colnames(rna)

  obj <- Seurat::CreateSeuratObject(counts = rna)
  obj[["ADT"]] <- Seurat::CreateAssayObject(counts = adt)
  obj$group <- factor(sample(c("A", "B"), 10, replace = TRUE))

  tmp <- tempfile(fileext = ".h5mu")
  on.exit(unlink(tmp), add = TRUE)

  # Save (native writer)
  writeH5MU(obj, tmp, overwrite = TRUE, verbose = FALSE)
  expect_true(file.exists(tmp))

  # Load back (native reader)
  loaded <- readH5MU(tmp, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), 10)

  # Check assays are present
  loaded_assays <- Seurat::Assays(loaded)
  expect_true(length(loaded_assays) >= 2)
  expect_true("RNA" %in% loaded_assays)
  expect_true("ADT" %in% loaded_assays)

  # Check data preservation
  orig_counts <- as.matrix(Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"))
  loaded_counts <- as.matrix(Seurat::GetAssayData(loaded, assay = "RNA", layer = "counts"))
  expect_equal(dim(orig_counts), dim(loaded_counts))
  expect_true(all.equal(orig_counts, loaded_counts, check.attributes = FALSE))

  # Check metadata preserved
  expect_true("group" %in% colnames(loaded[[]]))
})
