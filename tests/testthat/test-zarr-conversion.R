# Tests for Zarr-based AnnData conversion (LoadZarr / SaveZarr)

library(scConvert)

# All zarr tests require jsonlite for metadata I/O
skip_if_not_installed("jsonlite")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper: Create synthetic zarr v2 AnnData stores for testing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a minimal zarr v2 AnnData store with dense matrix
#' @keywords internal
create_test_zarr_dense <- function(dir, n_cells = 20, n_genes = 10, seed = 42) {
  set.seed(seed)
  counts <- matrix(rpois(n_cells * n_genes, lambda = 5),
                   nrow = n_cells, ncol = n_genes)
  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))

  compressor <- list(id = "zlib", level = 1L)

  # Root group
  scConvert:::.zarr_create_group(dir, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # X: dense matrix (cells x genes)
  scConvert:::.zarr_write_numeric(
    dir = file.path(dir, "X"),
    data = counts,
    dtype = "<f8",
    compressor = compressor,
    attrs = list(
      `encoding-type` = "array",
      `encoding-version` = "0.2.0"
    )
  )

  # obs
  scConvert:::.zarr_create_group(file.path(dir, "obs"), attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = list()
  ))
  scConvert:::.zarr_write_strings(
    dir = file.path(dir, "obs", "_index"),
    strings = cell_names,
    compressor = compressor,
    attrs = list(`encoding-type` = "string-array", `encoding-version` = "0.2.0")
  )

  # var
  scConvert:::.zarr_create_group(file.path(dir, "var"), attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = list()
  ))
  scConvert:::.zarr_write_strings(
    dir = file.path(dir, "var", "_index"),
    strings = gene_names,
    compressor = compressor,
    attrs = list(`encoding-type` = "string-array", `encoding-version` = "0.2.0")
  )

  # Empty groups for layers, obsm, obsp, uns
  scConvert:::.zarr_create_group(file.path(dir, "layers"))
  scConvert:::.zarr_create_group(file.path(dir, "obsm"))
  scConvert:::.zarr_create_group(file.path(dir, "obsp"))
  scConvert:::.zarr_create_group(file.path(dir, "uns"))

  list(counts = counts, cell_names = cell_names, gene_names = gene_names)
}

#' Create a zarr v2 AnnData store with sparse CSR matrix
#' @keywords internal
create_test_zarr_sparse <- function(dir, n_cells = 15, n_genes = 8, seed = 123) {
  set.seed(seed)
  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))

  # Create sparse matrix (cells x genes for storage)
  mat <- Matrix::rsparsematrix(nrow = n_cells, ncol = n_genes, density = 0.3)
  mat@x <- abs(mat@x) * 10

  # Convert to CSR
  rsp <- as(mat, "RsparseMatrix")

  compressor <- list(id = "zlib", level = 1L)

  # Root
  scConvert:::.zarr_create_group(dir, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # X: sparse CSR
  scConvert:::.zarr_create_group(file.path(dir, "X"), attrs = list(
    `encoding-type` = "csr_matrix",
    `encoding-version` = "0.1.0",
    shape = list(as.integer(n_cells), as.integer(n_genes))
  ))
  scConvert:::.zarr_write_numeric(file.path(dir, "X", "data"),
                                 data = rsp@x, dtype = "<f8", compressor = compressor)
  scConvert:::.zarr_write_numeric(file.path(dir, "X", "indices"),
                                 data = as.integer(rsp@j), dtype = "<i4",
                                 compressor = compressor)
  scConvert:::.zarr_write_numeric(file.path(dir, "X", "indptr"),
                                 data = as.integer(rsp@p), dtype = "<i4",
                                 compressor = compressor)

  # obs
  scConvert:::.zarr_create_group(file.path(dir, "obs"), attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = list()
  ))
  scConvert:::.zarr_write_strings(file.path(dir, "obs", "_index"),
                                 strings = cell_names, compressor = compressor,
                                 attrs = list(`encoding-type` = "string-array",
                                              `encoding-version` = "0.2.0"))

  # var
  scConvert:::.zarr_create_group(file.path(dir, "var"), attrs = list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = list()
  ))
  scConvert:::.zarr_write_strings(file.path(dir, "var", "_index"),
                                 strings = gene_names, compressor = compressor,
                                 attrs = list(`encoding-type` = "string-array",
                                              `encoding-version` = "0.2.0"))

  scConvert:::.zarr_create_group(file.path(dir, "layers"))
  scConvert:::.zarr_create_group(file.path(dir, "obsm"))
  scConvert:::.zarr_create_group(file.path(dir, "obsp"))
  scConvert:::.zarr_create_group(file.path(dir, "uns"))

  list(mat = mat, cell_names = cell_names, gene_names = gene_names)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AnnDataEncoding.R unit tests
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("DecodeCategorical handles basic encoding", {
  codes <- c(0L, 1L, 2L, 0L, 1L)
  categories <- c("A", "B", "C")
  result <- scConvert:::DecodeCategorical(codes, categories)
  expect_s3_class(result, "factor")
  expect_equal(levels(result), categories)
  expect_equal(as.character(result), c("A", "B", "C", "A", "B"))
})

test_that("DecodeCategorical handles -1 as NA", {
  codes <- c(0L, -1L, 1L)
  categories <- c("X", "Y")
  result <- scConvert:::DecodeCategorical(codes, categories)
  expect_true(is.na(result[2]))
  expect_equal(as.character(result[c(1, 3)]), c("X", "Y"))
})

test_that("EncodeCategorical roundtrips with DecodeCategorical", {
  original <- factor(c("cat", "dog", "cat", NA, "bird"), levels = c("bird", "cat", "dog"))
  encoded <- scConvert:::EncodeCategorical(original)
  decoded <- scConvert:::DecodeCategorical(encoded$codes, encoded$categories)
  expect_equal(levels(decoded), levels(original))
  expect_equal(as.character(decoded), as.character(original))
  expect_equal(is.na(decoded), is.na(original))
})

test_that("ReconstructSparseCSR creates correct dgCMatrix", {
  # Simple 3x3 identity matrix in CSR format (3x3, row-compressed)
  data <- c(1.0, 1.0, 1.0)
  indices <- c(0L, 1L, 2L)
  indptr <- c(0L, 1L, 2L, 3L)
  shape <- c(3L, 3L)

  # Without transpose
  mat <- scConvert:::ReconstructSparseCSR(data, indices, indptr, shape, transpose = FALSE)
  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(3, 3))
  expect_equal(sum(mat), 3)

  # With transpose
  mat_t <- scConvert:::ReconstructSparseCSR(data, indices, indptr, shape, transpose = TRUE)
  expect_equal(dim(mat_t), c(3, 3))
})

test_that("DeconstructSparseCSR inverts ReconstructSparseCSR", {
  set.seed(42)
  # Create a genes x cells sparse matrix (5 genes, 10 cells)
  original <- Matrix::rsparsematrix(5, 10, density = 0.4)
  original@x <- abs(original@x)

  csr <- scConvert:::DeconstructSparseCSR(original)
  expect_equal(csr$shape, c(10L, 5L))  # cells x genes

  # Reconstruct (with transpose = TRUE to get back genes x cells)
  reconstructed <- scConvert:::ReconstructSparseCSR(
    csr$data, csr$indices, csr$indptr, csr$shape, transpose = TRUE
  )
  expect_equal(dim(reconstructed), dim(original))
  expect_equal(sum(abs(reconstructed - original)), 0)
})

test_that("AnnDataReductionKey maps correctly", {
  expect_equal(scConvert:::AnnDataReductionKey("pca"), "PC_")
  expect_equal(scConvert:::AnnDataReductionKey("umap"), "UMAP_")
  expect_equal(scConvert:::AnnDataReductionKey("tsne"), "tSNE_")
  expect_equal(scConvert:::AnnDataReductionKey("custom"), "CUSTOM_")
})

test_that("vlen-utf8 encode/decode roundtrips", {
  strings <- c("hello", "world", "", "unicode\u00e9")
  encoded <- scConvert:::.encode_vlen_utf8(strings)
  decoded <- scConvert:::.decode_vlen_utf8(encoded, length(strings))
  expect_equal(decoded, strings)
})

test_that("vlen-utf8 handles empty and single strings", {
  # Empty
  encoded <- scConvert:::.encode_vlen_utf8(character(0))
  decoded <- scConvert:::.decode_vlen_utf8(encoded, 0)
  expect_equal(decoded, character(0))

  # Single
  encoded <- scConvert:::.encode_vlen_utf8("test")
  decoded <- scConvert:::.decode_vlen_utf8(encoded, 1)
  expect_equal(decoded, "test")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Zarr I/O low-level tests
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("zarr v2 numeric array write/read roundtrips (1D)", {
  tmp <- tempfile(pattern = "zarr_test_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  values <- c(1.5, 2.5, 3.5, 4.5)
  scConvert:::.zarr_write_numeric(tmp, values, dtype = "<f8",
                                 compressor = list(id = "zlib", level = 1L))
  result <- scConvert:::.zarr_read_numeric(tmp, "")

  # Read from the directory we wrote to - need to adjust path
  result <- scConvert:::.zarr_read_numeric(dirname(tmp), basename(tmp))
  expect_equal(result, values, tolerance = 1e-10)
})

test_that("zarr v2 numeric array write/read roundtrips (2D)", {
  tmp_store <- tempfile(pattern = "zarr_store_")
  on.exit(unlink(tmp_store, recursive = TRUE), add = TRUE)

  scConvert:::.zarr_create_group(tmp_store)
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  scConvert:::.zarr_write_numeric(file.path(tmp_store, "arr"), mat, dtype = "<f8",
                                 compressor = list(id = "zlib", level = 1L))
  result <- scConvert:::.zarr_read_numeric(tmp_store, "arr")

  expect_equal(dim(result), c(3, 4))
  expect_equal(as.vector(result), as.vector(mat))
})

test_that("zarr v2 string array write/read roundtrips", {
  tmp_store <- tempfile(pattern = "zarr_store_")
  on.exit(unlink(tmp_store, recursive = TRUE), add = TRUE)

  scConvert:::.zarr_create_group(tmp_store)
  strings <- c("AAACCTGAGCTAACTC-1", "AAACCTGAGTGGGATC-1", "AAACCTGTCAATCACG-1")
  scConvert:::.zarr_write_strings(file.path(tmp_store, "barcodes"), strings,
                                 compressor = list(id = "zlib", level = 1L))
  result <- scConvert:::.zarr_read_strings(tmp_store, "barcodes")

  expect_equal(result, strings)
})

test_that("zarr v2 integer array roundtrips", {
  tmp_store <- tempfile(pattern = "zarr_store_")
  on.exit(unlink(tmp_store, recursive = TRUE), add = TRUE)

  scConvert:::.zarr_create_group(tmp_store)
  ints <- c(0L, 1L, 2L, 0L, 1L, -1L)
  scConvert:::.zarr_write_numeric(file.path(tmp_store, "codes"), ints, dtype = "<i4",
                                 compressor = list(id = "zlib", level = 1L))
  result <- scConvert:::.zarr_read_numeric(tmp_store, "codes")
  expect_equal(result, ints)
})

test_that("zarr_store_version detects v2 correctly", {
  tmp <- tempfile(pattern = "zarr_v2_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  scConvert:::.zarr_create_group(tmp)
  expect_equal(scConvert:::.zarr_store_version(tmp), 2L)
})

test_that("zarr_node_type identifies groups and arrays", {
  tmp <- tempfile(pattern = "zarr_node_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  scConvert:::.zarr_create_group(tmp)
  scConvert:::.zarr_create_group(file.path(tmp, "grp"))
  scConvert:::.zarr_write_numeric(file.path(tmp, "arr"), c(1.0, 2.0),
                                 compressor = list(id = "zlib", level = 1L))

  expect_equal(scConvert:::.zarr_node_type(tmp, "grp"), "group")
  expect_equal(scConvert:::.zarr_node_type(tmp, "arr"), "array")
  expect_equal(scConvert:::.zarr_node_type(tmp, "nonexistent"), "missing")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LoadZarr tests
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("LoadZarr function exists and is exported", {
  expect_true(exists("LoadZarr"))
  expect_true(is.function(LoadZarr))
})

test_that("LoadZarr validates directory existence", {
  expect_error(LoadZarr("nonexistent_dir.zarr"), "not found")
})

test_that("LoadZarr creates Seurat object from dense zarr store", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  ref <- create_test_zarr_dense(tmp, n_cells = 20, n_genes = 10)

  result <- LoadZarr(tmp, verbose = FALSE)

  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), 20)
  expect_equal(nrow(result), 10)
  expect_identical(colnames(result), ref$cell_names)
  expect_equal(length(rownames(result)), length(ref$gene_names))
})

test_that("LoadZarr reads sparse CSR zarr correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  ref <- create_test_zarr_sparse(tmp)

  result <- LoadZarr(tmp, verbose = FALSE)

  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), 15)
  expect_equal(nrow(result), 8)
  expect_identical(colnames(result), ref$cell_names)
})

test_that("LoadZarr preserves categorical metadata", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  ref <- create_test_zarr_dense(tmp, n_cells = 10, n_genes = 5)
  compressor <- list(id = "zlib", level = 1L)

  # Add categorical column to obs
  categories <- c("TypeA", "TypeB", "TypeC")
  codes <- as.integer(c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0))

  scConvert:::.zarr_create_group(file.path(tmp, "obs", "cluster"), attrs = list(
    `encoding-type` = "categorical",
    `encoding-version` = "0.2.0",
    ordered = FALSE
  ))
  scConvert:::.zarr_write_numeric(file.path(tmp, "obs", "cluster", "codes"),
                                 data = codes, dtype = "<i1", compressor = compressor)
  scConvert:::.zarr_write_strings(file.path(tmp, "obs", "cluster", "categories"),
                                 strings = categories, compressor = compressor,
                                 attrs = list(`encoding-type` = "string-array",
                                              `encoding-version` = "0.2.0"))

  # Update obs column-order
  scConvert:::.zarr_write_attrs(file.path(tmp, "obs"), list(
    `encoding-type` = "dataframe",
    `encoding-version` = "0.2.0",
    `_index` = "_index",
    `column-order` = list("cluster")
  ))

  result <- LoadZarr(tmp, verbose = FALSE)

  expect_true("cluster" %in% colnames(result@meta.data))
  expect_s3_class(result$cluster, "factor")
  expect_equal(levels(result$cluster), categories)
  expect_equal(as.character(result$cluster), categories[codes + 1])
})

test_that("LoadZarr preserves dimensional reductions", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  n_cells <- 10
  ref <- create_test_zarr_dense(tmp, n_cells = n_cells, n_genes = 5)
  compressor <- list(id = "zlib", level = 1L)

  # Add PCA and UMAP to obsm
  set.seed(99)
  pca_embed <- matrix(rnorm(n_cells * 5), nrow = n_cells, ncol = 5)
  scConvert:::.zarr_write_numeric(file.path(tmp, "obsm", "X_pca"),
                                 data = pca_embed, dtype = "<f8",
                                 compressor = compressor)

  umap_embed <- matrix(rnorm(n_cells * 2), nrow = n_cells, ncol = 2)
  scConvert:::.zarr_write_numeric(file.path(tmp, "obsm", "X_umap"),
                                 data = umap_embed, dtype = "<f8",
                                 compressor = compressor)

  result <- LoadZarr(tmp, verbose = FALSE)

  expect_true("pca" %in% names(result@reductions))
  expect_true("umap" %in% names(result@reductions))

  pca <- Seurat::Embeddings(result, "pca")
  expect_equal(nrow(pca), n_cells)
  expect_equal(ncol(pca), 5)

  umap <- Seurat::Embeddings(result, "umap")
  expect_equal(nrow(umap), n_cells)
  expect_equal(ncol(umap), 2)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SaveZarr tests
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("SaveZarr function exists and is exported", {
  expect_true(exists("SaveZarr"))
  expect_true(is.function(SaveZarr))
})

test_that("SaveZarr validates inputs", {
  expect_error(SaveZarr("not_a_seurat", "out.zarr"), "must be a Seurat")
})

test_that("SaveZarr creates valid zarr v2 store", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(42)
  mat <- matrix(rpois(100, 5), nrow = 10, ncol = 10,
                dimnames = list(paste0("Gene", 1:10), paste0("Cell", 1:10)))
  obj <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  SaveZarr(obj, tmp, verbose = FALSE)

  expect_true(dir.exists(tmp))
  expect_true(file.exists(file.path(tmp, ".zgroup")))
  expect_true(file.exists(file.path(tmp, ".zattrs")))
  expect_true(dir.exists(file.path(tmp, "X")))
  expect_true(dir.exists(file.path(tmp, "obs")))
  expect_true(dir.exists(file.path(tmp, "var")))
})

test_that("SaveZarr respects overwrite parameter", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  mat <- matrix(1, nrow = 3, ncol = 3,
                dimnames = list(paste0("G", 1:3), paste0("C", 1:3)))
  obj <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  SaveZarr(obj, tmp, verbose = FALSE)
  expect_error(SaveZarr(obj, tmp, verbose = FALSE), "exists")
  expect_silent(SaveZarr(obj, tmp, overwrite = TRUE, verbose = FALSE))
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Roundtrip tests
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Seurat -> SaveZarr -> LoadZarr roundtrip preserves data", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(42)
  n_genes <- 20
  n_cells <- 50
  mat <- matrix(rpois(n_genes * n_cells, 5), nrow = n_genes, ncol = n_cells,
                dimnames = list(paste0("Gene", 1:n_genes), paste0("Cell", 1:n_cells)))
  original <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  # Add metadata
  original$cell_type <- factor(sample(c("A", "B", "C"), n_cells, replace = TRUE))
  original$score <- runif(n_cells)

  SaveZarr(original, tmp, verbose = FALSE)
  loaded <- LoadZarr(tmp, verbose = FALSE)

  expect_equal(ncol(loaded), ncol(original))
  expect_equal(nrow(loaded), nrow(original))
  expect_equal(colnames(loaded), colnames(original))

  # Check metadata roundtrip
  expect_true("cell_type" %in% colnames(loaded@meta.data))
  expect_s3_class(loaded$cell_type, "factor")
  expect_equal(as.character(loaded$cell_type), as.character(original$cell_type))

  expect_true("score" %in% colnames(loaded@meta.data))
  expect_equal(loaded$score, original$score, tolerance = 1e-6)
})

test_that("Roundtrip preserves sparse matrix values", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(123)
  n_genes <- 30
  n_cells <- 40
  mat <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.2)
  mat@x <- abs(mat@x) * 10
  dimnames(mat) <- list(paste0("Gene", 1:n_genes), paste0("Cell", 1:n_cells))

  original <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  SaveZarr(original, tmp, verbose = FALSE)
  loaded <- LoadZarr(tmp, verbose = FALSE)

  # Compare matrix values
  orig_counts <- Seurat::GetAssayData(original, layer = "counts")
  loaded_counts <- Seurat::GetAssayData(loaded, layer = "counts")

  expect_equal(dim(loaded_counts), dim(orig_counts))
  expect_equal(sum(loaded_counts), sum(orig_counts), tolerance = 1e-6)
})

test_that("Roundtrip preserves dimensional reductions", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(42)
  n_genes <- 20
  n_cells <- 30
  mat <- matrix(rpois(n_genes * n_cells, 5), nrow = n_genes, ncol = n_cells,
                dimnames = list(paste0("Gene", 1:n_genes), paste0("Cell", 1:n_cells)))
  original <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  # Add PCA reduction
  pca_embed <- matrix(rnorm(n_cells * 5), nrow = n_cells, ncol = 5,
                      dimnames = list(colnames(original), paste0("PC_", 1:5)))
  original[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = pca_embed, key = "PC_", assay = "RNA"
  )

  SaveZarr(original, tmp, verbose = FALSE)
  loaded <- LoadZarr(tmp, verbose = FALSE)

  expect_true("pca" %in% names(loaded@reductions))
  loaded_pca <- Seurat::Embeddings(loaded, "pca")
  expect_equal(dim(loaded_pca), c(n_cells, 5))
  expect_equal(as.vector(loaded_pca), as.vector(pca_embed), tolerance = 1e-6)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cross-format conversion tests (via scConvert hub)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("zarr format is registered in FormatRegistry", {
  formats <- scConvert:::ListFormats()
  expect_true("zarr" %in% formats$loaders)
  expect_true("zarr" %in% formats$savers)
})

test_that("scConvert Seurat -> zarr works via hub", {
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(42)
  mat <- matrix(rpois(100, 5), nrow = 10, ncol = 10,
                dimnames = list(paste0("Gene", 1:10), paste0("Cell", 1:10)))
  obj <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  result <- scConvert(obj, dest = tmp, overwrite = TRUE, verbose = FALSE)
  expect_true(dir.exists(tmp))

  # Verify we can load it back
  loaded <- LoadZarr(tmp, verbose = FALSE)
  expect_s4_class(loaded, "Seurat")
  expect_equal(ncol(loaded), 10)
})
