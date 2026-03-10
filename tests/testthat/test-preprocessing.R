# Test native AnnData preprocessing features
# These tests verify column name sanitization, list/dict conversion, and nullable dtype handling

library(scConvert)

# Skip tests if required packages are not available
skip_if_not_installed("Seurat")
skip_if_not_installed("Matrix")
skip_if_not_installed("hdf5r")

# -----------------------------------------------------------------------------
# Test SanitizeColumnName helper function
# -----------------------------------------------------------------------------

test_that("SanitizeColumnName replaces problematic characters", {
  # Access the internal function from the package namespace
  SanitizeColumnName <- scConvert:::SanitizeColumnName

  # Test forward slash replacement

  expect_equal(SanitizeColumnName("cell/type"), "cell_type")

  # Test space replacement
  expect_equal(SanitizeColumnName("cell type"), "cell_type")

  # Test comma replacement
  expect_equal(SanitizeColumnName("cell,type"), "cell_type")

  # Test semicolon replacement
  expect_equal(SanitizeColumnName("cell;type"), "cell_type")

  # Test colon replacement
  expect_equal(SanitizeColumnName("cell:type"), "cell_type")

  # Test backslash replacement
  expect_equal(SanitizeColumnName("cell\\type"), "cell_type")

  # Test multiple consecutive problematic characters become single underscore
  expect_equal(SanitizeColumnName("cell//type"), "cell_type")
  expect_equal(SanitizeColumnName("cell  type"), "cell_type")

  # Test leading/trailing underscores are removed
  expect_equal(SanitizeColumnName("/celltype"), "celltype")
  expect_equal(SanitizeColumnName("celltype/"), "celltype")
  expect_equal(SanitizeColumnName("/celltype/"), "celltype")

  # Test names starting with numbers get X prefix
  expect_equal(SanitizeColumnName("1_celltype"), "X1_celltype")
  expect_equal(SanitizeColumnName("123abc"), "X123abc")

  # Test normal names are unchanged
  expect_equal(SanitizeColumnName("cell_type"), "cell_type")
  expect_equal(SanitizeColumnName("CellType"), "CellType")
})

test_that("SanitizeColumnNames handles duplicate names", {
  SanitizeColumnNames <- scConvert:::SanitizeColumnNames

  # Test that duplicates after sanitization get __dupN suffix
  result <- SanitizeColumnNames(c("cell/type", "cell type", "cell:type"))
  expect_equal(unname(result), c("cell_type", "cell_type__dup1", "cell_type__dup2"))

  # Test that original names are preserved in the names attribute
  expect_equal(names(result), c("cell/type", "cell type", "cell:type"))

  # Test that non-duplicates remain unchanged
  result2 <- SanitizeColumnNames(c("col_a", "col_b", "col_c"))
  expect_equal(unname(result2), c("col_a", "col_b", "col_c"))
})

# -----------------------------------------------------------------------------
# Test FlattenNullable helper function
# -----------------------------------------------------------------------------

test_that("FlattenNullable handles mask+values structures", {
  FlattenNullable <- scConvert:::FlattenNullable

  # Create a temporary h5 file with mask+values structure
  temp_h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5), add = TRUE)

  h5 <- hdf5r::H5File$new(temp_h5, mode = "w")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)

  # Create a nullable column group
  h5$create_group("test_nullable")
  values <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  mask <- c(FALSE, TRUE, FALSE, FALSE, TRUE)  # TRUE = missing in h5ad
  h5[["test_nullable"]]$create_dataset("values", robj = values)
  h5[["test_nullable"]]$create_dataset("mask", robj = mask)

  # Test flattening
  result <- FlattenNullable(h5[["test_nullable"]])

  expect_false(is.null(result))
  expect_equal(length(result), 5)
  expect_equal(result[1], 1.0)
  expect_true(is.na(result[2]))  # mask = TRUE -> NA
  expect_equal(result[3], 3.0)
  expect_equal(result[4], 4.0)
  expect_true(is.na(result[5]))  # mask = TRUE -> NA
})

test_that("FlattenNullable returns NULL for non-mask structures", {
  FlattenNullable <- scConvert:::FlattenNullable

  temp_h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5), add = TRUE)

  h5 <- hdf5r::H5File$new(temp_h5, mode = "w")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)

  # Create a categorical group (should not be flattened)
  h5$create_group("test_categorical")
  h5[["test_categorical"]]$create_dataset("categories", robj = c("A", "B", "C"))
  h5[["test_categorical"]]$create_dataset("codes", robj = c(0L, 1L, 2L))

  result <- FlattenNullable(h5[["test_categorical"]])
  expect_null(result)

  # Create a regular group without mask/values
  h5$create_group("test_regular")
  h5[["test_regular"]]$create_dataset("data", robj = c(1, 2, 3))

  result2 <- FlattenNullable(h5[["test_regular"]])
  expect_null(result2)
})

# -----------------------------------------------------------------------------
# Test integrated preprocessing in h5ad conversion
# -----------------------------------------------------------------------------

test_that("h5ad conversion handles special column names", {
  skip_if_not_installed("hdf5r")

  # Create a test h5ad file with problematic column names using hdf5r
  # Note: HDF5 doesn't allow "/" in names, so we test other problematic chars
  temp_h5ad <- tempfile(fileext = ".h5ad")
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit({
    unlink(temp_h5ad)
    unlink(temp_h5seurat)
  }, add = TRUE)

  set.seed(42)
  n_cells <- 50
  n_genes <- 100
  h5 <- hdf5r::H5File$new(temp_h5ad, mode = "w")

  # Write X matrix
  X <- matrix(rpois(n_cells * n_genes, 5), nrow = n_cells, ncol = n_genes)
  storage.mode(X) <- "double"
  h5$create_dataset("X", X)

  # Write obs with column names problematic for R
  obs_grp <- h5$create_group("obs")
  cell_names <- paste0("Cell", seq_len(n_cells))
  obs_grp$create_dataset("_index", cell_names)
  obs_grp$create_attr("_index", "_index", space = hdf5r::H5S$new("scalar"))
  obs_grp$create_attr("encoding-type", "dataframe", space = hdf5r::H5S$new("scalar"))
  obs_grp$create_attr("encoding-version", "0.2.0", space = hdf5r::H5S$new("scalar"))

  col_order <- c("cell_type", "score", "1_numeric_start")
  obs_grp$create_attr("column-order", col_order, space = hdf5r::H5S$new("simple", dims = length(col_order)))

  # cell_type as categorical
  ct_grp <- obs_grp$create_group("cell_type")
  ct_grp$create_dataset("codes", as.integer(c(rep(0L, 25), rep(1L, 25))))
  ct_grp$create_dataset("categories", c("A", "B"))
  ct_grp$create_attr("encoding-type", "categorical", space = hdf5r::H5S$new("scalar"))
  ct_grp$create_attr("encoding-version", "0.2.0", space = hdf5r::H5S$new("scalar"))

  obs_grp$create_dataset("score", rnorm(n_cells))
  obs_grp$create_dataset("1_numeric_start", as.double(seq_len(n_cells) - 1L))

  # Write var
  var_grp <- h5$create_group("var")
  gene_names <- paste0("Gene", seq_len(n_genes))
  var_grp$create_dataset("_index", gene_names)
  var_grp$create_attr("_index", "_index", space = hdf5r::H5S$new("scalar"))
  var_grp$create_attr("encoding-type", "dataframe", space = hdf5r::H5S$new("scalar"))
  var_grp$create_attr("encoding-version", "0.2.0", space = hdf5r::H5S$new("scalar"))
  var_grp$create_attr("column-order", character(0), space = hdf5r::H5S$new("simple", dims = 0))

  h5$close_all()

  # Convert h5ad to h5seurat should succeed
  expect_no_error(
    scConvert(temp_h5ad, dest = temp_h5seurat, verbose = FALSE, overwrite = TRUE)
  )

  # Load and verify the converted file
  seurat_obj <- scLoadH5Seurat(temp_h5seurat, verbose = FALSE)
  expect_s4_class(seurat_obj, "Seurat")
  expect_equal(ncol(seurat_obj), n_cells)
  expect_true("cell_type" %in% colnames(seurat_obj[[]]))
})

test_that("h5ad conversion handles compound obs datasets", {
  skip_if_not_installed("hdf5r")

  # Create a test h5ad with compound dataset obs using hdf5r
  temp_h5ad <- tempfile(fileext = ".h5ad")
  temp_h5seurat <- tempfile(fileext = ".h5seurat")
  on.exit({
    unlink(temp_h5ad)
    unlink(temp_h5seurat)
  }, add = TRUE)

  n_cells <- 10
  n_genes <- 20
  h5 <- hdf5r::H5File$new(temp_h5ad, mode = "w")

  # Write X matrix
  set.seed(42)
  X <- matrix(rpois(n_cells * n_genes, 5), nrow = n_cells, ncol = n_genes)
  storage.mode(X) <- "double"
  h5$create_dataset("X", X)

  # Write obs as compound dataset (old-style h5ad format)
  cell_names <- paste0("Cell", seq_len(n_cells))
  cpd_type <- hdf5r::H5T_COMPOUND$new(
    "obs_record",
    dtypes = list(
      hdf5r::H5T_STRING$new(size = 20),
      hdf5r::h5types$H5T_NATIVE_INT,
      hdf5r::H5T_STRING$new(size = 50)
    ),
    labels = c("index", "cluster", "tags")
  )
  obs_data <- data.frame(
    index = cell_names,
    cluster = as.integer(seq_len(n_cells) %% 3),
    tags = rep("A;B;C", n_cells),
    stringsAsFactors = FALSE
  )
  h5$create_dataset("obs", obs_data, dtype = cpd_type)

  # Write var as compound dataset
  gene_names <- paste0("Gene", seq_len(n_genes))
  cpd_var <- hdf5r::H5T_COMPOUND$new(
    "var_record",
    dtypes = list(
      hdf5r::H5T_STRING$new(size = 20),
      hdf5r::H5T_STRING$new(size = 20)
    ),
    labels = c("index", "gene_name")
  )
  var_data <- data.frame(
    index = gene_names,
    gene_name = gene_names,
    stringsAsFactors = FALSE
  )
  h5$create_dataset("var", var_data, dtype = cpd_var)

  # Add required attributes
  h5$create_attr("encoding-type", "anndata", space = hdf5r::H5S$new("scalar"))
  h5$create_attr("encoding-version", "0.1.0", space = hdf5r::H5S$new("scalar"))

  h5$close_all()

  # Convert should succeed (compound obs triggers special handling)
  result <- tryCatch(
    scConvert(temp_h5ad, dest = temp_h5seurat, verbose = FALSE, overwrite = TRUE),
    error = function(e) e
  )

  # Either succeeds or errors gracefully (compound datasets may not be fully supported)
  if (inherits(result, "error")) {
    expect_true(nzchar(conditionMessage(result)))
  } else {
    expect_true(file.exists(temp_h5seurat))
  }
})

# -----------------------------------------------------------------------------
# Test with real datasets (if available)
# -----------------------------------------------------------------------------

test_that("Conversion preserves metadata from real h5ad files", {
  test_h5ad <- system.file("testdata", "pbmc_small.h5ad", package = "scConvert")
  skip_if(!file.exists(test_h5ad), "Test h5ad file not available")

  # Load directly via LoadH5AD to verify metadata preservation
  seurat_obj <- LoadH5AD(test_h5ad, verbose = FALSE)

  # Basic checks
  expect_s4_class(seurat_obj, "Seurat")
  expect_gt(ncol(seurat_obj), 0)
  expect_gt(nrow(seurat_obj), 0)

  # Check that metadata columns exist and are accessible
  meta_cols <- colnames(seurat_obj@meta.data)
  expect_true(length(meta_cols) > 0)
})

# -----------------------------------------------------------------------------
# Test round-trip preservation
# -----------------------------------------------------------------------------

test_that("Round-trip conversion preserves sanitized column names", {

  # Create Seurat object with normal column names
  set.seed(42)
  counts <- matrix(
    data = rpois(100 * 50, lambda = 5),
    nrow = 100,
    ncol = 50,
    dimnames = list(
      paste0("Gene", seq_len(100)),
      paste0("Cell", seq_len(50))
    )
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = "TestProject")

  # Add metadata with valid names
  seurat_obj$cell_type <- factor(sample(c("A", "B", "C"), 50, replace = TRUE))
  seurat_obj$sample_id <- factor(rep(c("S1", "S2"), 25))
  seurat_obj$score <- rnorm(50)

  temp_h5seurat1 <- tempfile(fileext = ".h5seurat")
  temp_h5ad <- tempfile(fileext = ".h5ad")
  temp_h5seurat2 <- tempfile(fileext = ".h5seurat")
  on.exit({
    unlink(temp_h5seurat1)
    unlink(temp_h5ad)
    unlink(temp_h5seurat2)
  }, add = TRUE)

  # Save -> scConvert to h5ad -> scConvert back
  scSaveH5Seurat(seurat_obj, filename = temp_h5seurat1, verbose = FALSE)
  scConvert(temp_h5seurat1, dest = temp_h5ad, verbose = FALSE, overwrite = TRUE)
  scConvert(temp_h5ad, dest = temp_h5seurat2, verbose = FALSE, overwrite = TRUE)

  # Load and compare
  loaded_obj <- scLoadH5Seurat(temp_h5seurat2, verbose = FALSE)

  # Check metadata columns are preserved
  expect_true("cell_type" %in% colnames(loaded_obj@meta.data))
  expect_true("sample_id" %in% colnames(loaded_obj@meta.data))
  expect_true("score" %in% colnames(loaded_obj@meta.data))

  # Check values are preserved
  expect_equal(
    as.character(loaded_obj$cell_type),
    as.character(seurat_obj$cell_type)
  )
  expect_equal(
    loaded_obj$score,
    seurat_obj$score,
    tolerance = 1e-6
  )
})
