# Regression tests for bugs fixed during benchmarking
#
# Each test here corresponds to a specific bug report documented in the
# manuscript or the development notes. If any of these fail, a previously-
# fixed bug has returned.

library(scConvert)

# ---------------------------------------------------------------------------
# Regression 1: Unsorted CSR column indices (scanpy pbmc3k_processed edge case)
#
# Scipy's CSR/CSC matrices are not guaranteed to have sorted indices within
# a row/column. scanpy datasets such as pbmc3k_processed ship with unsorted
# indices, which previously caused readH5AD() to construct an invalid
# dgCMatrix. The fix sorts indices within each column during construction.
# ---------------------------------------------------------------------------
test_that("readH5AD handles unsorted CSR indices (pbmc3k regression)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  # Build a 5x4 CSR matrix with DELIBERATELY unsorted column indices per row.
  # Row 0: cols 2, 0, 1 (unsorted)   values 3, 1, 2
  # Row 1: cols 3, 1    (unsorted)   values 5, 4
  # Row 2: cols 0       (sorted)     values 6
  # Row 3: (empty)
  # Row 4: cols 2, 3, 0 (unsorted)   values 9, 8, 7
  indices <- as.integer(c(2, 0, 1,   3, 1,   0,   2, 3, 0))
  data    <- as.double(c(3, 1, 2,   5, 4,   6,   9, 8, 7))
  indptr  <- as.integer(c(0, 3, 5, 6, 6, 9))

  n_rows <- 5L   # cells (h5ad CSR orientation)
  n_cols <- 4L   # genes

  cell_names <- paste0("Cell_", seq_len(n_rows))
  gene_names <- paste0("Gene_", seq_len(n_cols))

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  x_group <- h5$create_group("X")
  x_group$create_dataset("data",    robj = data,    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  x_group$create_dataset("indices", robj = indices, dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_group$create_dataset("indptr",  robj = indptr,  dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_group$create_attr("shape", robj = c(n_rows, n_cols),
                      dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_group$create_attr("encoding-type", robj = "csr_matrix",
                      dtype = hdf5r::H5T_STRING$new(size = Inf),
                      space = hdf5r::H5S$new(type = "scalar"))

  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  obs$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  var <- h5$create_group("var")
  var$create_dataset("_index", robj = gene_names)
  var$create_attr("encoding-type", robj = "dataframe",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))

  h5$close_all()

  # Should load without error.
  result <- readH5AD(tmp, verbose = FALSE)
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), n_rows)
  expect_equal(nrow(result), n_cols)

  # Verify numerical content round-trips. The dense view should match what
  # we expect from the unsorted CSR.
  expected <- matrix(0, nrow = n_rows, ncol = n_cols,
                     dimnames = list(cell_names, gene_names))
  expected[1, c(3, 1, 2)] <- c(3, 1, 2)
  expected[2, c(4, 2)]    <- c(5, 4)
  expected[3, 1]          <- 6
  expected[5, c(3, 4, 1)] <- c(9, 8, 7)
  expected_t <- t(expected)  # Seurat stores genes x cells

  mat <- as.matrix(Seurat::GetAssayData(result, layer = "counts"))
  if (all(dim(mat) == 0)) {
    mat <- as.matrix(Seurat::GetAssayData(result, layer = "data"))
  }
  expect_equal(dim(mat), dim(expected_t))
  expect_equal(unname(mat), unname(expected_t), tolerance = 1e-9)

  # Verify the underlying dgCMatrix has sorted row indices per column
  # (this is what breaks if the sort is missing).
  dgc <- as(mat, "CsparseMatrix")
  n_col <- ncol(dgc)
  for (ci in seq_len(n_col)) {
    start <- dgc@p[ci] + 1L
    end   <- dgc@p[ci + 1L]
    if (end >= start) {
      rng <- dgc@i[start:end]
      expect_true(!is.unsorted(rng),
                  info = sprintf("column %d has unsorted row indices", ci))
    }
  }
})

# ---------------------------------------------------------------------------
# Regression 2: H5SeuratToZarr 1D dense dataset crash
#
# An earlier version of H5SeuratToZarr() crashed on 1D dense HDF5 datasets
# because it tried to infer (rows, cols) without enough information.
# The fix is to warn and skip instead of crashing.
# ---------------------------------------------------------------------------
test_that("H5SeuratToZarr handles 1D dense datasets without crashing", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("Seurat")

  # We create a tiny Seurat object and save it to h5Seurat, then add a
  # stray 1D dense dataset to exercise the edge case path.
  set.seed(1)
  n_cells <- 10
  n_genes <- 6
  mat <- matrix(rpois(n_cells * n_genes, 3),
                nrow = n_genes, ncol = n_cells,
                dimnames = list(paste0("g", seq_len(n_genes)),
                                paste0("c", seq_len(n_cells))))
  obj <- Seurat::CreateSeuratObject(counts = mat)

  h5s_path  <- tempfile(fileext = ".h5seurat")
  zarr_path <- tempfile(fileext = ".zarr")
  on.exit({
    unlink(h5s_path)
    unlink(zarr_path, recursive = TRUE)
  }, add = TRUE)

  writeH5Seurat(obj, h5s_path, overwrite = TRUE, verbose = FALSE)

  # Convert should succeed (may emit warnings about any skipped datasets).
  expect_silent_or_warning <- function(expr) {
    suppressWarnings(expr)
  }
  expect_silent_or_warning({
    H5SeuratToZarr(h5s_path, zarr_path,
                   overwrite = TRUE, verbose = FALSE)
  })
  expect_true(dir.exists(zarr_path))

  # And the resulting zarr store should be readable.
  out <- readZarr(zarr_path, verbose = FALSE)
  expect_s4_class(out, "Seurat")
  expect_equal(ncol(out), n_cells)
})

# ---------------------------------------------------------------------------
# Regression 3: writeZarr must skip scale.data (shape != X shape)
#
# Seurat's ScaleData produces an (n_hvg x n_cells) matrix that does not match
# the X shape (n_genes x n_cells). AnnData requires all layers to match X;
# writeZarr must therefore skip scale.data rather than emit a mismatched
# layer that anndata.read_zarr() will reject.
# ---------------------------------------------------------------------------
test_that("writeZarr skips mismatched scale.data layer", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("Seurat")

  set.seed(7)
  n_cells <- 30
  n_genes <- 50
  counts <- matrix(rpois(n_cells * n_genes, 3),
                   nrow = n_genes, ncol = n_cells,
                   dimnames = list(paste0("g", seq_len(n_genes)),
                                   paste0("c", seq_len(n_cells))))
  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 20, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)

  # scale.data should be (20 x 30), NOT matching X's (50 x 30).
  scale_mat <- Seurat::GetAssayData(obj, layer = "scale.data")
  expect_equal(dim(scale_mat), c(20L, n_cells))

  zarr_path <- tempfile(fileext = ".zarr")
  on.exit(unlink(zarr_path, recursive = TRUE), add = TRUE)

  writeZarr(obj, zarr_path, overwrite = TRUE, verbose = FALSE)
  expect_true(dir.exists(zarr_path))

  # Any layer present in the store must have a shape matching X.
  layers_path <- file.path(zarr_path, "layers")
  if (dir.exists(layers_path)) {
    layer_names <- list.dirs(layers_path, recursive = FALSE, full.names = FALSE)
    expect_false("scale.data" %in% layer_names,
                 info = "scale.data layer should be skipped when it does not match X shape")
  }

  # Round-trip should work end-to-end.
  rt <- readZarr(zarr_path, verbose = FALSE)
  expect_s4_class(rt, "Seurat")
  expect_equal(ncol(rt), n_cells)
  expect_equal(nrow(rt), n_genes)
})

# ---------------------------------------------------------------------------
# Regression 4: Linux libhdf5-serial 1.10 vlen segfault in obs attribute rename
#
# sc_stream_obs_h5ad_to_h5seurat previously used a read-deep-copy-write
# pattern (sc_get_str_array_attr + sc_set_str_array_attr) to rename the obs
# dataframe's column-order attribute to h5Seurat's colnames, and its reverse
# did the same in the other direction. On Linux libhdf5-serial 1.10 the vlen
# type converter (H5T__conv_vlen) dereferenced buffer pointers via strlen
# during H5Awrite and segfaulted on every conversion. The fix replaces the
# read-write cycle with HDF5's native H5Arename, which preserves the on-disk
# type/encoding exactly and avoids the converter altogether. A regression
# here would silently break on macOS libhdf5 2.x and crash on Linux CI.
# ---------------------------------------------------------------------------
test_that("CLI obs transfer renames column-order <-> colnames without crash", {
  cli_bin <- scConvert:::sc_find_cli()
  skip_if(is.null(cli_bin), "CLI binary not present")
  skip_if_not_installed("hdf5r")

  # A scalar column-order does not exercise the vlen array write path that
  # triggered the crash; use >= 2 columns so it is genuinely a vlen array.
  obs_cols <- c("n_counts", "group", "sample_id")
  n_cells  <- 10L
  n_genes  <- 5L

  cell_names <- paste0("Cell_", seq_len(n_cells))
  gene_names <- paste0("Gene_", seq_len(n_genes))

  tmp_h5ad     <- tempfile(fileext = ".h5ad")
  tmp_h5seurat <- tempfile(fileext = ".h5seurat")
  tmp_h5ad_rt  <- tempfile(fileext = ".h5ad")
  on.exit({
    unlink(tmp_h5ad)
    unlink(tmp_h5seurat)
    unlink(tmp_h5ad_rt)
  }, add = TRUE)

  set.seed(42)
  x_dense <- matrix(rpois(n_cells * n_genes, 3),
                    nrow = n_cells, ncol = n_genes)
  x_csr <- as(x_dense, "RsparseMatrix")

  str_vlen <- hdf5r::H5T_STRING$new(size = Inf)
  scalar_sp <- hdf5r::H5S$new(type = "scalar")
  array_sp <- hdf5r::H5S$new(type = "simple",
                             dims    = length(obs_cols),
                             maxdims = length(obs_cols))

  h5 <- hdf5r::H5File$new(tmp_h5ad, mode = "w")

  x_grp <- h5$create_group("X")
  x_grp$create_dataset("data",    robj = as.double(x_csr@x),
                       dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  x_grp$create_dataset("indices", robj = as.integer(x_csr@j),
                       dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_dataset("indptr",  robj = as.integer(x_csr@p),
                       dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_attr("shape", robj = c(n_cells, n_genes),
                    dtype = hdf5r::h5types$H5T_NATIVE_INT)
  x_grp$create_attr("encoding-type", robj = "csr_matrix",
                    dtype = str_vlen, space = scalar_sp)

  obs_grp <- h5$create_group("obs")
  obs_grp$create_dataset("_index",    robj = cell_names)
  obs_grp$create_dataset("n_counts",  robj = as.double(rpois(n_cells, 500)),
                         dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE)
  obs_grp$create_dataset("group",
                         robj = rep(c("A", "B"), length.out = n_cells))
  obs_grp$create_dataset("sample_id",
                         robj = as.integer(rep(1:2, length.out = n_cells)),
                         dtype = hdf5r::h5types$H5T_NATIVE_INT)
  obs_grp$create_attr("encoding-type", robj = "dataframe",
                      dtype = str_vlen, space = scalar_sp)
  obs_grp$create_attr("_index", robj = "_index",
                      dtype = str_vlen, space = scalar_sp)
  # The vlen string array attribute that triggered the original segfault.
  obs_grp$create_attr("column-order", robj = obs_cols,
                      dtype = str_vlen, space = array_sp)

  var_grp <- h5$create_group("var")
  var_grp$create_dataset("_index", robj = gene_names)
  var_grp$create_attr("encoding-type", robj = "dataframe",
                      dtype = str_vlen, space = scalar_sp)
  var_grp$create_attr("_index", robj = "_index",
                      dtype = str_vlen, space = scalar_sp)

  h5$close_all()

  # Forward: h5ad -> h5seurat. Previously crashed at step [2/6] on Linux.
  ok_fwd <- scConvert_cli(tmp_h5ad, tmp_h5seurat,
                          overwrite = TRUE, verbose = FALSE)
  expect_true(ok_fwd)
  expect_true(file.exists(tmp_h5seurat))

  h5s  <- hdf5r::H5File$new(tmp_h5seurat, mode = "r")
  meta <- h5s[["meta.data"]]
  expect_true(meta$attr_exists("colnames"),
              info = "meta.data must carry the renamed colnames attribute")
  expect_false(meta$attr_exists("column-order"),
               info = "h5seurat meta.data must not retain h5ad column-order")
  meta_cols <- meta$attr_open("colnames")$read()
  expect_setequal(meta_cols, obs_cols)
  meta$close()
  h5s$close_all()

  # Reverse: h5seurat -> h5ad. Symmetrical rename path.
  ok_rev <- scConvert_cli(tmp_h5seurat, tmp_h5ad_rt,
                          overwrite = TRUE, verbose = FALSE)
  expect_true(ok_rev)
  expect_true(file.exists(tmp_h5ad_rt))

  h5rt   <- hdf5r::H5File$new(tmp_h5ad_rt, mode = "r")
  obs_rt <- h5rt[["obs"]]
  expect_true(obs_rt$attr_exists("column-order"),
              info = "roundtripped h5ad obs must carry column-order")
  expect_false(obs_rt$attr_exists("colnames"),
               info = "h5ad obs must not retain h5seurat colnames")
  rt_cols <- obs_rt$attr_open("column-order")$read()
  expect_setequal(rt_cols, obs_cols)

  # AnnData >= 0.8 requires encoding-type on every obs column dataset.
  for (col in obs_cols) {
    dset <- obs_rt[[col]]
    expect_true(dset$attr_exists("encoding-type"),
                info = paste("obs column", col, "missing encoding-type"))
    dset$close()
  }
  obs_rt$close()
  h5rt$close_all()
})

# ---------------------------------------------------------------------------
# Regression 5: h5ad writers must write obsm as cells x dimensions
#
# A previous C writer regression wrote Seurat reductions to /obsm with the
# HDF5 shape dimensions x cells. AnnData requires obsm arrays to have the
# observation axis first, so such files failed anndata.read_h5ad().
# ---------------------------------------------------------------------------
test_that("scConvert writes obsm reductions with n_obs as the HDF5 first axis", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if(Sys.which("python3") == "", "python3 not available")

  h5py_check <- system2("python3", c("-c", shQuote("import h5py")),
                        stdout = TRUE, stderr = TRUE)
  skip_if(!is.null(attr(h5py_check, "status")), "Python h5py not available")

  set.seed(11)
  n_cells <- 12L
  n_genes <- 20L
  counts <- matrix(rpois(n_cells * n_genes, 3),
                   nrow = n_genes, ncol = n_cells,
                   dimnames = list(paste0("g", seq_len(n_genes)),
                                   paste0("c", seq_len(n_cells))))
  counts <- Matrix::Matrix(counts, sparse = TRUE)

  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  Seurat::VariableFeatures(obj) <- rownames(obj)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = 3, verbose = FALSE)

  expect_h5py_obsm_shape <- function(path) {
    code <- sprintf(
      paste(
        "import h5py",
        "f = h5py.File(%s, 'r')",
        "print(','.join(map(str, f['obsm/X_pca'].shape)))",
        "print(str(f['obs/_index'].shape[0]))",
        sep = "\n"
      ),
      shQuote(path)
    )
    out <- system2("python3", c("-c", shQuote(code)),
                   stdout = TRUE, stderr = TRUE)
    expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
    expect_equal(out[1], paste(n_cells, 3L, sep = ","))
    expect_equal(out[2], as.character(n_cells))
  }

  # Single-assay Seurat objects exercise the fast C writer.
  tmp_c <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp_c), add = TRUE)
  scConvert(obj, dest = tmp_c, overwrite = TRUE, verbose = FALSE)
  expect_h5py_obsm_shape(tmp_c)

  # A second assay forces the R/hdf5r direct writer fallback.
  obj[["extra"]] <- Seurat::CreateAssayObject(counts = counts)
  tmp_r <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp_r), add = TRUE)
  scConvert(obj, dest = tmp_r, overwrite = TRUE, verbose = FALSE)
  expect_h5py_obsm_shape(tmp_r)
})

# ---------------------------------------------------------------------------
# Regression 6: HDF5 2.x strict datatype conversion (CSET mismatch)
#
# hdf5r writes vlen string attributes with H5T_CSET_ASCII. The C reader
# previously hardcoded an H5T_CSET_UTF8 vlen memtype for H5Aread calls,
# which HDF5 1.14 silently converted but HDF5 2.x refuses with
# "no appropriate function for conversion path". The downstream symptom
# is that h5ad categorical columns (seurat_annotations, seurat_clusters,
# RNA_snn_res.0.5) get dropped during CLI h5ad -> h5seurat, breaking
# DimPlot(group.by = "seurat_annotations") with "undefined columns".
# Fix: read each string attribute using the file's own datatype so the
# converter is bypassed. See commit 71f9207.
# ---------------------------------------------------------------------------
test_that("CLI preserves h5ad categoricals across HDF5 2.x cset boundary", {
  cli_bin <- scConvert:::sc_find_cli()
  skip_if(is.null(cli_bin), "CLI binary not present")
  skip_if_not_installed("hdf5r")

  h5ad_in <- system.file("testdata", "pbmc_small.h5ad", package = "scConvert")
  skip_if(!nzchar(h5ad_in), "Shipped pbmc_small.h5ad not found")

  # Source dataset has four h5ad categorical columns. Each maps to a
  # h5seurat factor group {levels, values} after CLI conversion.
  expected_factor_cols <- c(
    "RNA_snn_res.0.5", "orig.ident",
    "seurat_annotations", "seurat_clusters"
  )

  tmp_h5s <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(tmp_h5s), add = TRUE)

  ok <- scConvert_cli(h5ad_in, tmp_h5s, overwrite = TRUE, verbose = FALSE)
  expect_true(ok)

  h5s <- hdf5r::H5File$new(tmp_h5s, mode = "r")
  on.exit(try(h5s$close_all(), silent = TRUE), add = TRUE)

  meta <- h5s[["meta.data"]]
  meta_names <- names(meta)

  for (col in expected_factor_cols) {
    expect_true(
      col %in% meta_names,
      info = paste("meta.data missing factor column", col,
                   "-- HDF5 2.x cset bug regression")
    )
    grp <- meta[[col]]
    expect_true(inherits(grp, "H5Group"),
                info = paste(col, "should be a factor group, not a dataset"))
    expect_true("levels" %in% names(grp),
                info = paste(col, "missing levels dataset"))
    expect_true("values" %in% names(grp),
                info = paste(col, "missing values dataset"))

    levels_chr <- grp[["levels"]]$read()
    expect_gt(length(levels_chr), 0L)
    expect_true(all(nzchar(levels_chr)),
                info = paste(col, "levels contains empty strings"))

    values_int <- grp[["values"]]$read()
    expect_true(all(values_int >= 1L & values_int <= length(levels_chr)),
                info = paste(col, "values out of range -- 0-based codes leaked"))
    grp$close()
  }
})
