# Tests for native Stereo-seq .gef reading (pure R, no reticulate).
#
# Builds minimal synthetic GEF files directly with hdf5r so the tests ship
# without downloading real data. The schema follows the STOmics specification:
#   /geneExp/bin{N}/gene        {geneID, geneName, offset, count}
#   /geneExp/bin{N}/expression  {x, y, count}

library(scConvert)

skip_unless_hdf5 <- function() {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("Seurat")
}

# Write a minimal square-bin .gef with `genes` and `spots` at the requested
# bin size. Returns a list with the file path and the dense ground-truth
# expression matrix (gene x spot).
.build_synthetic_squarebin_gef <- function(bin_size = 50) {
  genes <- data.frame(
    geneID   = c("ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004",
                 "ENSG000005"),
    geneName = c("Gpr1", "Tpt1", "Actb", "Gapdh", "Malat1"),
    stringsAsFactors = FALSE
  )
  # Three spots at (x, y) coordinates; gene counts vary per spot.
  # Expression table below lists them in row-major (gene, spot) order:
  #   gene1: spots 1,2       counts 3, 2
  #   gene2: spots 1,3       counts 1, 4
  #   gene3: spot 2          count 5
  #   gene4: spots 1,2,3     counts 1, 1, 1
  #   gene5: (none)
  exp_x <- as.integer(c(10, 20,   10, 30,   20,   10, 20, 30))
  exp_y <- as.integer(c( 5, 15,    5, 25,   15,    5, 15, 25))
  exp_c <- as.integer(c( 3,  2,    1,  4,    5,    1,  1,  1))
  genes$offset <- as.integer(c(0, 2, 4, 5, 8))
  genes$count  <- as.integer(c(2, 2, 1, 3, 0))

  tmp <- tempfile(fileext = ".gef")
  h5  <- hdf5r::H5File$new(tmp, mode = "w")

  # Root-level BGI attributes.
  h5$create_attr("bin_type", robj = "bin",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  h5$create_attr("sn", robj = "TEST001",
                  dtype = hdf5r::H5T_STRING$new(size = Inf),
                  space = hdf5r::H5S$new(type = "scalar"))
  h5$create_attr("resolution", robj = 500L,
                  dtype = hdf5r::h5types$H5T_NATIVE_UINT32,
                  space = hdf5r::H5S$new(type = "scalar"))

  bin_path <- sprintf("geneExp/bin%d", bin_size)
  h5$create_group("geneExp")
  h5[["geneExp"]]$create_group(sprintf("bin%d", bin_size))
  bin_grp <- h5[[bin_path]]

  # Gene compound dataset.
  gene_ct <- hdf5r::H5T_COMPOUND$new(
    labels = c("geneID", "geneName", "offset", "count"),
    dtypes = list(
      hdf5r::H5T_STRING$new(size = 32L),
      hdf5r::H5T_STRING$new(size = 32L),
      hdf5r::h5types$H5T_NATIVE_UINT32,
      hdf5r::h5types$H5T_NATIVE_UINT32
    )
  )
  bin_grp$create_dataset(
    name = "gene",
    robj = genes,
    dtype = gene_ct
  )

  # Expression compound dataset.
  expr_df <- data.frame(
    x     = exp_x,
    y     = exp_y,
    count = as.integer(exp_c),
    stringsAsFactors = FALSE
  )
  expr_ct <- hdf5r::H5T_COMPOUND$new(
    labels = c("x", "y", "count"),
    dtypes = list(
      hdf5r::h5types$H5T_NATIVE_INT32,
      hdf5r::h5types$H5T_NATIVE_INT32,
      hdf5r::h5types$H5T_NATIVE_UINT32
    )
  )
  bin_grp$create_dataset(
    name = "expression",
    robj = expr_df,
    dtype = expr_ct
  )
  h5$close_all()

  # Ground truth: dense gene x spot matrix in the same spot order as the
  # reader will produce (unique (x,y) pairs in first-seen order).
  spot_key <- paste(exp_x, exp_y, sep = "_")
  uniq_keys <- unique(spot_key)
  n_spot <- length(uniq_keys)
  dense <- matrix(0, nrow = nrow(genes), ncol = n_spot,
                  dimnames = list(genes$geneName, uniq_keys))
  gene_idx <- rep.int(seq_len(nrow(genes)), times = genes$count)
  spot_idx <- match(spot_key, uniq_keys)
  for (k in seq_along(gene_idx)) {
    dense[gene_idx[k], spot_idx[k]] <- exp_c[k]
  }

  list(path = tmp, dense = dense, n_genes = nrow(genes), n_spots = n_spot)
}

test_that("LoadStereoSeqGef reads minimal square-bin .gef end-to-end", {
  skip_unless_hdf5()

  bundle <- .build_synthetic_squarebin_gef(bin_size = 50)
  on.exit(unlink(bundle$path), add = TRUE)

  obj <- LoadStereoSeqGef(bundle$path, bin_size = 50, verbose = FALSE)

  expect_s4_class(obj, "Seurat")
  expect_equal(nrow(obj), bundle$n_genes)
  expect_equal(ncol(obj), bundle$n_spots)
  expect_equal(obj@misc$spatial_technology, "StereoSeq")
  expect_equal(obj@misc$stereo_seq$bin_size, 50)
  expect_true("spatial_x" %in% colnames(obj@meta.data))
  expect_true("spatial_y" %in% colnames(obj@meta.data))

  # Values must round-trip bit-exactly (reordering spots by content is OK —
  # compare sorted cell-sums and total counts instead of per-index equality
  # since the spot factor ordering depends on the bit-packed key).
  counts_mat <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
  expect_equal(sum(counts_mat), sum(bundle$dense))
  expect_equal(sort(unname(colSums(counts_mat))),
               sort(unname(colSums(bundle$dense))))
  expect_equal(sort(unname(rowSums(counts_mat))),
               sort(unname(rowSums(bundle$dense))))
  expect_setequal(rownames(counts_mat), rownames(bundle$dense))
})

test_that("LoadStereoSeqGef rejects inconsistent gene$count vs expression", {
  skip_unless_hdf5()

  # Build a file where the gene table claims more records than the
  # expression table contains. Should abort with a helpful message.
  bundle <- .build_synthetic_squarebin_gef(bin_size = 50)
  on.exit(unlink(bundle$path), add = TRUE)

  # Patch the gene$count field to deliberately mismatch.
  h5 <- hdf5r::H5File$new(bundle$path, mode = "r+")
  h5[["geneExp/bin50"]]$link_delete("gene")
  genes_bad <- data.frame(
    geneID   = c("ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004",
                 "ENSG000005"),
    geneName = c("Gpr1", "Tpt1", "Actb", "Gapdh", "Malat1"),
    offset   = as.integer(c(0, 2, 4, 5, 8)),
    count    = as.integer(c(2, 2, 1, 3, 99)),
    stringsAsFactors = FALSE
  )
  gene_ct <- hdf5r::H5T_COMPOUND$new(
    labels = c("geneID", "geneName", "offset", "count"),
    dtypes = list(
      hdf5r::H5T_STRING$new(size = 32L),
      hdf5r::H5T_STRING$new(size = 32L),
      hdf5r::h5types$H5T_NATIVE_UINT32,
      hdf5r::h5types$H5T_NATIVE_UINT32
    )
  )
  h5[["geneExp/bin50"]]$create_dataset("gene", robj = genes_bad, dtype = gene_ct)
  h5$close_all()

  expect_error(
    LoadStereoSeqGef(bundle$path, bin_size = 50, verbose = FALSE),
    regexp = "Corrupt Stereo-seq GEF"
  )
})

test_that("LoadStereoSeqGef errors for missing bin size with available list", {
  skip_unless_hdf5()

  bundle <- .build_synthetic_squarebin_gef(bin_size = 50)
  on.exit(unlink(bundle$path), add = TRUE)

  expect_error(
    LoadStereoSeqGef(bundle$path, bin_size = 200, verbose = FALSE),
    regexp = "bin200|Available bins"
  )
})
