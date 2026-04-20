# Property-based roundtrip tests.
#
# Each seed produces a randomized Seurat object with varying layer, metadata,
# embedding, and graph content. The object is pushed through a format-specific
# roundtrip (Seurat -> file -> Seurat) and asserted equal to the original at
# every axis we claim to preserve. Failures here indicate loss of fidelity in
# the conversion layer, not in user code.

suppressPackageStartupMessages({
  library(scConvert)
  library(Seurat)
  library(Matrix)
})

# --- Generator -------------------------------------------------------------

gen_random_seurat <- function(seed,
                              n_cells = NULL,
                              n_genes = NULL,
                              with_data = TRUE,
                              with_reductions = TRUE,
                              with_graphs = TRUE,
                              with_varp = FALSE) {
  set.seed(seed)
  if (is.null(n_cells)) n_cells <- sample(20:80, 1)
  if (is.null(n_genes)) n_genes <- sample(30:120, 1)

  nnz <- max(1L, round(n_cells * n_genes * 0.15))
  counts <- sparseMatrix(
    i = sample.int(n_genes, nnz, replace = TRUE),
    j = sample.int(n_cells, nnz, replace = TRUE),
    x = as.double(rpois(nnz, 5) + 1L),
    dims = c(n_genes, n_cells),
    dimnames = list(paste0("G", seq_len(n_genes)),
                    paste0("C", seq_len(n_cells)))
  )
  obj <- CreateSeuratObject(counts = counts)
  if (with_data) obj <- NormalizeData(obj, verbose = FALSE)

  # Heterogeneous obs columns: numeric, integer, factor, logical, character
  obj$num_col  <- rnorm(n_cells)
  obj$int_col  <- sample.int(100, n_cells, replace = TRUE)
  obj$fct_col  <- factor(sample(c("A", "B", "C"), n_cells, replace = TRUE),
                         levels = c("A", "B", "C", "D_unused"))
  obj$bool_col <- sample(c(TRUE, FALSE), n_cells, replace = TRUE)
  obj$chr_col  <- paste0("tag_", sample.int(5, n_cells, replace = TRUE))

  # var metadata: one numeric, one logical (mirrors "highly_variable")
  feat <- GetAssay(obj, "RNA")
  feat_meta <- data.frame(
    mean_expr = rowMeans(as.matrix(counts)),
    is_hvg    = sample(c(TRUE, FALSE), n_genes, replace = TRUE),
    row.names = rownames(counts)
  )
  obj[["RNA"]] <- AddMetaData(obj[["RNA"]], feat_meta)

  if (with_reductions) {
    n_pc <- 5L
    emb <- matrix(rnorm(n_cells * n_pc), nrow = n_cells, ncol = n_pc,
                  dimnames = list(colnames(counts), paste0("PC_", seq_len(n_pc))))
    load <- matrix(rnorm(n_genes * n_pc), nrow = n_genes, ncol = n_pc,
                   dimnames = list(rownames(counts), paste0("PC_", seq_len(n_pc))))
    obj[["pca"]] <- CreateDimReducObject(embeddings = emb, loadings = load,
                                         key = "PC_", assay = "RNA")
  }

  if (with_graphs) {
    nnz_g <- round(n_cells * 3)
    g <- sparseMatrix(
      i = sample.int(n_cells, nnz_g, replace = TRUE),
      j = sample.int(n_cells, nnz_g, replace = TRUE),
      x = runif(nnz_g),
      dims = c(n_cells, n_cells),
      dimnames = list(colnames(counts), colnames(counts))
    )
    obj[["RNA_nn"]] <- as.Graph(g)
  }

  if (with_varp) {
    nnz_v <- round(n_genes * 2)
    vp <- sparseMatrix(
      i = sample.int(n_genes, nnz_v, replace = TRUE),
      j = sample.int(n_genes, nnz_v, replace = TRUE),
      x = runif(nnz_v),
      dims = c(n_genes, n_genes),
      dimnames = list(rownames(counts), rownames(counts))
    )
    Misc(obj, slot = "__varp__") <- list(gene_corr = vp)
  }

  obj
}

# --- Comparators -----------------------------------------------------------

expect_matrix_equal <- function(a, b, label) {
  a <- as(a, "CsparseMatrix"); b <- as(b, "CsparseMatrix")
  expect_equal(dim(a), dim(b), info = paste(label, "dim"))
  expect_equal(rownames(a), rownames(b), info = paste(label, "rownames"))
  expect_equal(colnames(a), colnames(b), info = paste(label, "colnames"))
  expect_equal(a@i, b@i, info = paste(label, "row indices"))
  expect_equal(a@p, b@p, info = paste(label, "col pointers"))
  expect_equal(a@x, b@x, info = paste(label, "values"), tolerance = 1e-10)
}

expect_obs_equal <- function(a, b) {
  # Ignore orig.ident/nCount_RNA/nFeature_RNA (recomputed or renamed across formats);
  # focus on the user-added columns we generated.
  cols <- intersect(colnames(a), colnames(b))
  cols <- setdiff(cols, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                          "n_counts", "n_genes"))
  for (col in cols) {
    av <- a[[col]]; bv <- b[[col]]
    if (is.factor(av) || is.factor(bv)) {
      expect_equal(as.character(av), as.character(bv),
                   info = paste("obs col", col))
    } else if (is.logical(av) || is.logical(bv)) {
      expect_equal(as.logical(av), as.logical(bv),
                   info = paste("obs col", col))
    } else {
      expect_equal(av, bv, info = paste("obs col", col))
    }
  }
}

# --- h5ad roundtrip --------------------------------------------------------

test_that("h5ad roundtrip preserves counts, data, and obs metadata", {
  for (seed in 1:5) {
    obj <- gen_random_seurat(seed, with_graphs = FALSE, with_reductions = FALSE)
    tmp <- tempfile(fileext = ".h5ad")
    on.exit(unlink(tmp), add = TRUE)

    writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
    back <- readH5AD(tmp, verbose = FALSE)

    expect_matrix_equal(
      GetAssayData(obj,  assay = "RNA", layer = "counts"),
      GetAssayData(back, assay = "RNA", layer = "counts"),
      label = paste0("seed=", seed, " counts")
    )
    expect_matrix_equal(
      GetAssayData(obj,  assay = "RNA", layer = "data"),
      GetAssayData(back, assay = "RNA", layer = "data"),
      label = paste0("seed=", seed, " data")
    )
    expect_obs_equal(obj[[]], back[[]])
  }
})

test_that("h5ad roundtrip preserves obsm reductions", {
  for (seed in 11:13) {
    obj <- gen_random_seurat(seed, with_reductions = TRUE, with_graphs = FALSE)
    tmp <- tempfile(fileext = ".h5ad")
    on.exit(unlink(tmp), add = TRUE)

    writeH5AD(obj, tmp, overwrite = TRUE, verbose = FALSE)
    back <- readH5AD(tmp, verbose = FALSE)

    expect_true("pca" %in% Reductions(back),
                info = paste0("seed=", seed))
    expect_equal(
      unname(Embeddings(obj,  reduction = "pca")),
      unname(Embeddings(back, reduction = "pca")),
      tolerance = 1e-10,
      info = paste0("seed=", seed, " embeddings")
    )
  }
})

# --- h5seurat roundtrip ----------------------------------------------------

test_that("h5seurat roundtrip preserves counts and obs metadata", {
  for (seed in 21:23) {
    obj <- gen_random_seurat(seed, with_data = FALSE,
                             with_graphs = FALSE, with_reductions = FALSE)
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)

    writeH5Seurat(obj, tmp, overwrite = TRUE, verbose = FALSE)
    back <- readH5Seurat(tmp, verbose = FALSE)

    expect_matrix_equal(
      GetAssayData(obj,  assay = "RNA", layer = "counts"),
      GetAssayData(back, assay = "RNA", layer = "counts"),
      label = paste0("seed=", seed, " counts")
    )
    expect_obs_equal(obj[[]], back[[]])
  }
})

# --- Zarr roundtrip --------------------------------------------------------

test_that("zarr roundtrip preserves counts, data, and obs metadata", {
  for (seed in 31:33) {
    obj <- gen_random_seurat(seed, with_graphs = FALSE, with_reductions = FALSE)
    tmp <- tempfile(fileext = ".zarr")
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

    writeZarr(obj, tmp, overwrite = TRUE, verbose = FALSE)
    back <- readZarr(tmp, verbose = FALSE)

    expect_matrix_equal(
      GetAssayData(obj,  assay = "RNA", layer = "counts"),
      GetAssayData(back, assay = "RNA", layer = "counts"),
      label = paste0("seed=", seed, " counts")
    )
    expect_matrix_equal(
      GetAssayData(obj,  assay = "RNA", layer = "data"),
      GetAssayData(back, assay = "RNA", layer = "data"),
      label = paste0("seed=", seed, " data")
    )
    expect_obs_equal(obj[[]], back[[]])
  }
})

# --- Multi-hop: h5ad -> h5seurat -> h5ad (byte-level on terminal file) -----

test_that("multi-hop h5ad -> h5seurat -> h5ad preserves matrix content", {
  for (seed in 41:43) {
    obj <- gen_random_seurat(seed, with_graphs = FALSE, with_reductions = FALSE)
    f1 <- tempfile(fileext = ".h5ad")
    f2 <- tempfile(fileext = ".h5seurat")
    f3 <- tempfile(fileext = ".h5ad")
    on.exit(unlink(c(f1, f2, f3)), add = TRUE)

    writeH5AD(obj, f1, overwrite = TRUE, verbose = FALSE)
    scConvert(f1, f2, overwrite = TRUE, verbose = FALSE)
    scConvert(f2, f3, overwrite = TRUE, verbose = FALSE)
    back <- readH5AD(f3, verbose = FALSE)

    expect_matrix_equal(
      GetAssayData(obj,  assay = "RNA", layer = "counts"),
      GetAssayData(back, assay = "RNA", layer = "counts"),
      label = paste0("seed=", seed, " multi-hop counts")
    )
  }
})
