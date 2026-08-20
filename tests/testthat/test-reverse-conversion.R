# Reverse-conversion (h5ad -> Seurat) integrity fixes:
#   * version/layout-aware counts resolution (layers/counts, raw/X, X,
#     /uns/scConvert stamp)
#   * non-integer counts warning
#   * post-read shape/orientation + name-identity verification
#   * duplicate barcode/gene handling
#   * reduction key collision handling + remapping
#   * categorical level order + ordered flag
#   * gene identity preservation (orig_var_index)
#   * read/write provenance records

library(scConvert)

# ---- helpers -----------------------------------------------------------------

scalar_str_attr <- function(obj, name, val) {
  invisible(obj$create_attr(name, robj = val,
                            dtype = hdf5r::H5T_STRING$new(size = Inf),
                            space = hdf5r::H5S$new(type = "scalar")))
}

write_csr <- function(parent, name, m) {
  # m is a dense (cells x genes) matrix; store as h5ad CSR
  csr <- as(Matrix::t(as(m, "CsparseMatrix")), "CsparseMatrix")
  g <- parent$create_group(name)
  g$create_dataset("data", robj = as.numeric(csr@x))
  g$create_dataset("indices", robj = as.integer(csr@i))
  g$create_dataset("indptr", robj = as.integer(csr@p))
  invisible(g$create_attr("shape", robj = as.integer(c(nrow(m), ncol(m)))))
  scalar_str_attr(g, "encoding-type", "csr_matrix")
  invisible(g)
}

# Build a synthetic h5ad; returns the ground-truth matrices/names
make_h5ad <- function(tmp, n_cells = 12, n_genes = 6, lognorm_x = TRUE,
                      layers_counts = TRUE, obsm = NULL, obs_cat = NULL,
                      stamp = NULL, gene_names = NULL, cell_names = NULL) {
  set.seed(7)
  counts <- matrix(rpois(n_cells * n_genes, 4), nrow = n_cells)
  xmat <- if (lognorm_x) log1p(counts / rowSums(counts) * 1e4) else counts
  if (is.null(cell_names)) cell_names <- paste0("Cell-", seq_len(n_cells))
  if (is.null(gene_names)) gene_names <- paste0("Gene-", seq_len(n_genes))
  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  write_csr(h5, "X", xmat)
  if (layers_counts) {
    write_csr(h5$create_group("layers"), "counts", counts)
  }
  obs <- h5$create_group("obs")
  obs$create_dataset("_index", robj = cell_names)
  scalar_str_attr(obs, "encoding-type", "dataframe")
  if (!is.null(obs_cat)) {
    cg <- obs$create_group("celltype")
    cg$create_dataset("codes", robj = obs_cat$codes)
    cg$create_dataset("categories", robj = obs_cat$categories)
    scalar_str_attr(cg, "encoding-type", "categorical")
    invisible(cg$create_attr("ordered", robj = isTRUE(obs_cat$ordered)))
    # anndata always records column-order; the C reader enumerates via it
    invisible(obs$create_attr("column-order", robj = "celltype"))
  }
  var <- h5$create_group("var")
  var$create_dataset("_index", robj = gene_names)
  scalar_str_attr(var, "encoding-type", "dataframe")
  var$create_dataset("gene_ids", robj = paste0("ENSG", seq_len(n_genes)))
  invisible(var$create_attr("column-order", robj = "gene_ids"))
  if (!is.null(obsm)) {
    og <- h5$create_group("obsm")
    for (nm in names(obsm)) og$create_dataset(nm, robj = obsm[[nm]])
  }
  if (!is.null(stamp)) {
    ug <- h5$create_group("uns")
    sg <- ug$create_group("scConvert")
    sg$create_dataset("counts_location", robj = stamp)
    sg$create_dataset("version", robj = "0.3.0")
  }
  h5$close_all()
  invisible(list(counts = counts, x = xmat, genes = gene_names,
                 cells = cell_names))
}

# ---- counts-layer resolution (item: counts-layer ambiguity) ------------------

test_that("layers/counts becomes the counts layer and X becomes data (both readers)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  info <- make_h5ad(tmp, lognorm_x = TRUE, layers_counts = TRUE)

  for (usec in c(FALSE, TRUE)) {
    obj <- suppressWarnings(readH5AD(tmp, verbose = FALSE, use.c = usec))
    cts <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
    dat <- as.matrix(Seurat::GetAssayData(obj, layer = "data"))
    # counts layer holds the integer counts from /layers/counts
    expect_true(all(cts == round(cts)), info = paste("use.c =", usec))
    expect_true(all(abs(t(cts) - info$counts) < 1e-8), info = paste("use.c =", usec))
    # data layer holds the (log-normalized) /X values
    expect_true(all(abs(t(dat) - info$x) < 1e-6), info = paste("use.c =", usec))
    # provenance records what happened
    prov <- obj@misc$scConvert_read
    expect_identical(prov$counts_source, "layers/counts")
    expect_identical(prov$x_mapped_to, "data")
    expect_true(isTRUE(prov$verified))
  }
})

test_that("layers/counts is honored even when `components` excludes layers", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  info <- make_h5ad(tmp, lognorm_x = TRUE, layers_counts = TRUE)

  for (usec in c(FALSE, TRUE)) {
    obj <- suppressWarnings(
      readH5AD(tmp, verbose = FALSE, use.c = usec, components = c("X", "obs"))
    )
    cts <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
    expect_true(all(cts == round(cts)), info = paste("use.c =", usec))
    expect_true(all(abs(t(cts) - info$counts) < 1e-8), info = paste("use.c =", usec))
  }
})

test_that("non-integer counts raise a classed scConvert_counts_warning", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  make_h5ad(tmp, lognorm_x = TRUE, layers_counts = FALSE)

  expect_warning(
    readH5AD(tmp, verbose = FALSE, use.c = FALSE),
    class = "scConvert_counts_warning"
  )

  # Integer-valued X without a dedicated counts slot stays quiet
  tmp2 <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp2), add = TRUE)
  make_h5ad(tmp2, lognorm_x = FALSE, layers_counts = FALSE)
  got <- character(0)
  withCallingHandlers(
    readH5AD(tmp2, verbose = FALSE, use.c = FALSE),
    warning = function(w) {
      got <<- c(got, class(w)[1])
      invokeRestart("muffleWarning")
    }
  )
  expect_false("scConvert_counts_warning" %in% got)
})

test_that("the /uns/scConvert counts_location stamp overrides layout heuristics", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  # File has layers/counts, but the stamp pins counts to X: the stamp wins.
  make_h5ad(tmp, lognorm_x = TRUE, layers_counts = TRUE, stamp = "X")

  obj <- suppressWarnings(readH5AD(tmp, verbose = FALSE, use.c = FALSE))
  expect_identical(obj@misc$scConvert_read$counts_source, "X")
})

# ---- post-read verification (items: transpose correctness, read-back) --------

test_that(".h5ad_verify_read catches orientation, identity, and duplicate faults", {
  skip_if_not_installed("Seurat")

  m <- matrix(rpois(50, 3), nrow = 5, ncol = 10,
              dimnames = list(paste0("g", 1:5), paste0("c", 1:10)))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(as(m, "CsparseMatrix")))

  # Correct expectation passes
  expect_true(scConvert:::.h5ad_verify_read(obj, paste0("c", 1:10),
                                            paste0("g", 1:5), "f.h5ad"))
  # Swapped dims = transpose/orientation fault
  err <- tryCatch(
    scConvert:::.h5ad_verify_read(obj, paste0("c", 1:5), paste0("g", 1:10), "f.h5ad"),
    error = function(e) e
  )
  expect_s3_class(err, "scConvert_data_error")
  expect_match(conditionMessage(err), "transpose/orientation")

  # Same dims, different cell order = identity fault
  err2 <- tryCatch(
    scConvert:::.h5ad_verify_read(obj, rev(paste0("c", 1:10)),
                                  paste0("g", 1:5), "f.h5ad"),
    error = function(e) e
  )
  expect_s3_class(err2, "scConvert_data_error")
  expect_match(conditionMessage(err2), "cell names")

  # Underscore mangling is tolerated (documented Seurat behavior)
  expect_true(scConvert:::.h5ad_verify_read(obj, paste0("c", 1:10),
                                            paste0("g", 1:5), "f.h5ad"))
  m2 <- m
  rownames(m2) <- paste0("g_", 1:5)
  obj2 <- suppressWarnings(Seurat::CreateSeuratObject(as(m2, "CsparseMatrix")))
  expect_true(scConvert:::.h5ad_verify_read(obj2, paste0("c", 1:10),
                                            paste0("g_", 1:5), "f.h5ad"))
})

# ---- duplicate names (items: cell order/identity, read-back verification) ----

test_that("duplicate feature names warn, dedupe, and preserve original identity", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  info <- make_h5ad(tmp, gene_names = c("G_1", "G_1", "G_2", "G_3", "G_4", "G_5"))

  for (usec in c(FALSE, TRUE)) {
    got <- character(0)
    obj <- withCallingHandlers(
      readH5AD(tmp, verbose = FALSE, use.c = usec),
      warning = function(w) {
        got <<- c(got, class(w)[1])
        invokeRestart("muffleWarning")
      }
    )
    expect_true("scConvert_names_warning" %in% got, info = paste("use.c =", usec))
    expect_equal(anyDuplicated(rownames(obj)), 0L)
    # counts relocation still lands despite underscore->dash mangling
    cts <- as.matrix(Seurat::GetAssayData(obj, layer = "counts"))
    expect_true(all(abs(t(cts) - info$counts) < 1e-8), info = paste("use.c =", usec))
    # original var index preserved as feature metadata
    fmeta <- obj[["RNA"]][[]]
    expect_true("orig_var_index" %in% colnames(fmeta))
    expect_identical(as.character(fmeta$orig_var_index), info$genes)
    # var columns still positionally aligned after mangling
    expect_identical(as.character(fmeta$gene_ids), paste0("ENSG", 1:6))
    expect_true(isTRUE(obj@misc$scConvert_read$duplicate_feature_names_renamed))
  }
})

test_that("duplicate cell barcodes warn, dedupe, and are recorded", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  cells <- c("AAA", "AAA", paste0("C", 3:12))
  make_h5ad(tmp, cell_names = cells)

  got <- character(0)
  obj <- withCallingHandlers(
    readH5AD(tmp, verbose = FALSE, use.c = FALSE),
    warning = function(w) {
      got <<- c(got, class(w)[1])
      invokeRestart("muffleWarning")
    }
  )
  expect_true("scConvert_names_warning" %in% got)
  expect_equal(anyDuplicated(colnames(obj)), 0L)
  expect_identical(colnames(obj), make.unique(cells))
  expect_true(isTRUE(obj@misc$scConvert_read$duplicate_cell_names_renamed))
})

# ---- reduction keys (item: reduction key collisions) -------------------------

test_that("colliding obsm keys are kept under distinct names with a warning", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  make_h5ad(tmp, obsm = list(X_pca = matrix(rnorm(36), 12, 3),
                             pca = matrix(rnorm(24), 12, 2)))

  for (usec in c(FALSE, TRUE)) {
    got <- character(0)
    obj <- withCallingHandlers(
      readH5AD(tmp, verbose = FALSE, use.c = usec),
      warning = function(w) {
        got <<- c(got, class(w)[1])
        invokeRestart("muffleWarning")
      }
    )
    expect_true("scConvert_reduction_warning" %in% got, info = paste("use.c =", usec))
    # Both embeddings survive (previously the second silently overwrote the first)
    expect_length(names(obj@reductions), 2L)
  }
})

test_that("`reductions` selects and renames obsm keys", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  make_h5ad(tmp, obsm = list(X_pca = matrix(rnorm(36), 12, 3),
                             X_umap = matrix(rnorm(24), 12, 2)))

  for (usec in c(FALSE, TRUE)) {
    obj <- suppressWarnings(
      readH5AD(tmp, verbose = FALSE, use.c = usec,
               reductions = c(scvi = "X_pca"))
    )
    expect_identical(names(obj@reductions), "scvi")
    expect_equal(ncol(Seurat::Embeddings(obj, "scvi")), 3L)
  }

  # Requesting a missing key warns with the reduction class
  expect_warning(
    suppressWarnings(
      readH5AD(tmp, verbose = FALSE, use.c = FALSE,
               reductions = c("X_pca", "X_missing")),
      classes = "updatedKeyWarning"
    ),
    class = "scConvert_reduction_warning"
  )
})

test_that(".h5ad_plan_reductions maps, excludes spatial, and resolves collisions", {
  plan <- scConvert:::.h5ad_plan_reductions(c("X_pca", "X_umap", "spatial"))
  expect_identical(unname(plan), c("X_pca", "X_umap"))
  expect_identical(names(plan), c("pca", "umap"))

  plan2 <- scConvert:::.h5ad_plan_reductions(c("X_pca", "X_umap", "spatial"),
                                             exclude_spatial = FALSE)
  expect_true("spatial" %in% names(plan2))

  expect_warning(
    plan3 <- scConvert:::.h5ad_plan_reductions(c("X_pca", "pca")),
    class = "scConvert_reduction_warning"
  )
  expect_equal(anyDuplicated(names(plan3)), 0L)
  expect_length(plan3, 2L)

  plan4 <- scConvert:::.h5ad_plan_reductions(c("X_pca", "X_umap"),
                                             reductions = c(lat = "X_pca"))
  expect_identical(plan4, stats::setNames("X_pca", "lat"))
})

# ---- categorical order + ordered flag (item: factor/category level order) ----

test_that("categorical level order and the ordered flag survive the read", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  # Deliberately non-alphabetical category order
  cats <- c("zeta", "alpha", "mid")
  make_h5ad(tmp, obs_cat = list(codes = as.integer(rep(c(0, 1, 2), 4)),
                                categories = cats, ordered = TRUE))

  for (usec in c(FALSE, TRUE)) {
    obj <- suppressWarnings(readH5AD(tmp, verbose = FALSE, use.c = usec))
    ct <- obj[["celltype"]][, 1]
    expect_s3_class(ct, "factor")
    expect_identical(levels(ct), cats, info = paste("use.c =", usec))
    expect_true(is.ordered(ct), info = paste("use.c =", usec))
  }

  # ordered = FALSE round-trips as an unordered factor
  tmp2 <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp2), add = TRUE)
  make_h5ad(tmp2, obs_cat = list(codes = as.integer(rep(c(0, 1, 2), 4)),
                                 categories = cats, ordered = FALSE))
  obj2 <- suppressWarnings(readH5AD(tmp2, verbose = FALSE, use.c = FALSE))
  ct2 <- obj2[["celltype"]][, 1]
  expect_identical(levels(ct2), cats)
  expect_false(is.ordered(ct2))
})

test_that("DecodeCategorical preserves category order and honors ordered=", {
  f <- scConvert:::DecodeCategorical(c(0L, 2L, 1L), c("z", "a", "m"))
  expect_identical(levels(f), c("z", "a", "m"))
  expect_false(is.ordered(f))
  f2 <- scConvert:::DecodeCategorical(c(0L, 2L, 1L), c("z", "a", "m"),
                                      ordered = TRUE)
  expect_true(is.ordered(f2))
  expect_identical(as.character(f2), c("z", "m", "a"))
})

# ---- writer provenance stamp (items: counts ambiguity, version pinning) ------

test_that("writeH5AD stamps /uns/scConvert and the round trip restores layers", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  set.seed(11)
  m <- matrix(rpois(60, 3), 10, 6,
              dimnames = list(paste0("G", 1:10), paste0("C", 1:6)))
  sobj <- suppressWarnings(
    Seurat::CreateSeuratObject(counts = as(m, "CsparseMatrix"))
  )
  sobj <- suppressWarnings(Seurat::NormalizeData(sobj, verbose = FALSE))

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  suppressWarnings(suppressMessages(writeH5AD(sobj, tmp, verbose = FALSE)))

  # Stamp present and correct: counts live in layers/counts
  h5 <- hdf5r::H5File$new(tmp, mode = "r")
  expect_true(h5$exists("uns"))
  expect_true(h5[["uns"]]$exists("scConvert"))
  expect_identical(
    as.character(h5[["uns"]][["scConvert"]][["counts_location"]]$read()),
    "layers/counts"
  )
  expect_identical(
    as.character(h5[["uns"]][["scConvert"]][["version"]]$read()),
    as.character(utils::packageVersion("scConvert"))
  )
  h5$close_all()

  # Round trip: integer counts back in counts, log-normalized values in data
  back <- suppressWarnings(readH5AD(tmp, verbose = FALSE))
  cts <- as.matrix(Seurat::GetAssayData(back, layer = "counts"))
  dat <- as.matrix(Seurat::GetAssayData(back, layer = "data"))
  expect_true(all(cts == round(cts)))
  expect_equal(cts, as.matrix(Seurat::GetAssayData(sobj, layer = "counts")),
               ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(dat, as.matrix(Seurat::GetAssayData(sobj, layer = "data")),
               ignore_attr = TRUE, tolerance = 1e-6)
  expect_identical(back@misc$scConvert_read$counts_source, "layers/counts")
})
