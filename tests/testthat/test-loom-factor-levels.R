# Test that scConvert's Loom writer preserves factor-column level ordering
# via the /col_attrs/__factor_levels__/{name} extension group, and that the
# reader restores them on load. Closes part of the 8/10 Loom fidelity gap.

library(scConvert)

test_that("Loom roundtrip preserves factor levels with non-alphabetical order", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  set.seed(42)
  n_cells <- 10L
  n_genes <- 5L
  counts <- matrix(as.integer(rpois(n_cells * n_genes, 3)),
                    nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell_", seq_len(n_cells))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts,
                                                       assay = "RNA"))

  # A factor with a deliberately non-alphabetical level ordering.
  ct_levels <- c("T_cell", "B_cell", "macrophage")
  obj$cell_type <- factor(
    sample(ct_levels, n_cells, replace = TRUE),
    levels = ct_levels
  )
  # A second factor, ordered (should round-trip as an ordered factor).
  sev_levels <- c("low", "medium", "high")
  obj$severity <- factor(
    sample(sev_levels, n_cells, replace = TRUE),
    levels = sev_levels,
    ordered = TRUE
  )

  original_ct  <- obj$cell_type
  original_sev <- obj$severity

  f <- tempfile(fileext = ".loom")
  on.exit(unlink(f), add = TRUE)
  writeLoom(obj, filename = f, overwrite = TRUE, verbose = FALSE)

  # scConvert_extensions group should hold the factor levels (outside
  # col_attrs so the Loom spec is not violated).
  h5 <- hdf5r::H5File$new(f, mode = "r")
  expect_true(h5$exists("scConvert_extensions/col_factor_levels"))
  expect_true(h5$exists("scConvert_extensions/col_factor_levels/cell_type"))
  expect_true(h5$exists("scConvert_extensions/col_factor_levels/severity"))
  h5$close_all()

  obj2 <- readLoom(f, verbose = FALSE)
  # Access via @meta.data directly: Seurat's [[<-/[[ accessors strip the
  # 'ordered' class attribute on unwrap, which would mask the round-trip.
  md2 <- obj2@meta.data

  expect_true(is.factor(md2$cell_type))
  expect_equal(levels(md2$cell_type), ct_levels)
  # Character content must round-trip too (Loom stores the raw values as
  # character; we only restore factor structure).
  expect_equal(as.character(md2$cell_type), as.character(original_ct))

  expect_true(is.factor(md2$severity))
  expect_true(is.ordered(md2$severity))
  expect_equal(levels(md2$severity), sev_levels)
  expect_equal(as.character(md2$severity), as.character(original_sev))
})
