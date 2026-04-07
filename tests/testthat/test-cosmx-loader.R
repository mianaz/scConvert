# Tests for LoadCosMx() — the thin Seurat::LoadNanostring wrapper.

library(scConvert)

# Build a minimal synthetic CosMx bundle. NanoString's AtoMx export produces
# five files; we write the subset that Seurat::ReadNanostring requires.
.build_synthetic_cosmx_bundle <- function() {
  bundle_dir <- tempfile("cosmx_")
  dir.create(bundle_dir)

  # Expression matrix: first column is fov, second cell, remaining = genes.
  # Negative probes are in columns starting with 'SystemControl' or 'Negative'.
  n_cells   <- 6L
  gene_names <- c("GeneA", "GeneB", "GeneC", "GeneD")
  cell_ids   <- seq_len(n_cells)
  fov_ids    <- rep(1L, n_cells)
  counts <- matrix(as.integer(runif(n_cells * length(gene_names), 0, 20)),
                   nrow = n_cells, ncol = length(gene_names))
  expr_df <- data.frame(
    fov    = fov_ids,
    cell_ID = cell_ids,
    counts,
    stringsAsFactors = FALSE
  )
  colnames(expr_df)[3:(2 + length(gene_names))] <- gene_names
  write.csv(expr_df,
            file.path(bundle_dir, "synth_exprMat_file.csv"),
            row.names = FALSE, quote = FALSE)

  # Metadata file: fov, cell_ID, x/y centroid columns, Area, Total counts.
  meta_df <- data.frame(
    fov        = fov_ids,
    cell_ID    = cell_ids,
    CenterX_global_px = seq_len(n_cells) * 10,
    CenterY_global_px = seq_len(n_cells) * 15,
    Area       = rep(100, n_cells),
    Mean.MembraneStain = rep(1.0, n_cells),
    Mean.DAPI  = rep(1.0, n_cells),
    Mean.CD298_B2M = rep(1.0, n_cells),
    Mean.PanCK = rep(1.0, n_cells),
    Mean.CD45  = rep(1.0, n_cells),
    Max.MembraneStain = rep(1.0, n_cells),
    Max.DAPI   = rep(1.0, n_cells),
    Max.CD298_B2M = rep(1.0, n_cells),
    Max.PanCK  = rep(1.0, n_cells),
    Max.CD45   = rep(1.0, n_cells),
    stringsAsFactors = FALSE
  )
  write.csv(meta_df,
            file.path(bundle_dir, "synth_metadata_file.csv"),
            row.names = FALSE, quote = FALSE)

  # FOV positions file.
  fov_df <- data.frame(
    fov       = 1L,
    x_global_px = 0,
    y_global_px = 0,
    stringsAsFactors = FALSE
  )
  write.csv(fov_df,
            file.path(bundle_dir, "synth_fov_positions_file.csv"),
            row.names = FALSE, quote = FALSE)

  # Transcript (molecule) file.
  tx_df <- data.frame(
    fov        = rep(1L, 8),
    cell_ID    = rep(cell_ids[1:4], 2),
    x_local_px = runif(8) * 50,
    y_local_px = runif(8) * 50,
    x_global_px = runif(8) * 50,
    y_global_px = runif(8) * 50,
    z          = rep(0L, 8),
    target     = rep(gene_names[1:4], 2),
    CellComp   = rep("Nuclear", 8),
    stringsAsFactors = FALSE
  )
  write.csv(tx_df,
            file.path(bundle_dir, "synth_tx_file.csv"),
            row.names = FALSE, quote = FALSE)

  bundle_dir
}

test_that(".is_cosmx_dir detects canonical bundles and rejects unrelated dirs", {
  bundle <- .build_synthetic_cosmx_bundle()
  on.exit(unlink(bundle, recursive = TRUE), add = TRUE)

  expect_true(scConvert:::.is_cosmx_dir(bundle))

  empty_dir <- tempfile("empty_")
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)
  expect_false(scConvert:::.is_cosmx_dir(empty_dir))

  expect_false(scConvert:::.is_cosmx_dir("/nonexistent/path/xyz"))
})

test_that("FileType() recognises a CosMx bundle directory as 'cosmx'", {
  bundle <- .build_synthetic_cosmx_bundle()
  on.exit(unlink(bundle, recursive = TRUE), add = TRUE)

  expect_equal(scConvert:::FileType(bundle), "cosmx")
})

test_that("LoadCosMx errors on a non-CosMx directory with a helpful message", {
  empty_dir <- tempfile("empty_")
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)

  expect_error(
    LoadCosMx(empty_dir, verbose = FALSE),
    regexp = "does not look like a CosMx bundle"
  )
})
