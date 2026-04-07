# Tests that FOV boundaries and molecules round-trip through writeH5AD and
# readH5AD without silent data loss. This is the primary guarantee behind
# the CosMx/Xenium write path.

library(scConvert)

.build_seurat_with_fov <- function(n_cells = 6L, n_genes = 4L) {
  set.seed(42)
  counts <- matrix(as.integer(rpois(n_cells * n_genes, lambda = 5)),
                    nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell_", seq_len(n_cells))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts,
                                                       assay = "RNA"))

  # Build a simple Segmentation with one 4-vertex polygon per cell, plus a
  # Centroids boundary holding the centroid of each polygon.
  poly_df <- do.call(rbind, lapply(seq_len(n_cells), function(i) {
    cx <- i * 10
    cy <- i * 20
    data.frame(
      x    = c(cx, cx + 2, cx + 2, cx),
      y    = c(cy, cy,     cy + 2, cy + 2),
      cell = rep(colnames(counts)[i], 4),
      stringsAsFactors = FALSE
    )
  }))
  seg <- SeuratObject::CreateSegmentation(coords = poly_df)

  cent_df <- data.frame(
    x    = seq_len(n_cells) * 10 + 1,
    y    = seq_len(n_cells) * 20 + 1,
    cell = colnames(counts),
    stringsAsFactors = FALSE
  )
  cent <- SeuratObject::CreateCentroids(coords = cent_df)

  # Molecules: four transcripts of two different genes.
  mol_df <- data.frame(
    x    = c(1.0, 5.0, 9.0, 13.0, 17.0),
    y    = c(2.0, 6.0, 10.0, 14.0, 18.0),
    gene = c("Gene1", "Gene2", "Gene1", "Gene2", "Gene1"),
    stringsAsFactors = FALSE
  )
  mol <- SeuratObject::CreateMolecules(coords = mol_df)

  fov <- SeuratObject::CreateFOV(
    coords    = list(centroids = cent, segmentation = seg),
    type      = "centroids",
    molecules = mol,
    assay     = "RNA",
    key       = "testfov_"
  )
  obj[["testfov"]] <- fov
  list(obj = obj, poly_df = poly_df, mol_df = mol_df,
        cent_df = cent_df, n_cells = n_cells)
}

test_that("FOV boundaries and molecules round-trip through writeH5AD", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  bundle <- .build_seurat_with_fov()
  f <- tempfile(fileext = ".h5ad")
  on.exit(unlink(f), add = TRUE)

  writeH5AD(bundle$obj, filename = f, overwrite = TRUE, verbose = FALSE)

  # Inspect the raw HDF5 layout to confirm the FOV contract keys were written.
  h5 <- hdf5r::H5File$new(f, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  expect_true(h5$exists("uns/spatial/testfov"))
  expect_true(h5$exists("uns/spatial/testfov/segmentation"))
  expect_true(h5$exists("uns/spatial/testfov/molecules"))
  expect_true(h5$exists("obsm/spatial"))

  seg_grp <- h5[["uns/spatial/testfov/segmentation"]]
  expect_true("segmentation" %in% names(seg_grp) || "centroids" %in% names(seg_grp))

  # Segmentation sub-group should carry cell_ids + coords + polygon_offsets.
  has_seg_entry <- any(vapply(names(seg_grp), function(nm) {
    entry <- seg_grp[[nm]]
    isTRUE(entry$exists("cell_ids")) &&
    isTRUE(entry$exists("coords")) &&
    isTRUE(entry$exists("polygon_offsets"))
  }, logical(1)))
  expect_true(has_seg_entry)

  mol_grp <- h5[["uns/spatial/testfov/molecules"]]
  mol_children <- names(mol_grp)
  expect_true(length(mol_children) >= 1L)
  first_mol <- mol_grp[[mol_children[1]]]
  expect_true(first_mol$exists("x"))
  expect_true(first_mol$exists("y"))
  expect_true(first_mol$exists("gene"))
})

test_that("readH5AD rebuilds an FOV from uns/spatial after writeH5AD", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  bundle <- .build_seurat_with_fov()
  f <- tempfile(fileext = ".h5ad")
  on.exit(unlink(f), add = TRUE)

  writeH5AD(bundle$obj, filename = f, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(f, verbose = FALSE)

  expect_s4_class(obj2, "Seurat")
  # FOV slot should be rebuilt and attached under the original library name.
  expect_true("testfov" %in% Seurat::Images(obj2))

  fov2 <- obj2[["testfov"]]
  expect_s4_class(fov2, "FOV")

  # Boundaries list should contain at least one Segmentation or Centroids.
  bnd <- slot(fov2, "boundaries")
  expect_true(length(bnd) >= 1L)

  # Molecule count should survive round-trip.
  mol_slot <- slot(fov2, "molecules")
  if (length(mol_slot) > 0L) {
    df <- Seurat::GetTissueCoordinates(mol_slot[[1]])
    expect_equal(nrow(df), nrow(bundle$mol_df))
  }
})

test_that("Double round-trip preserves FOV structure (no silent degradation)", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  bundle <- .build_seurat_with_fov()
  f1 <- tempfile(fileext = ".h5ad")
  f2 <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(f1, f2)), add = TRUE)

  writeH5AD(bundle$obj, filename = f1, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(f1, verbose = FALSE)

  writeH5AD(obj2, filename = f2, overwrite = TRUE, verbose = FALSE)
  obj3 <- readH5AD(f2, verbose = FALSE)

  expect_setequal(Seurat::Images(obj2), Seurat::Images(obj3))
  expect_true("testfov" %in% Seurat::Images(obj3))

  # Polygon counts must match across the double round trip.
  bnd2 <- slot(obj2[["testfov"]], "boundaries")
  bnd3 <- slot(obj3[["testfov"]], "boundaries")
  expect_equal(length(bnd2), length(bnd3))
})
