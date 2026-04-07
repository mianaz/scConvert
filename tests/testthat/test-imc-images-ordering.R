# Regression test for the IMC double-roundtrip 14/15 -> 11/13 degradation.
#
# Root cause: SeuratSpatialToH5AD() previously processed only Images()[1],
# silently dropping additional IMC panels. Additionally, HDF5 group iteration
# order is unstable, so cross-roundtrip image selection reshuffled between
# rounds. Fix: deterministic sort + loop over all images.

library(scConvert)

.create_multi_image_seurat <- function(image_names,
                                        n_cells = 8L, n_genes = 4L) {
  set.seed(42)
  counts <- matrix(as.integer(rpois(n_cells * n_genes, lambda = 5)),
                    nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell_", seq_len(n_cells))
  obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts,
                                                       assay = "RNA"))

  # For each requested image name build a Centroids-only FOV.
  for (img_name in image_names) {
    df <- data.frame(
      x    = seq_len(n_cells) * 10,
      y    = seq_len(n_cells) * 20,
      cell = colnames(counts),
      stringsAsFactors = FALSE
    )
    cent <- SeuratObject::CreateCentroids(coords = df)
    fov <- SeuratObject::CreateFOV(
      coords = cent, type = "centroids",
      assay = "RNA", key = paste0(img_name, "_")
    )
    obj[[img_name]] <- fov
  }
  obj
}

test_that("Multi-image IMC roundtrip preserves every image across two cycles", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  # Deliberately non-alphabetical insertion order so sorted selection differs
  # from insertion order on the write side.
  img_names <- c("imc_panel_b", "imc_panel_a", "imc_panel_c")
  obj <- .create_multi_image_seurat(img_names)
  expect_equal(length(Seurat::Images(obj)), 3L)

  f1 <- tempfile(fileext = ".h5ad")
  f2 <- tempfile(fileext = ".h5ad")
  on.exit(unlink(c(f1, f2)), add = TRUE)

  writeH5AD(obj, filename = f1, overwrite = TRUE, verbose = FALSE)
  obj2 <- readH5AD(f1, verbose = FALSE)

  # All three images must survive the first round-trip.
  expect_setequal(Seurat::Images(obj2), img_names)

  writeH5AD(obj2, filename = f2, overwrite = TRUE, verbose = FALSE)
  obj3 <- readH5AD(f2, verbose = FALSE)

  # Second round-trip must NOT degrade image count or identity.
  expect_setequal(Seurat::Images(obj3), Seurat::Images(obj2))
  expect_equal(length(Seurat::Images(obj3)), length(Seurat::Images(obj2)))
})

test_that("writeH5AD writes uns/spatial subgroups in sorted order", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")

  img_names <- c("zimage", "aimage", "mimage")
  obj <- .create_multi_image_seurat(img_names)

  f <- tempfile(fileext = ".h5ad")
  on.exit(unlink(f), add = TRUE)
  writeH5AD(obj, filename = f, overwrite = TRUE, verbose = FALSE)

  h5 <- hdf5r::H5File$new(f, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  expect_true(h5$exists("uns/spatial"))
  sp_names <- names(h5[["uns/spatial"]])
  expect_setequal(sp_names, img_names)
})
