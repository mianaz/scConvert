library(scConvert)

test_that("SOMA saver lambda accepts filename= (dispatcher contract)", {
  skip_if_not_installed("Seurat")
  saver <- scConvert:::GetSaver(ext = "soma")
  expect_false(is.null(saver))
  expect_true("filename" %in% names(formals(saver)) ||
                "..." %in% names(formals(saver)))

  # Observational check: dispatcher passes filename=, not dest=.
  # Verify the saver function signature accepts filename positionally or by name.
  fn_args <- names(formals(saver))
  expect_true("filename" %in% fn_args,
              info = paste("saver args:", paste(fn_args, collapse = ",")))
})

test_that("SpatialData saver lambda accepts filename= (dispatcher contract)", {
  saver <- scConvert:::GetSaver(ext = "spatialdata.zarr")
  expect_false(is.null(saver))
  fn_args <- names(formals(saver))
  expect_true("filename" %in% fn_args,
              info = paste("saver args:", paste(fn_args, collapse = ",")))
})

test_that("scConvert(Seurat, x.soma) dispatches without 'dest is missing'", {
  skip_if_not_installed("Seurat")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".soma")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # We only care that the dispatch layer does not error with
  # 'argument "dest" is missing, with no default'. If tiledbsoma is
  # not installed, writeSOMA throws its own require-namespace error,
  # which is fine -- that proves dispatch succeeded.
  err <- tryCatch(
    scConvert(obj, tmp, verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("\"dest\" is missing", err, fixed = TRUE),
               info = paste("error was:", err))
})

test_that("scConvert(Seurat, x.spatialdata.zarr) dispatches without 'dest is missing'", {
  skip_if_not_installed("Seurat")

  counts <- matrix(1:12, nrow = 3)
  rownames(counts) <- paste0("Gene", 1:3)
  colnames(counts) <- paste0("Cell", 1:4)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  tmp <- tempfile(fileext = ".spatialdata.zarr")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  err <- tryCatch(
    scConvert(obj, tmp, verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("\"dest\" is missing", err, fixed = TRUE),
               info = paste("error was:", err))
})
