# Streaming (out-of-memory) h5ad <-> h5Seurat conversion must give exactly
# the result of the in-memory conversion.

library(scConvert)

.q <- function(expr) suppressWarnings(suppressMessages(expr))

.expect_same_object <- function(a, b, label) {
  expect_identical(SeuratObject::Assays(a), SeuratObject::Assays(b), info = label)
  for (as in SeuratObject::Assays(a)) {
    la <- SeuratObject::Layers(a[[as]])
    expect_setequal(la, SeuratObject::Layers(b[[as]]))
    for (l in la) {
      ma <- SeuratObject::LayerData(a, assay = as, layer = l)
      mb <- SeuratObject::LayerData(b, assay = as, layer = l)
      expect_identical(dimnames(ma), dimnames(mb), info = paste(label, as, l))
      expect_equal(as.matrix(ma), as.matrix(mb), info = paste(label, as, l))
    }
    expect_identical(SeuratObject::VariableFeatures(a[[as]]),
                     SeuratObject::VariableFeatures(b[[as]]), info = label)
  }
  expect_equal(a[[]], b[[]], info = paste(label, "meta.data"))
  expect_setequal(SeuratObject::Reductions(a), SeuratObject::Reductions(b))
  for (r in SeuratObject::Reductions(a)) {
    expect_equal(SeuratObject::Embeddings(a, r), SeuratObject::Embeddings(b, r), info = label)
  }
  expect_setequal(SeuratObject::Graphs(a), SeuratObject::Graphs(b))
  for (g in SeuratObject::Graphs(a)) {
    expect_equal(as.matrix(a[[g]]), as.matrix(b[[g]]), info = paste(label, g))
  }
}

.both_ways <- function(src, ext, reader) {
  out <- lapply(c(FALSE, TRUE), function(st) {
    d <- tempfile(fileext = paste0(".", ext))
    withr::defer(unlink(d), envir = parent.frame(3))
    .q(scConvert(src, dest = d, overwrite = TRUE, verbose = FALSE, stream = st))
    .q(reader(d, verbose = FALSE))
  })
  out
}

test_that("streamed h5ad -> h5Seurat equals the in-memory conversion", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  withr::local_options(scConvert.stream_chunk_bytes = 2048)   # many small blocks
  o <- .q(SeuratObject::UpdateSeuratObject(SeuratObject::pbmc_small))
  o <- .q(Seurat::FindNeighbors(o, dims = 1:5, k.param = 5, verbose = FALSE))
  src <- tempfile(fileext = ".h5ad")
  on.exit(unlink(src), add = TRUE)
  .q(writeH5AD(o, src, overwrite = TRUE, verbose = FALSE))
  for (f in c(src, system.file("extdata", "spatial_demo.h5ad", package = "scConvert"))) {
    if (!nzchar(f)) next
    r <- .both_ways(f, "h5seurat", readH5Seurat)
    .expect_same_object(r[[1]], r[[2]], basename(f))
  }
})

test_that("streaming handles raw/X with its own features, dense and CSC layers", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  withr::local_options(scConvert.stream_chunk_bytes = 2048)
  set.seed(7)
  o <- .q(SeuratObject::UpdateSeuratObject(SeuratObject::pbmc_small))
  src <- tempfile(fileext = ".h5ad")
  on.exit(unlink(src), add = TRUE)
  .q(writeH5AD(o, src, overwrite = TRUE, verbose = FALSE))
  h <- hdf5r::H5File$new(src, mode = "r+")
  genes <- h[["var/_index"]][]
  n <- length(h[["obs/_index"]][])
  raw_genes <- sample(c(genes, paste0("EXTRA_", 1:10)))
  rawm <- methods::as(Matrix::rsparsematrix(n, length(raw_genes), 0.2), "RsparseMatrix")
  rawm@x <- round(abs(rawm@x) * 10) + 1
  sc <- function(o, nm, v) o$create_attr(nm, v, space = hdf5r::H5S$new("scalar"))
  rg <- h$create_group("raw")
  xg <- rg$create_group("X")
  xg[["data"]] <- rawm@x
  xg[["indices"]] <- rawm@j
  xg$create_dataset("indptr", robj = rawm@p, dtype = hdf5r::h5types$H5T_STD_I64LE)
  sc(xg, "encoding-type", "csr_matrix"); sc(xg, "encoding-version", "0.1.0")
  xg$create_attr("shape", as.integer(c(n, length(raw_genes))))
  vg <- rg$create_group("var")
  vg[["_index"]] <- raw_genes
  sc(vg, "_index", "_index"); sc(vg, "encoding-type", "dataframe"); sc(vg, "encoding-version", "0.2.0")
  vg$create_attr("column-order", character(0), space = hdf5r::H5S$new(dims = 0))
  if (h$exists("layers/counts")) h[["layers"]]$link_delete("counts")
  if (h$exists("uns/scConvert")) h[["uns"]]$link_delete("scConvert")
  h[["layers"]]$create_dataset("dense", robj = matrix(runif(length(genes) * n), nrow = length(genes)))
  sc(h[["layers/dense"]], "encoding-type", "array"); sc(h[["layers/dense"]], "encoding-version", "0.2.0")
  cm <- Matrix::rsparsematrix(n, length(genes), 0.1)
  cg <- h[["layers"]]$create_group("csc")
  cg[["data"]] <- cm@x; cg[["indices"]] <- cm@i; cg[["indptr"]] <- cm@p
  sc(cg, "encoding-type", "csc_matrix"); sc(cg, "encoding-version", "0.1.0")
  cg$create_attr("shape", as.integer(c(n, length(genes))))
  h$close_all()

  r <- .both_ways(src, "h5seurat", readH5Seurat)
  .expect_same_object(r[[1]], r[[2]], "raw/X + dense + csc")
  counts <- as.matrix(SeuratObject::LayerData(r[[2]], layer = "counts"))
  expect_equal(unname(counts), unname(as.matrix(Matrix::t(rawm))[match(genes, raw_genes), ]))
})

test_that("streamed h5Seurat -> h5ad equals the in-memory conversion", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  withr::local_options(scConvert.stream_chunk_bytes = 2048)
  o <- .q(SeuratObject::UpdateSeuratObject(SeuratObject::pbmc_small))
  o[["RNA"]] <- methods::as(o[["RNA"]], "Assay5")
  src <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(src), add = TRUE)
  writeH5Seurat(o, src, overwrite = TRUE, verbose = FALSE)
  r <- .both_ways(src, "h5ad", readH5AD)
  .expect_same_object(r[[1]], r[[2]], "h5seurat -> h5ad")
})

test_that("stream = 'auto' follows the size threshold", {
  f <- tempfile()
  writeBin(raw(100), f)
  on.exit(unlink(f), add = TRUE)
  ns <- asNamespace("scConvert")
  withr::local_options(scConvert.stream_threshold = 50)
  expect_true(ns$.use_streaming("auto", f))
  withr::local_options(scConvert.stream_threshold = 500)
  expect_false(ns$.use_streaming("auto", f))
  expect_true(ns$.use_streaming(TRUE, f))
  expect_error(ns$.use_streaming("yes", f))
})

test_that("scConvertBatch converts several files and reports failures", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  o <- .q(SeuratObject::UpdateSeuratObject(SeuratObject::pbmc_small))
  dir <- tempfile("batch")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  srcs <- file.path(dir, paste0("s", 1:3, ".h5ad"))
  for (f in srcs) .q(writeH5AD(o, f, overwrite = TRUE, verbose = FALSE))
  bad <- file.path(dir, "missing.h5ad")
  res <- .q(scConvertBatch(c(srcs, bad), dest = "h5seurat", workers = 2, stream = TRUE))
  expect_identical(res$ok, c(TRUE, TRUE, TRUE, FALSE))
  expect_true(all(file.exists(sub("\\.h5ad$", ".h5seurat", srcs))))
  expect_match(res$error[4], "not found")
  x <- .q(readH5Seurat(res$dest[1], verbose = FALSE))
  expect_equal(dim(x), dim(o))
})

test_that("every HDF5 pair streams to the same result as in memory", {
  skip_on_cran()
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("Seurat")
  withr::local_options(scConvert.stream_chunk_bytes = 4096)
  o <- .q(SeuratObject::UpdateSeuratObject(SeuratObject::pbmc_small))
  o[["RNA"]] <- methods::as(o[["RNA"]], "Assay5")
  rd <- list(h5ad = readH5AD, h5seurat = readH5Seurat, h5mu = readH5MU, loom = readLoom)
  wr <- list(h5ad = writeH5AD, h5seurat = writeH5Seurat, h5mu = writeH5MU, loom = writeLoom)
  dir <- withr::local_tempdir()
  for (s in names(rd)) {
    src <- file.path(dir, paste0("src.", s))
    .q(wr[[s]](o, src, overwrite = TRUE, verbose = FALSE))
    for (d in setdiff(names(rd), s)) {
      label <- paste(s, "->", d)
      fallback <- character(0)
      out <- vapply(c(FALSE, TRUE), function(st) {
        f <- file.path(dir, paste0("out_", st, ".", d))
        withCallingHandlers(
          suppressWarnings(scConvert(src, dest = f, overwrite = TRUE, verbose = TRUE, stream = st)),
          message = function(m) {
            if (grepl("Streaming not possible", conditionMessage(m))) fallback <<- conditionMessage(m)
            invokeRestart("muffleMessage")
          })
        f
      }, character(1))
      expect_identical(fallback, character(0), info = label)
      a <- .q(rd[[d]](out[1], verbose = FALSE))
      b <- .q(rd[[d]](out[2], verbose = FALSE))
      .expect_same_object(a, b, label)
    }
  }
})
