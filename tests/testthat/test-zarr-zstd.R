test_that(".zarr_resolve_compressor auto-picks zstd when zstdlite is available", {
  skip_if_not_installed("zstdlite")
  spec <- scConvert:::.zarr_resolve_compressor(NULL)
  expect_equal(spec$id, "zstd")
  expect_type(spec$level, "integer")
})

test_that(".zarr_resolve_compressor returns zlib when zstdlite is absent", {
  # We cannot reliably uninstall mid-test; mock requireNamespace. mockery::stub
  # reassigns the function by its local name, so bind it first -- stubbing the
  # `scConvert:::` expression directly leaves the namespace copy unmocked.
  skip_if_not_installed("mockery")
  fn <- scConvert:::.zarr_resolve_compressor
  mockery::stub(fn, "requireNamespace", function(...) FALSE)
  spec <- fn(NULL)
  expect_equal(spec$id, "zlib")
})

test_that(".zarr_resolve_compressor accepts strings and lists", {
  expect_equal(scConvert:::.zarr_resolve_compressor("zlib")$id, "zlib")
  expect_equal(scConvert:::.zarr_resolve_compressor("gzip")$id, "zlib")
  expect_equal(scConvert:::.zarr_resolve_compressor("none"), NULL)
  expect_equal(scConvert:::.zarr_resolve_compressor(FALSE), NULL)
  spec <- scConvert:::.zarr_resolve_compressor(list(id = "zstd", level = 9L))
  expect_equal(spec$id, "zstd")
  expect_equal(spec$level, 9L)
})

test_that(".zarr_resolve_compressor errors on unsupported codecs", {
  expect_error(scConvert:::.zarr_resolve_compressor("snappy"),
               "Unsupported zarr compressor")
})

test_that("zstd round-trip preserves bytes via .zarr_compress/.zarr_decompress", {
  skip_if_not_installed("zstdlite")
  raw <- as.raw(c(1:255, sample(0:255, 1000, TRUE)))
  spec <- list(id = "zstd", level = 3L)
  enc <- scConvert:::.zarr_compress(raw, spec)
  dec <- scConvert:::.zarr_decompress(enc, spec)
  expect_identical(dec, raw)
})

test_that("writeZarr defaults to zstd then roundtrips to readZarr", {
  skip_if_not_installed("zstdlite")
  obj <- suppressWarnings(SeuratObject::pbmc_small)
  dest <- tempfile(fileext = ".zarr")
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))

  writeZarr(obj, dest, verbose = FALSE)

  # Sparse X is stored as a group; the data chunk array carries the codec
  zarray <- jsonlite::fromJSON(file.path(dest, "X", "data", ".zarray"))
  expect_equal(zarray$compressor$id, "zstd")

  obj2 <- readZarr(dest, verbose = FALSE)
  expect_s4_class(obj2, "Seurat")
  expect_equal(ncol(obj), ncol(obj2))
  expect_equal(nrow(obj), nrow(obj2))
})

test_that("writeZarr respects compressor = 'zlib' override", {
  obj <- suppressWarnings(SeuratObject::pbmc_small)
  dest <- tempfile(fileext = ".zarr")
  on.exit(unlink(dest, recursive = TRUE, force = TRUE))

  writeZarr(obj, dest, compressor = "zlib", verbose = FALSE)

  zarray <- jsonlite::fromJSON(file.path(dest, "X", "data", ".zarray"))
  expect_equal(zarray$compressor$id, "zlib")
})
