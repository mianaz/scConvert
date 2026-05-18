test_that(".is_remote_zarr_url detects cloud URL schemes", {
  expect_true(scConvert:::.is_remote_zarr_url("s3://bucket/key.zarr"))
  expect_true(scConvert:::.is_remote_zarr_url("S3://bucket/key.zarr"))
  expect_true(scConvert:::.is_remote_zarr_url("gs://bucket/key.zarr"))
  expect_true(scConvert:::.is_remote_zarr_url("https://example.org/x.zarr"))
  expect_true(scConvert:::.is_remote_zarr_url("http://example.org/x.zarr"))
  expect_false(scConvert:::.is_remote_zarr_url("/tmp/x.zarr"))
  expect_false(scConvert:::.is_remote_zarr_url("./x.zarr"))
  expect_false(scConvert:::.is_remote_zarr_url("x.zarr"))
  expect_false(scConvert:::.is_remote_zarr_url(NA_character_))
  expect_false(scConvert:::.is_remote_zarr_url(c("s3://a/b", "s3://c/d")))
})

test_that(".zarr_translate_url maps s3:// to virtual-hosted https", {
  i <- scConvert:::.zarr_translate_url("s3://my-bucket/path/to/store.zarr")
  expect_equal(i$kind, "s3")
  expect_equal(i$bucket, "my-bucket")
  expect_equal(i$prefix, "path/to/store.zarr")
  expect_equal(i$base, "https://my-bucket.s3.amazonaws.com/path/to/store.zarr")
  expect_equal(i$list_url, "https://my-bucket.s3.amazonaws.com")
})

test_that(".zarr_translate_url maps gs:// to storage.googleapis.com", {
  i <- scConvert:::.zarr_translate_url("gs://my-bucket/x.zarr")
  expect_equal(i$kind, "gs")
  expect_equal(i$bucket, "my-bucket")
  expect_equal(i$prefix, "x.zarr")
  expect_equal(i$base, "https://storage.googleapis.com/my-bucket/x.zarr")
})

test_that(".zarr_translate_url passes plain https through", {
  i <- scConvert:::.zarr_translate_url("https://example.org/x.zarr/")
  expect_equal(i$kind, "http")
  expect_equal(i$base, "https://example.org/x.zarr")
})

test_that(".zarr_translate_url strips trailing slashes", {
  i <- scConvert:::.zarr_translate_url("s3://bucket/key.zarr/")
  expect_equal(i$base, "https://bucket.s3.amazonaws.com/key.zarr")
})

test_that(".zarr_list_remote_keys requires httr and xml2 errors clearly", {
  # Cannot reliably uninstall packages mid-test; just verify the error path
  # produces a sensible message when called with a raw http kind.
  info <- scConvert:::.zarr_translate_url("https://example.org/x.zarr")
  expect_error(scConvert:::.zarr_list_remote_keys(info),
               "does not enumerate arbitrary HTTP")
})

test_that(".scconvert_cache_dir returns a usable path", {
  d <- scConvert:::.scconvert_cache_dir()
  expect_type(d, "character")
  expect_length(d, 1L)
  expect_true(nzchar(d))
})

test_that("readZarr against s3:// URL: opt-in live test", {
  url <- Sys.getenv("SCCONVERT_TEST_S3_ZARR", "")
  skip_if(!nzchar(url),
          "Set SCCONVERT_TEST_S3_ZARR=s3://bucket/key.zarr to run.")
  skip_if_not_installed("httr")
  skip_if_not_installed("xml2")
  skip_if_offline()
  obj <- readZarr(url, verbose = FALSE, cache = FALSE)
  expect_s4_class(obj, "Seurat")
})
