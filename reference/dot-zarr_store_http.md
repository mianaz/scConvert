# HTTP-backed zarr store with per-chunk lazy fetch

Lazy: each `read_bytes()` issues a fresh HTTP GET (cached on disk if
`cache_dir` is non-NULL). `list_dirs()` uses S3 ListObjectsV2 with the
slash delimiter. Once a store is opened, the full object listing is
fetched up front and memoised; subsequent
[`exists()`](https://rdrr.io/r/base/exists.html) / `list_dirs()` checks
hit the memoised manifest.

## Usage

``` r
.zarr_store_http(url, cache_dir = NULL)
```
