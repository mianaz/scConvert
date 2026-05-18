# Mirror a remote zarr store to a local directory

Fetches every object under the prefix in parallel-friendly serial GETs,
preserving relative paths. Uses an on-disk cache keyed by URL when
\`cache_dir\` is non-NULL: a subsequent call against the same URL
returns the cached directory without re-downloading.

## Usage

``` r
.zarr_fetch_remote(url, cache_dir = NULL, verbose = TRUE)
```

## Arguments

- url:

  Remote URL (s3://, gs://).

- cache_dir:

  Cache directory, or NULL for tempdir (no persistence).

- verbose:

  Print progress.

## Value

Path to local directory containing the mirrored store.
