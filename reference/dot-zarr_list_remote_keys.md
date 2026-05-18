# List all object keys under a prefix in an S3/GCS bucket

Uses ListObjectsV2 (S3) or list_v1 (GCS) over anonymous HTTP. Paginates
via continuation tokens. Returns keys relative to the prefix.

## Usage

``` r
.zarr_list_remote_keys(info)
```

## Arguments

- info:

  Output of \`.zarr_translate_url()\`.

## Value

Character vector of keys relative to the prefix (no leading slash).
Empty character if listing failed.
