# Read a zarr v2 numeric array from disk

Read a zarr v2 numeric array from disk

## Usage

``` r
.zarr_read_numeric(store_path, rel_path, slice = NULL)
```

## Arguments

- store_path:

  Path to zarr store

- rel_path:

  Relative path to array within store

- slice:

  Optional named list of per-dimension index vectors:
  `list(d1 = 1L:1000L)` (1D), `list(d1 = ..., d2 = ...)` (2D). `NULL` on
  a dim means "all". When supplied, only the chunks that intersect the
  requested indices are fetched.

## Value

Numeric vector or matrix
