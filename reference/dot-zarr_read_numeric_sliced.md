# Chunk-selective read of a 1D or 2D zarr v2 numeric array

Chunk-selective read of a 1D or 2D zarr v2 numeric array

## Usage

``` r
.zarr_read_numeric_sliced(
  store,
  rel_path,
  meta,
  dtype_info,
  shape,
  chunks,
  compressor,
  slice
)
```

## Arguments

- slice:

  Named list with `d1`, optionally `d2`: each either `NULL` (full dim)
  or an integer vector of 1-based indices into that dim.
