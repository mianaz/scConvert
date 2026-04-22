# Write a zarr v2 numeric array

Write a zarr v2 numeric array

## Usage

``` r
.zarr_write_numeric(dir, data, dtype = NULL, compressor = NULL, attrs = NULL)
```

## Arguments

- dir:

  Directory for the array

- data:

  Numeric vector or matrix

- dtype:

  Zarr dtype string (default "\<f8" for float64)

- compressor:

  Compressor spec list (default gzip level 4)

- attrs:

  Optional attributes for .zattrs
