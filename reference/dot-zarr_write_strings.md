# Write a zarr v2 vlen-utf8 string array

Write a zarr v2 vlen-utf8 string array

## Usage

``` r
.zarr_write_strings(dir, strings, compressor = NULL, attrs = NULL)
```

## Arguments

- dir:

  Directory for the array

- strings:

  Character vector

- compressor:

  Compressor spec list (default gzip level 4)

- attrs:

  Optional attributes for .zattrs
