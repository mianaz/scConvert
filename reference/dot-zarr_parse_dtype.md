# Parse zarr v2 dtype string to R type info

Parse zarr v2 dtype string to R type info

## Usage

``` r
.zarr_parse_dtype(dtype)
```

## Arguments

- dtype:

  Zarr dtype string (e.g., "\<f8", "\<i4", "\|O")

## Value

List with r_type, size, endian, is_object
