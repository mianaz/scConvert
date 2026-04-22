# Read a vlen-utf8 encoded string array from zarr v2 store

AnnData zarr stores encode string arrays using the vlen-utf8 filter.
Each string in the chunk is stored as a 4-byte little-endian int32
length prefix followed by that many bytes of UTF-8 text.

## Usage

``` r
.zarr_read_strings(store_path, rel_path)
```

## Arguments

- store_path:

  Path to zarr store

- rel_path:

  Relative path to string array within store

## Value

Character vector
