# Read zarr node attributes

Reads .zattrs (v2) or zarr.json attributes field (v3).

## Usage

``` r
.zarr_read_attrs(store_path, rel_path = "")
```

## Arguments

- store_path:

  Path to zarr store

- rel_path:

  Relative path within store (empty string for root)

## Value

Named list of attributes
