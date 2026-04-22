# Read an AnnData DataFrame column from zarr store

Handles categorical, string-array, and numeric encodings.

## Usage

``` r
.zarr_read_anndata_column(store_path, col_path)
```

## Arguments

- store_path:

  Path to zarr store

- col_path:

  Full relative path to column (e.g., "obs/cell_type")

## Value

Vector (factor, character, numeric, or logical)
