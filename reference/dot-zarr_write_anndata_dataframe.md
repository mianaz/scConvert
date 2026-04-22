# Write an AnnData DataFrame (obs or var) to zarr store

Write an AnnData DataFrame (obs or var) to zarr store

## Usage

``` r
.zarr_write_anndata_dataframe(
  store_path,
  group_path,
  df,
  index,
  compressor,
  verbose = TRUE
)
```

## Arguments

- store_path:

  Path to zarr store

- group_path:

  Relative path for the group (e.g., "obs", "var")

- df:

  Data frame to write

- index:

  Character vector of index values (cell/gene names)

- compressor:

  Compressor spec

- verbose:

  Show progress
