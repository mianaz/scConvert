# Read an AnnData matrix (dense or sparse) from zarr store

Dispatches based on encoding-type attribute.

## Usage

``` r
.zarr_read_anndata_matrix(store_path, rel_path, transpose = TRUE)
```

## Arguments

- store_path:

  Path to zarr store

- rel_path:

  Relative path to matrix node

- transpose:

  If TRUE (default), transpose cells x genes to genes x cells

## Value

A matrix or dgCMatrix
