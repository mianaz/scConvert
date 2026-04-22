# Read observation embeddings from an H5AD file

Read the cell embeddings (obsm) from an H5AD file and return them as a
named list of matrices. Each entry corresponds to a dimensional
reduction (e.g., PCA, UMAP, t-SNE).

## Usage

``` r
readH5AD_obsm(file)
```

## Arguments

- file:

  Path to an H5AD file

## Value

A named list of matrices, where each matrix has cells as rows and
embedding dimensions as columns. Names are derived from the obsm keys
with the “X\_” prefix removed.
