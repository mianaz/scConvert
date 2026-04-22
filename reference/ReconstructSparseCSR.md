# Reconstruct sparse CSR/CSC matrix from array components

Reconstruct sparse CSR/CSC matrix from array components

## Usage

``` r
ReconstructSparseCSR(data, indices, indptr, shape, transpose = TRUE)
```

## Arguments

- data:

  Numeric vector of non-zero values

- indices:

  Integer vector of column (CSR) or row (CSC) indices (0-based)

- indptr:

  Integer vector of row (CSR) or column (CSC) pointers

- shape:

  Integer vector of length 2: c(n_rows, n_cols)

- transpose:

  If TRUE, transpose result (h5ad stores cells x genes, Seurat needs
  genes x cells)

## Value

A dgCMatrix
