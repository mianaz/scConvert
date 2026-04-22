# Deconstruct dgCMatrix to AnnData CSR components

AnnData stores expression matrices as cells x genes in CSR format. This
function converts a genes x cells dgCMatrix to CSR components for the
cells x genes layout.

## Usage

``` r
DeconstructSparseCSR(mat, coerce = TRUE)
```

## Arguments

- mat:

  A dgCMatrix (genes x cells)

- coerce:

  If TRUE (default), coerce mat to dgCMatrix first

## Value

A list with `data`, `indices` (0-based), `indptr`, and `shape` (cells,
genes)
