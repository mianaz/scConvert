# Read sparse matrix from HDF5

Internal function to read sparse matrix data from h5Seurat files,
handling both V4 and V5 formats

## Usage

``` r
ReadSparseMatrix(h5_group, verbose = FALSE)
```

## Arguments

- h5_group:

  The HDF5 group containing sparse matrix data

- verbose:

  Show progress updates

## Value

A matrix object, or NULL if reading fails
