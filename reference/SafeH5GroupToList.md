# Safely convert H5Group to list, handling 3D+ arrays

This function converts an HDF5 group to an R list, properly handling
datasets with 3 or more dimensions (which fail with standard as.list).
This is needed for Squidpy UMAP data that stores 3D arrays.

## Usage

``` r
SafeH5GroupToList(h5obj, recursive = TRUE)
```

## Arguments

- h5obj:

  An H5Group or H5D object

- recursive:

  Whether to recursively convert subgroups

## Value

A list representation of the HDF5 structure
