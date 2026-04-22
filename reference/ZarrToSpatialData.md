# Wrap a standard anndata zarr store as a SpatialData store

Takes an existing anndata zarr store and wraps it as a SpatialData store
by creating the SpatialData directory layout and placing the anndata
store as a table.

## Usage

``` r
ZarrToSpatialData(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input anndata .zarr directory

- dest:

  Path for output SpatialData .zarr directory

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
