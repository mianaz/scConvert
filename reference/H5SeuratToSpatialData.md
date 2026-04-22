# Convert an h5Seurat file to SpatialData zarr format

Reads an h5Seurat file and writes it as a SpatialData zarr store.

## Usage

``` r
H5SeuratToSpatialData(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .h5Seurat file

- dest:

  Path for output SpatialData .zarr directory

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
