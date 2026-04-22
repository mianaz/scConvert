# Convert a SpatialData zarr store to h5Seurat format

Reads a SpatialData store and writes it as an h5Seurat file.

## Usage

``` r
SpatialDataToH5Seurat(
  source,
  dest,
  table = "table",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- source:

  Path to input SpatialData .zarr directory

- dest:

  Path for output .h5Seurat file

- table:

  Name of the table within the SpatialData store (default: "table")

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
