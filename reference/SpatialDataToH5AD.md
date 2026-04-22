# Convert a SpatialData zarr store to h5ad format

Extracts the anndata table from a SpatialData store and writes it as an
h5ad file. Spatial coordinates and images are embedded in the h5ad using
standard scanpy/squidpy conventions (`obsm/spatial`, `uns/spatial`).

## Usage

``` r
SpatialDataToH5AD(
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

  Path for output .h5ad file

- table:

  Name of the table within the SpatialData store (default: "table")

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
