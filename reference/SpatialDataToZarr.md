# Convert a SpatialData zarr store to a standard anndata zarr store

Extracts the anndata table from a SpatialData store and writes it as a
standalone anndata zarr store (without the SpatialData wrapper).

## Usage

``` r
SpatialDataToZarr(
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

  Path for output .zarr directory

- table:

  Name of the table within the SpatialData store (default: "table")

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
