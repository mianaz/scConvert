# Read shapes from a SpatialData shapes directory

SpatialData shapes are stored as zarr DataFrames with x, y columns and
optionally a radius column. This follows the anndata on-disk DataFrame
encoding (vlen-utf8 string index + numeric columns).

## Usage

``` r
.read_spatialdata_shapes(path, verbose = TRUE)
```

## Arguments

- path:

  Path to the shapes zarr group (e.g., `shapes/spots`)

- verbose:

  Show progress

## Value

A list with `coords` (numeric matrix with x, y columns) and `radius`
(numeric scalar or NULL). Returns NULL if the shapes cannot be read.
