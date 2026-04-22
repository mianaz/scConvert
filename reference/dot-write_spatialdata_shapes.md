# Write spatial coordinates as a SpatialData shapes directory

Writes a zarr DataFrame with x, y, and radius columns following the
SpatialData shapes specification.

## Usage

``` r
.write_spatialdata_shapes(coords, radius = 1, path, verbose = TRUE)
```

## Arguments

- coords:

  Numeric matrix with x, y columns (rows = cells/spots)

- radius:

  Numeric scalar for spot radius

- path:

  Output path for the shapes zarr group

- verbose:

  Show progress
