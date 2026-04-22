# Write point coordinates as a SpatialData points directory

Writes a zarr DataFrame containing point annotations (e.g., transcript
positions) following the SpatialData points specification.

## Usage

``` r
.write_spatialdata_points(points, path, verbose = TRUE)
```

## Arguments

- points:

  Data frame or matrix with point data

- path:

  Output path for the points zarr group

- verbose:

  Show progress
