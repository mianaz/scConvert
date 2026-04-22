# Read points from a SpatialData points directory

Points are stored as zarr DataFrames. They typically represent
transcript locations with x, y, and gene columns. Because there can be
millions of points, these are stored in the Seurat object's `misc` slot
rather than as coordinates.

## Usage

``` r
.read_spatialdata_points(path, verbose = TRUE)
```

## Arguments

- path:

  Path to the points zarr group

- verbose:

  Show progress

## Value

A data.frame with columns read from the zarr DataFrame, or NULL.
