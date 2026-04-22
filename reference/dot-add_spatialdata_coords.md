# Add spatial coordinates from SpatialData to Seurat object

Add spatial coordinates from SpatialData to Seurat object

## Usage

``` r
.add_spatialdata_coords(
  obj,
  coords,
  radius = 1,
  region_name = "spots",
  cell_names = NULL,
  instance_key = NULL,
  verbose = TRUE
)
```

## Arguments

- obj:

  Seurat object

- coords:

  Coordinate matrix (n x 2) with x, y columns

- radius:

  Spot radius

- region_name:

  Name of the spatial region

- cell_names:

  Cell names in the Seurat object

- instance_key:

  Column in obs linking to shapes index

- verbose:

  Show progress

## Value

Modified Seurat object
