# Add Visium spatial data to Seurat object

Add Visium spatial data to Seurat object

## Usage

``` r
AddVisiumSpatialData(
  seurat_obj,
  coords,
  visium_data,
  h5ad,
  assay_name = "Spatial",
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object

- coords:

  Spatial coordinates

- visium_data:

  Visium-specific metadata

- assay_name:

  Assay name

- verbose:

  Print messages

## Value

Modified Seurat object
