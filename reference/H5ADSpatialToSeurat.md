# Convert spatial coordinates from h5ad to Seurat format

Convert spatial coordinates from h5ad to Seurat format

## Usage

``` r
H5ADSpatialToSeurat(
  h5ad_file,
  seurat_obj = NULL,
  assay_name = "Spatial",
  verbose = TRUE
)
```

## Arguments

- h5ad_file:

  H5AD file handle or path

- seurat_obj:

  Seurat object to add spatial data to

- assay_name:

  Name of the assay to associate spatial data with

- verbose:

  Print progress messages

## Value

Modified Seurat object with spatial data
