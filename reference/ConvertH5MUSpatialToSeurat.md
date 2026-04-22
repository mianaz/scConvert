# Convert h5mu spatial data to Seurat format for a specific modality

Convert h5mu spatial data to Seurat format for a specific modality

## Usage

``` r
ConvertH5MUSpatialToSeurat(
  h5mu_file,
  seurat_obj,
  modality,
  assay_name,
  verbose = FALSE
)
```

## Arguments

- h5mu_file:

  H5File connection to h5mu file

- seurat_obj:

  Seurat object

- modality:

  Modality name in h5mu file

- assay_name:

  Corresponding assay name in Seurat

- verbose:

  Logical

## Value

Modified Seurat object
