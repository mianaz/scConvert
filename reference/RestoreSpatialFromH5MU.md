# Restore spatial data from h5mu file to Seurat object

Restore spatial data from h5mu file to Seurat object

## Usage

``` r
RestoreSpatialFromH5MU(
  h5mu_file,
  seurat_obj,
  modalities,
  assay.names,
  verbose = FALSE
)
```

## Arguments

- h5mu_file:

  H5File connection to h5mu file

- seurat_obj:

  Seurat object to add spatial data to

- modalities:

  Character vector of modalities that were loaded

- assay.names:

  Named vector mapping modality to assay names

- verbose:

  Logical; print messages

## Value

Modified Seurat object with spatial data
