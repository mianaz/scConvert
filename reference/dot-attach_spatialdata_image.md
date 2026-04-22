# Attach an image array to a Seurat object

Attach an image array to a Seurat object

## Usage

``` r
.attach_spatialdata_image(obj, img_data, img_name = "image", verbose = TRUE)
```

## Arguments

- obj:

  Seurat object (must already have spatial coordinates)

- img_data:

  3D array (H, W, C) normalized to \[0, 1\]

- img_name:

  Name for the image

- verbose:

  Show progress

## Value

Modified Seurat object
