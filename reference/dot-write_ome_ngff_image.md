# Write an image as an OME-NGFF zarr store

Writes a 3D array (H x W x C) as an OME-NGFF zarr image with multiscales
metadata. This creates a single resolution level.

## Usage

``` r
.write_ome_ngff_image(image, path, name = "image", verbose = TRUE)
```

## Arguments

- image:

  Numeric 3D array (height, width, channels)

- path:

  Output path for the image zarr group

- name:

  Name for the multiscales metadata

- verbose:

  Show progress
