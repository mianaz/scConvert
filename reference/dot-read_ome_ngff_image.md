# Read an OME-NGFF image from a zarr directory

Reads the highest resolution level from an OME-NGFF zarr image. The
image is stored as a chunked zarr array with multiscales metadata
describing axes and resolution levels. Only the first dataset (highest
resolution) is read.

## Usage

``` r
.read_ome_ngff_image(path, verbose = TRUE)
```

## Arguments

- path:

  Path to the OME-NGFF zarr image group

- verbose:

  Show progress

## Value

A numeric array with dimensions (height, width, channels) normalized to
\\0, 1\\, suitable for a Seurat image slot. Returns NULL on failure.
