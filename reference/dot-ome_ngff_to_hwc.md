# Convert raw OME-NGFF array to height x width x channels format

Convert raw OME-NGFF array to height x width x channels format

## Usage

``` r
.ome_ngff_to_hwc(arr, shape, axes_names = NULL)
```

## Arguments

- arr:

  Numeric vector or matrix from zarr read

- shape:

  Original zarr shape

- axes_names:

  Character vector of axis names (e.g., c("c", "y", "x"))

## Value

3D array (H x W x C) or NULL
