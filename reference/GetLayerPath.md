# Get Layer Path

Get the correct path for a layer/slot in V4 vs V5 format

## Usage

``` r
GetLayerPath(h5_group, slot_name, verbose = FALSE)
```

## Arguments

- h5_group:

  The HDF5 assay group

- slot_name:

  Name of the slot/layer

- verbose:

  Show progress messages

## Value

The path to the data, or NULL if not found
