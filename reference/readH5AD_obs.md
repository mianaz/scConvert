# Read observation metadata from an H5AD file

Read the cell-level metadata (obs) from an H5AD file and return it as a
data.frame. Handles both categorical (factor) and standard data types.

## Usage

``` r
readH5AD_obs(file)
```

## Arguments

- file:

  Path to an H5AD file

## Value

A `data.frame` with cell metadata, where row names are cell barcodes
from the `_index` dataset and columns correspond to observation
annotations
