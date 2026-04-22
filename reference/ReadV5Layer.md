# Read V5 Layer Data

Read matrix data from V5 layer format

## Usage

``` r
ReadV5Layer(
  h5_group,
  layer_name,
  features = NULL,
  cells = NULL,
  verbose = FALSE
)
```

## Arguments

- h5_group:

  The HDF5 assay group

- layer_name:

  Name of the layer to read

- features:

  Expected feature names (for dimension validation)

- cells:

  Expected cell names (for dimension validation)

- verbose:

  Show progress messages

## Value

A matrix object, or NULL if layer doesn't exist
