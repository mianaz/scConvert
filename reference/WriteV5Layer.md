# Write V5 Layer Data

Write matrix data in V5 layer format (sparse or dense)

## Usage

``` r
WriteV5Layer(
  h5_group,
  layer_name,
  matrix_data,
  features = NULL,
  verbose = FALSE
)
```

## Arguments

- h5_group:

  The HDF5 group to write to

- layer_name:

  Name of the layer (e.g., "counts", "data", "scale.data")

- matrix_data:

  The matrix data to write

- features:

  Character vector of feature names

- verbose:

  Show progress messages

## Value

Invisible NULL
