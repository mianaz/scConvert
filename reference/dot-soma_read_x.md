# Read X matrix from a SOMA measurement

Read X matrix from a SOMA measurement

## Usage

``` r
.soma_read_x(ms, n_cells, n_features, cell_names, feature_names)
```

## Arguments

- ms:

  A SOMAMeasurement object

- n_cells:

  Expected number of cells

- n_features:

  Expected number of features

- cell_names:

  Character vector of cell names

- feature_names:

  Character vector of feature names

## Value

A sparse matrix (features x cells)
