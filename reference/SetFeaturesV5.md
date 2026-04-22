# Set V5 feature index in HDF5 group

Stores feature names in the V5 location (meta.data/\_index) if not
already present.

## Usage

``` r
SetFeaturesV5(h5_group, features, verbose = FALSE)
```

## Arguments

- h5_group:

  HDF5 group to store features in

- features:

  Character vector of feature names

- verbose:

  Print diagnostic messages

## Value

NULL (invisible)
