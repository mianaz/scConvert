# Get feature names with V5 fallback

Retrieves feature names from an HDF5 group, trying V5 location first
(meta.data/\_index), then falling back to direct features dataset.

## Usage

``` r
GetFeaturesV5Safe(h5_group, verbose = FALSE)
```

## Arguments

- h5_group:

  HDF5 group containing feature information

- verbose:

  Print diagnostic messages

## Value

Character vector of feature names, or NULL if not found
