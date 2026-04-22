# Convert FOV object to an h5ad-compatible list structure

Produced structure is consumed by
[`WriteFOVToH5AD`](https://mianaz.github.io/scConvert/reference/WriteFOVToH5AD.md).
Users do not normally call this directly.

## Usage

``` r
ConvertFOVToH5AD(fov_obj, verbose = TRUE)
```

## Arguments

- fov_obj:

  FOV object

- verbose:

  Print progress messages

## Value

named list with `centroids`, `segmentation`, `molecules`, `technology`,
`fov_name`
