# Create an HDF5 object path

Create an HDF5 object path

## Usage

``` r
H5Path(..., collapse = NULL)
```

## Arguments

- ...:

  one or more R objects, to be converted to character vectors.

- collapse:

  an optional character string to separate the results. Not
  [`NA_character_`](https://rdrr.io/r/base/NA.html). When `collapse` is
  a string, the result is always a string
  ([`character`](https://rdrr.io/r/base/character.html) of length 1).

## Value

A character vector with path ready for accessing data in an HDF5 file or
group
