# Determine the margin to use for a dataset

Determine the margin to use for a dataset

## Usage

``` r
GetMargin(dims, MARGIN = getOption(x = "scConvert.chunking.MARGIN"))
```

## Arguments

- dims:

  Dimensions of a dataset

- MARGIN:

  Either an integer value contained within `1:length(x = dims)` or one
  of the possible values of the `scConvert.chunking.MARGIN` option

## Value

An integer value with the `MARGIN`

## See also

[`scConvert-package`](https://mianaz.github.io/scConvert/reference/scConvert-package.md)
for package options

## Examples

``` r
# \donttest{
scConvert:::GetMargin(c(4, 10))
#> [1] 2
# }
```
