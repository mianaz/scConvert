# Get the proper HDF5 connection mode for writing depending on overwrite status

Get the proper HDF5 connection mode for writing depending on overwrite
status

## Usage

``` r
WriteMode(overwrite = FALSE)
```

## Arguments

- overwrite:

  Overwrite a file

## Value

`w` if `overwrite` else `w-`

## Examples

``` r
# \donttest{
scConvert:::WriteMode(TRUE)
#> [1] "w"
scConvert:::WriteMode(FALSE)
#> [1] "w-"
# }
```
