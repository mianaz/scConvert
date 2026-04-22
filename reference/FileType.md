# Determine a filetype based on its extension

Determine a filetype based on its extension

## Usage

``` r
FileType(file)
```

## Arguments

- file:

  Name of file

## Value

The extension, all lowercase

## Examples

``` r
# \donttest{
scConvert:::FileType('pbmc3k.h5Seurat')
#> [1] "h5seurat"
scConvert:::FileType('h5ad')
#> [1] "h5ad"
# }
```
