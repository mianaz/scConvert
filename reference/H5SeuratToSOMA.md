# Convert h5Seurat to SOMA format

Reads an h5Seurat file and writes it to TileDB-SOMA format. Routes
through a Seurat intermediate.

## Usage

``` r
H5SeuratToSOMA(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .h5seurat file

- dest:

  Output URI for the SOMA experiment

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readH5Seurat`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md),
[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`SOMAToH5Seurat`](https://mianaz.github.io/scConvert/reference/SOMAToH5Seurat.md)
