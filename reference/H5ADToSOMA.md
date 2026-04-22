# Convert h5ad to SOMA format

Reads an h5ad file and writes it to TileDB-SOMA format. Since SOMA is
not an HDF5 format, this conversion routes through a Seurat
intermediate.

## Usage

``` r
H5ADToSOMA(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .h5ad file

- dest:

  Output URI for the SOMA experiment

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readH5AD`](https://mianaz.github.io/scConvert/reference/readH5AD.md),
[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`SOMAToH5AD`](https://mianaz.github.io/scConvert/reference/SOMAToH5AD.md)
