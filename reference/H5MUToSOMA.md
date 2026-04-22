# Convert H5MU (MuData) to SOMA format

Reads an h5mu file and writes it to TileDB-SOMA format. Routes through a
Seurat intermediate.

## Usage

``` r
H5MUToSOMA(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .h5mu file

- dest:

  Output URI for the SOMA experiment

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readH5MU`](https://mianaz.github.io/scConvert/reference/readH5MU.md),
[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`SOMAToH5MU`](https://mianaz.github.io/scConvert/reference/SOMAToH5MU.md)
