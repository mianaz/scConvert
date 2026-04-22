# Convert Loom to SOMA format

Reads a Loom file and writes it to TileDB-SOMA format. Routes through a
Seurat intermediate.

## Usage

``` r
LoomToSOMA(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .loom file

- dest:

  Output URI for the SOMA experiment

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readLoom`](https://mianaz.github.io/scConvert/reference/readLoom.md),
[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`SOMAToLoom`](https://mianaz.github.io/scConvert/reference/SOMAToLoom.md)
