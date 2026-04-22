# Convert Zarr to SOMA format

Reads an AnnData Zarr store and writes it to TileDB-SOMA format. Routes
through a Seurat intermediate.

## Usage

``` r
ZarrToSOMA(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .zarr directory

- dest:

  Output URI for the SOMA experiment

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readZarr`](https://mianaz.github.io/scConvert/reference/readZarr.md),
[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`SOMAToZarr`](https://mianaz.github.io/scConvert/reference/SOMAToZarr.md)
