# Convert SOMA to Zarr format

Reads a SOMA experiment and writes it to AnnData Zarr format. Routes
through a Seurat intermediate.

## Usage

``` r
SOMAToZarr(
  source,
  dest,
  measurement = "RNA",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- source:

  URI of the input SOMA experiment

- dest:

  Path for the output .zarr directory

- measurement:

  Name of the measurement to export (default: `"RNA"`)

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readSOMA`](https://mianaz.github.io/scConvert/reference/readSOMA.md),
[`writeZarr`](https://mianaz.github.io/scConvert/reference/writeZarr.md),
[`ZarrToSOMA`](https://mianaz.github.io/scConvert/reference/ZarrToSOMA.md)
