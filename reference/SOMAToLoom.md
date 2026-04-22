# Convert SOMA to Loom format

Reads a SOMA experiment and writes it to Loom format. Routes through a
Seurat intermediate.

## Usage

``` r
SOMAToLoom(
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

  Path for the output .loom file

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
[`writeLoom`](https://mianaz.github.io/scConvert/reference/writeLoom.md),
[`LoomToSOMA`](https://mianaz.github.io/scConvert/reference/LoomToSOMA.md)
