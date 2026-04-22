# Convert SOMA to H5MU (MuData) format

Reads a SOMA experiment and writes it to h5mu format. Routes through a
Seurat intermediate. Best suited for SOMA experiments with multiple
measurements (assays).

## Usage

``` r
SOMAToH5MU(
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

  Path for the output .h5mu file

- measurement:

  Name of the primary measurement (default: `"RNA"`)

- overwrite:

  If `TRUE`, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`

## See also

[`readSOMA`](https://mianaz.github.io/scConvert/reference/readSOMA.md),
[`writeH5MU`](https://mianaz.github.io/scConvert/reference/writeH5MU.md),
[`H5MUToSOMA`](https://mianaz.github.io/scConvert/reference/H5MUToSOMA.md)
