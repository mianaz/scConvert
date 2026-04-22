# Convert SOMA to h5ad format

Reads a SOMA experiment and writes it to h5ad format. Since SOMA is not
an HDF5 format, this conversion routes through a Seurat intermediate.

## Usage

``` r
SOMAToH5AD(
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

  Path for the output .h5ad file

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
[`writeH5AD`](https://mianaz.github.io/scConvert/reference/writeH5AD.md),
[`H5ADToSOMA`](https://mianaz.github.io/scConvert/reference/H5ADToSOMA.md)
