# Convert SOMA to h5Seurat format

Reads a SOMA experiment and writes it to h5Seurat format. Routes through
a Seurat intermediate.

## Usage

``` r
SOMAToH5Seurat(
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

  Path for the output .h5seurat file

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
[`writeH5Seurat`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md),
[`H5SeuratToSOMA`](https://mianaz.github.io/scConvert/reference/H5SeuratToSOMA.md)
