# Write a Seurat object to TileDB-SOMA format

Writes a Seurat object to a TileDB-SOMA experiment. Uses the
`tiledbsoma::write_soma()` function from the tiledbsoma package, which
handles the full conversion of assays, metadata, embeddings, and graphs.

## Usage

``` r
writeSOMA(object, uri, measurement = NULL, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- object:

  A Seurat object

- uri:

  Output URI (local path or cloud URI such as `tiledb://...` or
  `s3://...`)

- measurement:

  Name for the primary measurement. If `NULL` (default), uses
  `DefaultAssay(object)`.

- overwrite:

  If `TRUE`, remove any existing SOMA at `uri` before writing. Default
  is `FALSE`.

- verbose:

  Show progress messages

## Value

Invisibly returns `uri`

## See also

[`readSOMA`](https://mianaz.github.io/scConvert/reference/readSOMA.md),
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(scConvert)
library(Seurat)

# Write a Seurat object to SOMA
writeSOMA(seurat_obj, uri = "output_soma")

# Overwrite existing
writeSOMA(seurat_obj, uri = "output_soma", overwrite = TRUE)
} # }
```
