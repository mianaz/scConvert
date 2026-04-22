# Direct Seurat-to-h5ad write (single-pass, no intermediate h5Seurat)

Writes a Seurat object directly to an h5ad file, bypassing the two-pass
Seurat -\> h5Seurat -\> h5ad pipeline. This avoids 3x I/O overhead and
provides significant speedup for in-memory Seurat objects.

## Usage

``` r
DirectSeuratToH5AD(
  object,
  filename,
  assay = DefaultAssay(object = object),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE
)
```

## Arguments

- object:

  A Seurat object

- filename:

  Output h5ad file path

- assay:

  Name of the assay to write (default: DefaultAssay)

- overwrite:

  Overwrite existing file (default: FALSE)

- verbose:

  Print progress messages (default: TRUE)

- standardize:

  Convert Seurat metadata names to scanpy conventions

## Value

Invisibly returns the output filename
