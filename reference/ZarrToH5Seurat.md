# Convert a zarr (AnnData) store to h5seurat format

Converts an AnnData zarr store to h5seurat format, translating AnnData
conventions (categories/codes 0-based, cells x genes CSR, obsm arrays)
to h5seurat layout (levels/values 1-based, genes x cells CSC,
reductions/cell.embeddings). By default uses streaming that avoids a
Seurat intermediate.

## Usage

``` r
ZarrToH5Seurat(
  source,
  dest,
  assay = "RNA",
  overwrite = FALSE,
  gzip = 4L,
  verbose = TRUE,
  stream = TRUE
)
```

## Arguments

- source:

  Path to input .zarr directory

- dest:

  Path for output .h5seurat file

- assay:

  Assay name (default "RNA")

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly. If FALSE, load as Seurat
  then save as h5seurat.

## Value

Invisibly returns `dest`
