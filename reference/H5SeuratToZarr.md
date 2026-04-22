# Convert an h5seurat file to zarr (AnnData) format

Converts an h5seurat file to AnnData zarr format, translating the
h5seurat layout (levels/values factors, genes x cells CSC, reductions as
cell.embeddings) to AnnData conventions (categories/codes, cells x genes
CSR, obsm arrays). By default uses streaming that avoids a Seurat
intermediate.

## Usage

``` r
H5SeuratToZarr(
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

  Path to input .h5seurat file

- dest:

  Path for output .zarr directory

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
  then save as zarr.

## Value

Invisibly returns `dest`
