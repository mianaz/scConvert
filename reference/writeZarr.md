# Save a Seurat object as an AnnData Zarr store

Writes a Seurat object to AnnData Zarr v2 format (.zarr directory). The
output is readable by Python `anndata.read_zarr()`, scanpy, and any tool
supporting the AnnData on-disk specification.

## Usage

``` r
writeZarr(
  object,
  filename,
  assay = DefaultAssay(object),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object

- filename:

  Path for the output .zarr directory

- assay:

  Name of assay to export (default: `DefaultAssay(object)`)

- overwrite:

  If TRUE, overwrite existing zarr store

- verbose:

  Show progress messages

- ...:

  Additional arguments (currently unused)

## Value

Invisibly returns `filename`
