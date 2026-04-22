# Load a NanoString CosMx SMI bundle into a Seurat object

Thin R-only wrapper around
[`Seurat::LoadNanostring()`](https://satijalab.org/seurat/reference/ReadNanostring.html)
that validates the bundle layout, sets the scConvert spatial-technology
tag, and returns a Seurat object ready for conversion to h5ad / h5Seurat
/ Zarr via
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md).

## Usage

``` r
LoadCosMx(data.dir, fov = "cosmx", assay = "Nanostring", verbose = TRUE)
```

## Arguments

- data.dir:

  Path to the CosMx bundle directory

- fov:

  FOV (field-of-view) name attached to the resulting Seurat object
  (default "cosmx")

- assay:

  Seurat assay name (default "Nanostring")

- verbose:

  Print progress messages

## Value

A Seurat object with an FOV slot containing cell centroids, segmentation
boundaries, and transcript molecules

## Details

The bundle directory must contain the canonical NanoString CosMx CSV
files (`*exprMat_file*.csv`, `*metadata_file*.csv`,
`*fov_positions_file*.csv`) and optionally `*tx_file*.csv` and polygon
files. These are the files produced by NanoString's AtoMx export.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- LoadCosMx("/path/to/cosmx_nsclc/")
scConvert(seurat_obj, "cosmx_nsclc.h5ad")
} # }
```
