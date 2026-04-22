# Load a Stereo-seq GEF file into a Seurat object

Pure-R native reader for BGI Stereo-seq `.gef` and `.cellbin.gef` HDF5
files. No Python or stereopy dependency.

## Usage

``` r
LoadStereoSeqGef(file, bin_size = 50, assay = "Spatial", verbose = TRUE)
```

## Arguments

- file:

  Path to a `.gef` or `.cellbin.gef` file

- bin_size:

  Bin size to read for square-bin files (default 50, ignored for
  cellbin)

- assay:

  Seurat assay name (default "Spatial")

- verbose:

  Print progress messages

## Value

A Seurat object with spatial coordinates in `meta.data$spatial_x/y` and
provenance in `@misc$stereo_seq`

## Details

The function auto-detects square-bin versus cell-bin layouts from the
root `bin_type` attribute. Square-bin files are read at the requested
bin size (typically 50 or 100 for downstream analysis; bin1 is raw
DNB-level data and can be very large). Cell-bin files are read into a
gene x cell sparse matrix using the cell segmentation table.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- LoadStereoSeqGef("mouse_embryo.gef", bin_size = 50)
scConvert(seurat_obj, "mouse_embryo.h5ad")
} # }
```
