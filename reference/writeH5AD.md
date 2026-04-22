# Direct Seurat to H5AD Conversion

Converts a Seurat object directly to an H5AD file, handling the
intermediate h5Seurat file automatically. The intermediate file is
created in a temporary location and removed after conversion.

## Usage

``` r
writeH5AD(
  object,
  filename = NULL,
  assay = DefaultAssay(object = object),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  gzip = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object to convert

- filename:

  Output H5AD filename. If not specified, uses the project name with
  .h5ad extension.

- assay:

  Name of assay to convert. Default is the default assay.

- overwrite:

  Logical; overwrite existing file. Default FALSE.

- verbose:

  Logical; show progress messages. Default TRUE.

- standardize:

  Logical; convert Seurat metadata names to scanpy conventions. Default
  FALSE.

- gzip:

  Integer gzip compression level (0-9), or NULL to use the package
  default.

- ...:

  Additional arguments passed to writeH5Seurat and Convert.

## Value

Invisibly returns the path to the created H5AD file.

## Details

This function provides a convenient one-step conversion from Seurat
objects to H5AD format (used by Python's scanpy/anndata). Internally,
it:

1.  Saves the Seurat object to a temporary h5Seurat file

2.  Converts the h5Seurat file to H5AD format

3.  Removes the intermediate h5Seurat file

This is useful when you want to export data for Python analysis without
keeping the intermediate h5Seurat file.

## See also

[`writeH5Seurat`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md),
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md),
[`readH5AD`](https://mianaz.github.io/scConvert/reference/readH5AD.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(scConvert)

# Load or create a Seurat object
pbmc <- pbmc_small

# Convert directly to H5AD
writeH5AD(pbmc, filename = "pbmc.h5ad")

# The file can now be loaded in Python:
# import scanpy as sc
# adata = sc.read_h5ad("pbmc.h5ad")
} # }
```
