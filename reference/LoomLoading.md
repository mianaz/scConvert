# Loom-file Loading

Version-specific loom-file loading functions

## Usage

``` r
LoadLoom0.1(
  file,
  assay = NULL,
  cells = "col_atts/CellID",
  features = "row_attrs/Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE
)

LoadLoom3.0(
  file,
  assay = NULL,
  cells = "col_attrs/CellID",
  features = "row_attrs/Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE
)
```

## Arguments

- assay:

  Name of assay to store expression data as. If `NULL`, will search for
  an HDF5 attribute named `SEURAT_ASSAY` or dataset at
  `/attrs/SEURAT_ASSAY` for the assay name. Defaults to `"RNA"` if not
  found.

- cells:

  Name of dataset in `/col_attrs` containing cell/sample names. If not
  found, uses column indices as cell names.

- features:

  Name of dataset in `/row_attrs` containing feature/gene names. If not
  found, uses row indices as feature names.

- normalized:

  Name of layer in `/layers` to load as normalized expression data;
  special value `"/matrix"` loads the main matrix as normalized data
  instead of counts

- scaled:

  Name of layer in `/layers` to load as scaled data

- filter:

  Logical; if `TRUE`, keep only cells and features marked as valid using
  `/col_attrs/Valid` and `/row_attrs/Valid` attributes

- verbose:

  Logical; if `TRUE` (default), show progress updates

## Value

A `Seurat` object containing the expression matrix and metadata from the
Loom file

## Details

`readLoom` will try to automatically fill slots of a `Seurat` object
based on data presence or absence in a given loom file. This method
varies by loom specification version. For version-specific details, see
sections below

## Loom 0.1 Loading

Loading data from loom files less than version 3.0.0 is not currently
supported

## Loom 3.0.0 Loading

For loom files version 3.0.0 and higher, `LoadLoom3.0` provides
comprehensive loading with support for filtering, dimensional
reductions, and cell graphs.
