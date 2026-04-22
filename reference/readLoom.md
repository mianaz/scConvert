# Load a Loom file as a Seurat object

Load gene expression data from Loom files (HDF5-based file format for
storing annotated matrices) into a Seurat object. Loom files are
commonly created by Python-based analysis tools and offer an alternative
to H5AD format for storing single-cell data.

## Usage

``` r
readLoom(
  file,
  assay = NULL,
  cells = "CellID",
  features = "Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE,
  ...
)

# S3 method for class 'character'
readLoom(file, ...)

# S3 method for class 'H5File'
readLoom(file, ...)

# S3 method for class 'loom'
readLoom(file, ...)

# S3 method for class 'loom'
as.Seurat(
  x,
  assay = NULL,
  cells = "CellID",
  features = "Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE,
  ...
)
```

## Arguments

- file, x:

  Name of Loom file (character path) or a `loom` object connection to
  load data from

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

- ...:

  Arguments passed to other methods

## Value

A `Seurat` object containing the expression matrix and metadata from the
Loom file

## Details

The Loom format stores data in a specific HDF5 structure:

- `/matrix`: Main expression matrix (features x cells)

- `/row_attrs`: Feature-level metadata (gene names, coordinates, etc.)

- `/col_attrs`: Cell/sample-level metadata

- `/layers`: Additional expression layers (normalized, scaled, etc.)

## Loom 0.1 Loading

Loading data from loom files less than version 3.0.0 is not currently
supported

## Loom 3.0.0 Loading

For loom files version 3.0.0 and higher, `LoadLoom3.0` provides
comprehensive loading with support for filtering, dimensional
reductions, and cell graphs.

## See also

[`writeLoom`](https://mianaz.github.io/scConvert/reference/writeLoom.md)
to save Seurat objects to Loom format
[`readH5Seurat`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
to load h5Seurat files
[`readH5AD`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
to load h5ad files [Loom
documentation](http://linnarssonlab.org/loompy/) [Loom file
conventions](http://linnarssonlab.org/loompy/conventions/index.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(scConvert)

# Load a basic Loom file
seurat_obj <- readLoom("data.loom")

# Load with specific assay name
seurat_obj <- readLoom("data.loom", assay = "RNA")

# Load and filter to valid cells/features only
seurat_obj <- readLoom("data.loom", filter = TRUE)

# Load specific layers
seurat_obj <- readLoom(
  "data.loom",
  normalized = "normalized",
  scaled = "scaled"
)
} # }
```
