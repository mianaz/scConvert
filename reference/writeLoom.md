# Save a Seurat object to a Loom file

Export a Seurat object to Loom format (HDF5-based file format optimized
for storing annotated matrices). Loom files are compatible with the
loompy Python package and other tools in the bioinformatics community.
This format is useful for sharing data with Python-based analysis
workflows or archiving analysis results.

## Usage

``` r
writeLoom(object, filename, overwrite = FALSE, verbose = TRUE, ...)

as.loom(x, ...)

# Default S3 method
writeLoom(object, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'Seurat'
writeLoom(
  object,
  filename = paste0(Project(object = object), ".loom"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# Default S3 method
as.loom(x, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'H5File'
as.loom(x, ...)

# S3 method for class 'Seurat'
as.loom(
  x,
  filename = paste0(Project(object = x), ".loom"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object to save

- filename:

  Name of file to save the object to. If not provided, defaults to
  `<ProjectName>.loom`. The .loom extension is added automatically if
  not present.

- overwrite:

  Logical; if `TRUE`, overwrite an existing file. Default is `FALSE`.

- verbose:

  Logical; if `TRUE` (default), show progress updates

- ...:

  Arguments passed to other methods

- x:

  A Seurat object to convert

## Value

`writeLoom`: Invisibly returns the filename of the saved file

`as.loom`: A
[`loom`](https://mianaz.github.io/scConvert/reference/loom-class.md)
object

## Details

The Loom format organizes data as follows:

- `/matrix`: Main expression matrix (features x cells)

- `/row_attrs`: Feature/gene-level annotations

- `/col_attrs`: Cell/sample-level metadata (cell names, cluster
  assignments, etc.)

- `/layers`: Additional expression layers if present

When saving a Seurat object:

- Default assay data becomes the main matrix

- Cell metadata and feature annotations are preserved

- Dimensional reductions are stored in col_attrs

- The SEURAT_ASSAY attribute stores the assay name for roundtrip loading

## See also

[`readLoom`](https://mianaz.github.io/scConvert/reference/readLoom.md)
to load Loom files back as Seurat objects
[`writeH5Seurat`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md)
to save in h5Seurat format
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md)
for converting between formats [Loom
documentation](http://linnarssonlab.org/loompy/)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(scConvert)

# Create a Seurat object (or use an existing one)
seurat_obj <- CreateSeuratObject(counts = pbmc_small$RNA@counts)

# Save to Loom format
writeLoom(seurat_obj, filename = "my_data.loom")

# Save with overwrite if needed
writeLoom(seurat_obj, filename = "my_data.loom", overwrite = TRUE)

# Load it back
loaded_obj <- readLoom("my_data.loom")

# For sharing with Python tools
writeLoom(seurat_obj, filename = "data_for_python.loom")
# Use in Python with: adata = loompy.connect("data_for_python.loom")
} # }
```
