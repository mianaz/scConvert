# Save a Seurat object to an h5Seurat file

Save a Seurat object to the efficient HDF5-based h5Seurat format. This
format is optimized for large, complex single-cell datasets including
multi-modal and spatial data. h5Seurat files can be rapidly converted to
other formats like h5ad (AnnData) or h5mu (MuData) for interoperability
with Python tools.

## Usage

``` r
writeH5Seurat(
  object,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  gzip = NULL,
  ...
)

as.h5Seurat(x, ...)

# Default S3 method
writeH5Seurat(object, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'Seurat'
writeH5Seurat(
  object,
  filename = paste0(Project(object = object), ".h5Seurat"),
  overwrite = FALSE,
  verbose = TRUE,
  gzip = NULL,
  ...
)

# Default S3 method
as.h5Seurat(x, filename, overwrite = FALSE, verbose = TRUE, ...)

# S3 method for class 'H5File'
as.h5Seurat(x, ...)

# S3 method for class 'Seurat'
as.h5Seurat(
  x,
  filename = paste0(Project(object = x), ".h5seurat"),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object, x:

  An object (typically a Seurat object)

- filename:

  Name of file to save the object to. If not provided, defaults to
  `<ProjectName>.h5seurat`. The .h5seurat extension is added
  automatically if not present.

- overwrite:

  Logical; if `TRUE`, overwrite an existing file with the same name.
  Default is `FALSE`.

- verbose:

  Show progress updates during save. Default is `TRUE`.

- gzip:

  Integer gzip compression level (0-9), or NULL to use the package
  default.

- ...:

  Arguments passed to other methods

## Value

`writeH5Seurat`: Invisibly returns the filename of the saved file

`as.h5Seurat`: An
[`h5Seurat`](https://mianaz.github.io/scConvert/reference/h5Seurat-class.md)
object

## Details

The h5Seurat format stores:

- All assays with their layers (counts, data, scale.data, etc.) for
  Seurat V5

- Dimensional reductions (PCA, UMAP, etc.)

- Nearest-neighbor graphs and similarity graphs

- Spatial images and coordinates (for spatial experiments)

- Cell metadata and feature annotations

- Cell identity classes

- Command history

- Miscellaneous and tool-specific data

The h5Seurat format is particularly useful for:

- Storing large datasets efficiently with HDF5 compression

- Rapid conversion to Python formats (h5ad, h5mu)

- Multi-modal and spatial transcriptomics experiments

- Preserving all Seurat V5 layer information

## Seurat V5 Layer Support

When saving Seurat V5 objects with multiple layers (e.g., counts, data,
scale.data), all layers are preserved and can be selectively loaded
using
[`readH5Seurat`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md).

## See also

[`readH5Seurat`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
to load a saved h5Seurat file `as.h5Seurat` for direct conversion
without object assignment
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md)
for converting to other formats (h5ad, h5mu)

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(scConvert)

# Create a simple example Seurat object
seurat_obj <- CreateSeuratObject(counts = GetAssayData(pbmc_small))

# Save to h5Seurat format
writeH5Seurat(seurat_obj, filename = "my_data.h5seurat")

# Save with overwrite if file already exists
writeH5Seurat(seurat_obj, filename = "my_data.h5seurat", overwrite = TRUE)

# Load the saved file back
seurat_obj <- readH5Seurat("my_data.h5seurat")

# For multimodal data (e.g., CITE-seq)
# writeH5Seurat automatically saves all assays
writeH5Seurat(citeseq_obj, filename = "multimodal_data.h5seurat")

# For spatial data (e.g., Visium)
writeH5Seurat(visium_obj, filename = "spatial_visium.h5seurat")
} # }
```
