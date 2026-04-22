# Load a Seurat object from an h5Seurat file

Load a previously saved Seurat object from an h5Seurat file. This
function supports flexible loading options, allowing you to load only
the components you need (e.g., specific assays, reductions) to minimize
memory usage on large datasets.

## Usage

``` r
readH5Seurat(file, ...)

# S3 method for class 'character'
readH5Seurat(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
)

# S3 method for class 'H5File'
readH5Seurat(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
readH5Seurat(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
as.Seurat(
  x,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = TRUE,
  tools = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- file, x:

  Name of an h5Seurat file path (character) or connected h5Seurat file
  to load

- ...:

  Arguments passed to other methods

- assays:

  One of:

  - `NULL` (default): Load all assays

  - A character vector with names of assays to load (e.g.,
    `c("RNA", "ADT")`)

  - A character vector specifying which data layers to load for all
    assays: `c("counts", "data")` loads only the counts and data layers,
    skipping scale.data

  - A named list for fine-grained control, e.g.,
    `list(RNA = "data", ADT = c("data", "scale.data"))`

- reductions:

  One of:

  - `NULL` (default): Load all reductions (PCA, UMAP, etc.)

  - A character vector with names of specific reductions (e.g.,
    `c("pca", "umap")`)

  - `NA`: Load only global (assay-independent) reductions

  - `FALSE`: Skip loading all reductions

  **Note**: Only reductions associated with a loaded assay or marked as
  global will be loaded.

- graphs:

  One of:

  - `NULL` (default): Load all graphs

  - A character vector with specific graph names (e.g.,
    `c("RNA_snn", "ADT_snn")`)

  - `FALSE`: Skip loading graphs

  **Note**: Only graphs associated with loaded assays will be available.

- neighbors:

  One of:

  - `NULL` (default): Load all neighbor information

  - A character vector with neighbor names

  - `FALSE`: Skip neighbors

- images:

  One of:

  - `NULL` (default): Load all images (for spatial experiments)

  - A character vector with image names

  - `NA`: Load only global images

  - `FALSE`: Skip images

- meta.data:

  Logical; if `TRUE` (default), load cell metadata

- commands:

  Logical; if `TRUE` (default), load command history. Commands are only
  loaded if their associated assays are loaded.

- misc:

  Logical; if `TRUE` (default when all assays loaded), load
  miscellaneous data

- tools:

  Logical; if `TRUE` (default when all assays loaded), load
  tool-specific information

- verbose:

  Logical; if `TRUE` (default), show progress messages

## Value

A `Seurat` object containing the requested components

## Details

The h5Seurat format is highly flexible for selective loading. This is
particularly useful when:

- Working with very large datasets where loading everything would exceed
  memory

- You only need specific assays or reductions for downstream analysis

- You want to quickly inspect object structure without full data loading

## Seurat V5 Layer Support

For Seurat V5 objects with multiple layers, you can selectively load
layers per assay. For example, use `assays = list(RNA = "data")` to load
only the normalized expression layer, skipping raw counts and scaled
data.

## See also

[`writeH5Seurat`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md)
to save a Seurat object to h5Seurat format
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md)
to convert to other formats

## Examples

``` r
if (FALSE) { # \dontrun{
library(scConvert)

# Load entire h5Seurat file
seurat_obj <- readH5Seurat("data.h5seurat")

# Load only specific assays
seurat_obj <- readH5Seurat("data.h5seurat", assays = c("RNA", "ADT"))

# Load only specific data layers (memory-efficient for large files)
seurat_obj <- readH5Seurat("data.h5seurat", assays = c("data"))  # Only normalized expression

# Load specific assays with different layers
seurat_obj <- readH5Seurat(
  "data.h5seurat",
  assays = list(RNA = c("data", "scale.data"), ADT = "data")
)

# Load without reductions (faster)
seurat_obj <- readH5Seurat("data.h5seurat", reductions = FALSE)

# Load UMAP and PCA reductions only
seurat_obj <- readH5Seurat("data.h5seurat", reductions = c("umap", "pca"))

# Load spatial data without graphs (for Visium experiments)
seurat_obj <- readH5Seurat("visium.h5seurat", images = TRUE, graphs = FALSE)
} # }
```
