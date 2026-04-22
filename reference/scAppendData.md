# Append data from an h5Seurat file to a preexisting [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html) object

Append data from an h5Seurat file to a preexisting
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
object

## Usage

``` r
scAppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'character'
scAppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'H5File'
scAppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
scAppendData(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = "commands",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- file:

  Name of an h5Seurat file path (character) or connected h5Seurat file

- object:

  A
  [`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
  object to append data to

- assays:

  One of:

  - A character vector with names of assays

  - A character vector with one or more of `counts`, `data`,
    `scale.data` describing which slots of **all assays** to load

  - A named list where each entry is either the name of an assay or a
    vector describing which slots (described above) to take from which
    assay

  - `NULL` for all assays

  - `FALSE` for no assays

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

- images:

  One of:

  - `NULL` (default): Load all images (for spatial experiments)

  - A character vector with image names

  - `NA`: Load only global images

  - `FALSE`: Skip images

- extras:

  Extra information to load; supports any combination of the following
  values:

  “commands”

  :   Load command logs. If `overwrite = TRUE`, replaces existing
      command logs

- overwrite:

  Overwrite existing data in `object` with data from `file`

- verbose:

  Logical; if `TRUE` (default), show progress messages

- ...:

  Arguments passed to other methods

## Value

`object` with the extra data requested
