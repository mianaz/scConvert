# Figure out which objects to load from an h5Seurat file

Figure out which objects to load from an h5Seurat file

## Usage

``` r
GetAssays(assays = NULL, index)

GetCommands(index, assays = NULL)

GetDimReducs(reductions, index, assays = NULL)

GetGraphs(graphs, index, assays = NULL)

GetImages(images, index, assays = NULL)

GetNeighbors(neighbors, index)
```

## Arguments

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

- index:

  An h5Seurat index
  ([`h5SI`](https://mianaz.github.io/scConvert/reference/h5SI.md))
  object

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

- neighbors:

  One of:

  - `NULL` (default): Load all neighbor information

  - A character vector with neighbor names

  - `FALSE`: Skip neighbors

## Value

`GetAssays`: A named list where each entry is a vector describing the
slots of an assay to load and the names are the assays to load

`GetCommands`: A vector of command log names that are derived from an
assay in `assay`

`GetDimReducs`: A vector of reduction names that are derived from an
assay in `assays` or global dimensional reductions

`GetGraphs`: A vector of graph names that are derived from an assay in
`assays`

`GetImages`: A vector of image names

`GetNeighbors`: A vector of neighbor names

## See also

[`readH5Seurat`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
