# Convert h5Seurat files to H5AD files

Convert h5Seurat files to H5AD files

## Usage

``` r
H5SeuratToH5AD(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE
)
```

## Arguments

- source:

  Source dataset: a Seurat object, SingleCellExperiment, loom
  connection, filename path, or H5File connection

- dest:

  Name/path of destination file or format. Supported formats: h5seurat,
  h5ad, h5mu, loom, rds. Also accepts `"sce"` to return a
  SingleCellExperiment object (in-memory, no file created).

- assay:

  For h5Seurat -\> other formats: name of assay to convert. For other
  formats -\> h5Seurat: name to assign to the assay. Default is "RNA".

- overwrite:

  Logical; if `TRUE`, overwrite an existing destination file. Default is
  `FALSE`.

- verbose:

  Logical; if `TRUE` (default), show progress updates

- standardize:

  Logical; if `TRUE`, convert Seurat-style metadata column names to
  scanpy/AnnData conventions when converting to h5ad format. For
  example, `nCount_RNA` becomes `n_counts`, `nFeature_RNA` becomes
  `n_genes`. Only applicable for conversions to h5ad format. Default is
  `FALSE`.

## Value

Returns a handle to `dest` as an
[`H5File`](http://hhoeflin.github.io/hdf5r/reference/H5File-class.md)
object

## h5Seurat to AnnData/H5AD

The h5Seurat to AnnData/H5AD conversion will try to automatically fill
in datasets based on data presence. Data presense is determined by the
h5Seurat index (`source$index()`). It works in the following manner:

### Assay data

- `X` will be filled with `scale.data` if `scale.data` is present;
  otherwise, it will be filled with `data`

- `var` will be filled with `meta.features` **only** for the features
  present in `X`; for example, if `X` is filled with `scale.data`, then
  `var` will contain only features that have been scaled

- `raw.X` will be filled with `data` if `X` is filled with `scale.data`;
  otherwise, it will be filled with `counts`. If `counts` is not
  present, then `raw` will not be filled

- `raw.var` will be filled with `meta.features` with the features
  present in `raw.X`; if `raw.X` is not filled, then `raw.var` will not
  be filled

### Cell-level metadata

Cell-level metadata is added to `obs`

### Dimensional reduction information

Only dimensional reductions associated with `assay` or marked as
[global](https://satijalab.github.io/seurat-object/reference/IsGlobal.html)
will be transfered to the H5AD file. For every reduction `reduc`:

- cell embeddings are placed in `obsm` and renamed to `X_reduc`

- feature loadings, if present, are placed in `varm` and renamed to
  either “PCs” if `reduc` is “pca” otherwise `reduc` in all caps

For example, if `reduc` is “ica”, then cell embeddings will be “X_ica”
in `obsm` and feature loaodings, if present, will be “ICA” in `varm`

### Nearest-neighbor graphs

If a nearest-neighbor graph is associated with `assay`, it will be added
to `uns/neighbors/distances`; if more than one graph is present, then
**only** the last graph according to the index will be added.

### Layers

Data from other assays can be added to `layers` if they have the same
shape as `X` (same number of cells and features). To determine this, the
shape of each alternate assays's `scale.data` and `data` slots are
determined. If they are the same shape as `X`, then that slot
(`scale.data` is given priority over `data`) will be added as a layer
named the name of the assay (eg. “SCT”). In addition, the features names
will be added to `var` as `assay_features` (eg. “SCT_features”).
