# Convert single-cell datasets between formats

Universal converter between single-cell file formats and object types.
Supports arbitrary source/destination pairs by routing through Seurat as
a hub format. Direct HDF5-level paths (h5ad \<-\> h5seurat) are used
when available for memory efficiency.

## Usage

``` r
scConvert(
  source,
  dest,
  assay,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
)

# S3 method for class 'character'
scConvert(
  source,
  dest,
  assay,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
)

# S3 method for class 'H5File'
scConvert(
  source,
  dest = "h5seurat",
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'h5Seurat'
scConvert(
  source,
  dest = "h5ad",
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
)

# S3 method for class 'Seurat'
scConvert(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
)

# S3 method for class 'loom'
scConvert(
  source,
  dest,
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
)

# S3 method for class 'SingleCellExperiment'
scConvert(
  source,
  dest,
  assay = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  standardize = FALSE,
  ...
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

- ...:

  Arguments passed to specific conversion methods

## Value

For file destinations, invisibly returns the destination filename. For
`dest = "sce"`, returns a SingleCellExperiment object.

## Details

**Supported Formats:**

- **R objects**: Seurat, SingleCellExperiment (requires
  SingleCellExperiment), loom (R6 connection)

- **File formats**: h5seurat, h5ad, h5mu, loom, rds

Any source format can be converted to any destination format.
Conversions without a direct path go through Seurat as a universal hub:
`Source -> Seurat -> Destination`.

**Direct Paths** (memory-efficient, no full dataset loading):

- `h5ad <-> h5seurat`: Direct HDF5-level copy

**Key Features:**

- Preserves expression matrices, metadata, and dimensional reductions

- For Visium/spatial data: reconstructs images with scale factors

- Handles multiple data layers (V5 compatibility)

## AnnData/H5AD to h5Seurat

The AnnData/H5AD to h5Seurat conversion will try to automatically fill
in datasets based on data presence. It works in the following manner:

### Expression data

The expression matrices `counts`, `data`, and `scale.data` are filled by
`/X` and `/raw/X` in the following manner:

- `counts` will be filled with `/raw/X` if present; otherwise, it will
  be filled with `/X`

- `data` will be filled with `/raw/X` if `/raw/X` is present and `/X` is
  dense; otherwise, it will be filled with `/X`

- `scale.data` will be filled with `/X` if it dense; otherwise, it will
  be empty

Feature names are taken from the feature-level metadata

### Feature-level metadata

Feature-level metadata is added to the `meta.features` datasets in each
assay. Feature names are taken from the dataset specified by the
“\_index” attribute, the “\_index” dataset, or the “index” dataset, in
that order. Metadata is populated with `/raw/var` if present, otherwise
with `/var`; if both `/raw/var` and `/var` are present, then
`meta.features` will be populated with `/raw/var` first, then `/var`
will be added to it. For columns present in both `/raw/var` and `/var`,
the values in `/var` will be used instead. **Note**: it is possible for
`/var` to have fewer features than `/raw/var`; if this is the case, then
only the features present in `/var` will be overwritten, with the
metadata for features *not* present in `/var` remaining as they were in
`/raw/var` or empty

### Cell-level metadata

Cell-level metadata is added to `meta.data`; the row names of the
metadata (as determined by the value of the “\_index” attribute, the
“\_index” dataset, or the “index” dataset, in that order) are added to
the “cell.names” dataset instead. If the “\_\_categories” dataset is
present, each dataset within “\_\_categories” will be stored as a factor
group. Cell-level metadata will be added as an HDF5 group unless factors
are **not** present and the `scConvert.dtypes.dataframe_as_group` option
is `FALSE`

### Dimensional reduction information:

Cell embeddings are taken from `/obsm`; dimensional reductions are named
based on their names from `obsm` by removing the preceding “X\_”.For
example, if a dimensional reduction is named “X_pca” in `/obsm`, the
resulting dimensional reduction information will be named “pca”. The key
will be set to one of the following:

- “PC\_” if “pca” is present in the dimensional reduction name
  (`grepl("pca", reduction.name, ignore.case = TRUE)`)

- “tSNE\_” if “tsne” is present in the dimensional reduction name
  (`grepl("tsne", reduction.name, ignore.case = TRUE)`)

- `reduction.name_` for all other reductions

Remember that the preceding “X\_” will be removed from the reduction
name before converting to a key. Feature loadings are taken from `/varm`
and placed in the associated dimensional reduction. The dimensional
reduction is determine from the loadings name in `/varm`:

- “PCs” will be added to a dimensional reduction named “pca”

- All other loadings in `/varm` will be added to a dimensional reduction
  named `tolower(loading)` (eg. a loading named “ICA” will be added to a
  dimensional reduction named “ica”)

If a dimensional reduction cannot be found according to the rules above,
the loading will not be taken from the AnnData/H5AD file. Miscellaneous
information will be taken from `/uns/reduction` where `reduction` is the
name of the reduction in `/obsm` without the preceding “X\_”; if no
dimensional reduction information present, then miscellaneous
information will not be taken from the AnnData/H5AD file. Standard
deviations are taken from a dataset `/uns/reduction/variance`; the
variances will be converted to standard deviations and added to the
`stdev` dataset of a dimensional reduction

### Nearest-neighbor graph

If a nearest neighbor graph is present in `/uns/neighbors/distances`, it
will be added as a graph dataset in the h5Seurat file and associated
with `assay`; if a value is present in `/uns/neighbors/params/method`,
the name of the graph will be `assay_method`, otherwise, it will be
`assay_anndata`

### Miscellaneous information

All groups and datasets from `/uns` will be copied to `misc` in the
h5Seurat file except for the following:

- Any group or dataset named the same as a dimensional reduction (eg.
  `/uns/pca`)

- `/uns/neighbors`

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

## See also

[`writeH5AD`](https://mianaz.github.io/scConvert/reference/writeH5AD.md)
for direct Seurat to h5ad convenience function
[`writeH5Seurat`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md)
to save Seurat objects
[`readH5Seurat`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
to load h5Seurat files
[`readH5AD`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
to directly load h5ad files
[`readH5MU`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
to load h5mu files
[`scConnect`](https://mianaz.github.io/scConvert/reference/scConnect.md)
to establish file connections

## Examples

``` r
if (FALSE) { # \dontrun{
library(scConvert)
library(Seurat)

# --- Any format to any format ---
scConvert("data.h5ad", dest = "data.h5seurat")     # h5ad -> h5seurat
scConvert("data.h5ad", dest = "data.rds")           # h5ad -> RDS
scConvert("data.h5ad", dest = "data.loom")           # h5ad -> loom
scConvert("data.h5mu", dest = "data.h5ad")           # h5mu -> h5ad
scConvert("data.loom", dest = "data.h5seurat")       # loom -> h5seurat
scConvert("data.rds",  dest = "data.h5ad")           # RDS -> h5ad

# --- From R objects ---
scConvert(seurat_obj, dest = "output.h5ad")          # Seurat -> h5ad
scConvert(seurat_obj, dest = "output.loom")           # Seurat -> loom
scConvert(seurat_obj, dest = "output.rds")            # Seurat -> RDS
scConvert(sce_obj, dest = "output.h5ad")              # SCE -> h5ad
sce <- scConvert(seurat_obj, dest = "sce")            # Seurat -> SCE
} # }
```
