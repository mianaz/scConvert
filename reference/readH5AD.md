# Load an AnnData H5AD file as a Seurat object

Direct conversion from H5AD format to Seurat object without intermediate
h5Seurat. Supports optional BPCells on-disk matrix loading for large
datasets that exceed available memory. When compiled C routines are
available, uses a fast native reader (typically 2-3x faster than the
pure-R path).

## Usage

``` r
readH5AD(
  file,
  assay.name = "RNA",
  use.bpcells = NULL,
  components = NULL,
  use.c = TRUE,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to H5AD file

- assay.name:

  Name for the primary assay (default: "RNA")

- use.bpcells:

  If not NULL, a directory path where BPCells will store the expression
  matrix on disk. Requires the BPCells package. The resulting Seurat
  object will reference the on-disk matrix instead of loading it into
  memory, enabling analysis of datasets larger than available RAM.

- components:

  Character vector of h5ad components to load. Default loads everything.
  Use `c("X")` for matrix-only (fastest), or any subset of
  `c("X", "obs", "var", "obsm", "obsp", "varp", "layers", "uns")`. The
  file path is stored in `misc[[".__h5ad_path__"]]` for deferred loading
  via
  [`scLoadMeta`](https://mianaz.github.io/scConvert/reference/scLoadMeta.md).

- use.c:

  Use compiled C reader when available (default: TRUE). Set to FALSE to
  force the pure-R hdf5r path.

- verbose:

  Show progress messages

## Value

A `Seurat` object. If `use.bpcells` is set, the count matrix is stored
on disk in BPCells format and the object uses minimal memory.
