# Read a TileDB-SOMA experiment as a Seurat object

Reads a TileDB-SOMA experiment (local or remote) and returns a Seurat
object. Supports cell/feature filtering via value queries, multiple
measurements (assays), embeddings (obsm), and neighbor graphs (obsp).

## Usage

``` r
readSOMA(
  uri,
  measurement = "RNA",
  obs_query = NULL,
  var_query = NULL,
  obs_column_names = NULL,
  var_column_names = NULL,
  verbose = TRUE
)
```

## Arguments

- uri:

  URI to a SOMA experiment. Can be a local path or a cloud URI (e.g.,
  `tiledb://...` or `s3://...`).

- measurement:

  Name of the measurement to use as the default assay (default:
  `"RNA"`).

- obs_query:

  Optional value filter string for cells, passed to
  `SOMADataFrame$read(value_filter = ...)`. For example,
  `"cell_type == 'T cell'"`.

- var_query:

  Optional value filter string for features.

- obs_column_names:

  Optional character vector of obs column names to read. If `NULL`
  (default), all columns are read.

- var_column_names:

  Optional character vector of var column names to read. If `NULL`
  (default), all columns are read.

- verbose:

  Show progress messages

## Value

A `Seurat` object

## Details

TileDB-SOMA is the cloud-native format underlying CELLxGENE Census (61M+
cells). This function requires the tiledbsoma R package.

The reader performs the following steps:

1.  Opens the SOMA experiment at `uri`

2.  Reads obs (cell metadata) with optional filtering

3.  Reads var (feature metadata) from the specified measurement

4.  Reads the X matrix (`data` layer) as a sparse matrix

5.  Reads any obsm entries as dimensional reductions

6.  Reads any obsp entries as neighbor graphs

7.  Reads additional measurements as extra Seurat assays

Cell names are determined from the following columns in priority order:
`obs_id`, `_index`, `index`, `barcode`, `cell_id`. If none are present,
synthetic names (`cell_1`, `cell_2`, ...) are generated.

## CELLxGENE Census workflow

CELLxGENE Census exposes 60M+ cells as a single SOMA experiment hosted
on public S3 (`s3://cellxgene-census-public-us-west-2/`). `readSOMA()`
can stream a filtered slice directly via tiledbsoma, which uses S3
byte-range requests — no full download. For programmatic Census
workflows (release version pinning, metadata discovery,
source-collection lookups), prefer the dedicated cellxgene.census R
package; `readSOMA()` converts the resulting SOMA experiment to Seurat
in a single step.

The hub dispatcher
([`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md))
recognises only file-extension formats; SOMA URIs do not have one and
therefore cannot be passed to `scConvert("s3://...", "out.h5ad")`
directly. Use the explicit two-step pattern instead:

    obj <- readSOMA("s3://...", obs_query = "...")
    writeH5AD(obj, "out.h5ad")

## See also

[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md),
[cellxgene.census](https://chanzuckerberg.github.io/cellxgene-census/r/)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read a local SOMA experiment
obj <- readSOMA("path/to/experiment")

# Stream a filtered slice from CELLxGENE Census public S3
census_uri <- paste0(
  "s3://cellxgene-census-public-us-west-2/",
  "cell-census/2024-07-01/soma/census_data/homo_sapiens"
)
obj <- readSOMA(
  uri = census_uri,
  obs_query = "cell_type == 'T cell' & tissue_general == 'blood'"
)

# Read a specific measurement
obj <- readSOMA("path/to/experiment", measurement = "ATAC")
} # }
```
