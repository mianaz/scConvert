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

## See also

[`writeSOMA`](https://mianaz.github.io/scConvert/reference/writeSOMA.md),
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read a local SOMA experiment
obj <- readSOMA("path/to/experiment")

# Read from CELLxGENE Census with filtering
obj <- readSOMA(
  uri = "s3://cellxgene-census/soma/...",
  obs_query = "cell_type == 'T cell' & tissue == 'blood'"
)

# Read a specific measurement
obj <- readSOMA("path/to/experiment", measurement = "ATAC")
} # }
```
