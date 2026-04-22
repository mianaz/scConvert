# Stream h5seurat meta.data to an h5ad-style DataFrame group

Reads meta.data columns (factors or plain datasets) and writes them as
AnnData DataFrame format (categoricals use categories/codes).

## Usage

``` r
.stream_h5seurat_md_to_h5ad_df(
  md_group,
  dst_parent,
  dst_name,
  cell_names,
  gzip,
  verbose = TRUE
)
```

## Arguments

- md_group:

  H5Group for meta.data

- dst_parent:

  H5Group or H5File to write into

- dst_name:

  Name for the DataFrame group (e.g., "obs")

- cell_names:

  Character vector of cell barcodes for \_index

- gzip:

  Integer compression level

- verbose:

  Logical

## Value

Invisible NULL
