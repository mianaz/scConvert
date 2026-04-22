# Stream an h5ad-style DataFrame group to h5seurat meta.data

Reads obs columns (categorical or plain) from an AnnData DataFrame group
and writes them to h5seurat meta.data format (factors use
levels/values).

## Usage

``` r
.stream_h5ad_df_to_h5seurat_md(src_group, dst_parent, gzip, verbose = TRUE)
```

## Arguments

- src_group:

  H5Group for the source DataFrame (obs)

- dst_parent:

  H5Group or H5File to write meta.data into

- gzip:

  Integer compression level

- verbose:

  Logical

## Value

Character vector of column names written
