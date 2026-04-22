# Run the scconvert CLI for file-to-file conversion

Performs file-to-file format conversion with three tiers:

1.  **C binary**: HDF5 format pairs (h5ad, h5seurat, h5mu)

2.  **Streaming**: zarr pairs (h5ad/h5seurat \\\leftrightarrow\\ zarr) –
    copies fields directly without creating a Seurat intermediate

3.  **R hub**: all other format pairs via
    [`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert.md)

## Usage

``` r
scConvert_cli(
  input,
  output,
  assay = "RNA",
  gzip = 4L,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- input:

  Input file path

- output:

  Output file path

- assay:

  Assay name (passed as –assay)

- gzip:

  Gzip compression level (0-9, passed as –gzip)

- overwrite:

  If TRUE, remove output file before conversion

- verbose:

  If TRUE, print CLI output

## Value

TRUE on success, FALSE on failure

## Details

Supported formats: h5ad, h5seurat, h5mu, loom, rds, zarr.
