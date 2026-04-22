# Register a file format with its load and save functions

Register a file format with its load and save functions

## Usage

``` r
RegisterFormat(ext, loader = NULL, saver = NULL)
```

## Arguments

- ext:

  Lowercase file extension (e.g., "h5seurat", "h5ad", "loom", "h5mu",
  "rds")

- loader:

  Function(file, assay, verbose, ...) -\> Seurat object

- saver:

  Function(object, dest, overwrite, verbose, ...) -\> invisible(dest)
