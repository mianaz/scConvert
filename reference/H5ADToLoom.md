# Convert an h5ad file to loom format

Reads the h5ad into a Seurat object and writes it as loom.

## Usage

``` r
H5ADToLoom(source, dest, overwrite = FALSE, gzip = 4L, verbose = TRUE)
```

## Arguments

- source:

  Path to input .h5ad file

- dest:

  Path for output .loom file

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
