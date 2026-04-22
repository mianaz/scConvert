# Fast C-based h5Seurat writer

Writes a Seurat object to h5Seurat format using native C HDF5 routines,
bypassing hdf5r's R6 method overhead (~10x faster).

## Usage

``` r
.writeH5Seurat_c(object, filename, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- object:

  Seurat object

- filename:

  Output file path

- overwrite:

  Allow overwriting existing file

- verbose:

  Show progress

## Value

TRUE on success, FALSE on failure
