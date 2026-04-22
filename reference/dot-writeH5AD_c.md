# Fast C-based h5ad writer

Writes a Seurat object directly to h5ad format using native C HDF5
routines. Exploits the zero-copy CSC\<-\>CSR reinterpretation:
dgCMatrix(genesxcells) CSC is identical to h5ad's cellsxgenes CSR with
relabeled arrays.

## Usage

``` r
.writeH5AD_c(object, filename, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- object:

  Seurat object

- filename:

  Output file path

- overwrite:

  Allow overwriting

- verbose:

  Show progress

## Value

TRUE on success, FALSE on failure
