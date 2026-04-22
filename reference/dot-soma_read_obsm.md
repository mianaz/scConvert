# Read obsm embeddings from a SOMA measurement

Read obsm embeddings from a SOMA measurement

## Usage

``` r
.soma_read_obsm(ms, obj, measurement, cell_names, verbose = TRUE)
```

## Arguments

- ms:

  A SOMAMeasurement object

- obj:

  A Seurat object to add reductions to

- measurement:

  Assay name

- cell_names:

  Cell names for rownames

- verbose:

  Verbosity

## Value

The Seurat object with embeddings added
