# Read obsp graphs from a SOMA measurement

Read obsp graphs from a SOMA measurement

## Usage

``` r
.soma_read_obsp(ms, obj, measurement, cell_names, verbose = TRUE)
```

## Arguments

- ms:

  A SOMAMeasurement object

- obj:

  A Seurat object to add graphs to

- measurement:

  Assay name

- cell_names:

  Cell names for row/column names

- verbose:

  Verbosity

## Value

The Seurat object with graphs added
