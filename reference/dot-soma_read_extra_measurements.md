# Read additional measurements as extra assays

Read additional measurements as extra assays

## Usage

``` r
.soma_read_extra_measurements(
  exp,
  obj,
  primary,
  n_cells,
  cell_names,
  verbose = TRUE
)
```

## Arguments

- exp:

  A SOMAExperiment object

- obj:

  A Seurat object to add assays to

- primary:

  Name of the primary measurement (already loaded)

- n_cells:

  Number of cells in the primary measurement

- cell_names:

  Cell names

- verbose:

  Verbosity

## Value

The Seurat object with extra assays added
