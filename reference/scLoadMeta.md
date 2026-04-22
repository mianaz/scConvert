# Load deferred h5ad components into a Seurat object

When `readH5AD` was called with a subset of `components`, this function
loads the remaining components from the stored file path.

## Usage

``` r
scLoadMeta(object, components = NULL, verbose = TRUE)
```

## Arguments

- object:

  A Seurat object created by `readH5AD` with partial components

- components:

  Character vector of additional components to load, or NULL for all

- verbose:

  Show progress messages

## Value

The Seurat object with additional components loaded
