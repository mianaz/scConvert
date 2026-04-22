# Register a direct conversion path between two file formats

Direct paths bypass the Seurat hub for efficiency (e.g., h5ad \<-\>
h5seurat can be converted via direct HDF5 operations without loading
into memory).

## Usage

``` r
RegisterDirectPath(stype, dtype, fn)
```

## Arguments

- stype:

  Source format extension string

- dtype:

  Destination format extension string

- fn:

  Conversion function(source, dest, assay, overwrite, verbose, ...)
