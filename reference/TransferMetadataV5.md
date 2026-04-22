# Transfer metadata from V5 h5Seurat to h5ad

Transfer metadata from V5 h5Seurat to h5ad

## Usage

``` r
TransferMetadataV5(source, dfile, dname = "obs", index = NULL, verbose = TRUE)
```

## Arguments

- source:

  Source H5 group (potentially wrapped in environment)

- dfile:

  Destination H5 file

- dname:

  Destination group name

- index:

  Cell index for subsetting

- verbose:

  Logical, print messages

## Value

NULL invisibly
