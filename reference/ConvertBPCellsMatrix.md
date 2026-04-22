# Convert BPCells matrix objects to dgCMatrix

Convert BPCells matrix objects to dgCMatrix

## Usage

``` r
ConvertBPCellsMatrix(mat, verbose = FALSE)
```

## Arguments

- mat:

  Matrix object (potentially BPCells IterableMatrix or RenameDims)

- verbose:

  Show conversion message

## Value

dgCMatrix or original matrix if not BPCells
