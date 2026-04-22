# Stream h5mu to h5seurat (direct HDF5-to-HDF5)

For each modality in /mod/, converts X, obs, var, obsm, obsp to h5seurat
layout. Uses the zero-copy CSR\<-\>CSC reinterpretation trick.

## Usage

``` r
.stream_h5mu_to_h5seurat(h5mu_path, h5seurat_path, gzip = 4L, verbose = TRUE)
```

## Arguments

- h5mu_path:

  Path to input h5mu file

- h5seurat_path:

  Path to output h5seurat file

- gzip:

  Integer gzip compression level

- verbose:

  Logical
