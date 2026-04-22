# Run two-step streaming conversion via a temp h5seurat file

Creates a temp h5seurat file, streams source -\> h5seurat, then h5seurat
-\> dest. The temp file is deleted on exit.

## Usage

``` r
.stream_composite_via_h5seurat(
  source,
  dest,
  step1_fn,
  step2_fn,
  step1_args = list(),
  step2_args = list(),
  gzip = 4L,
  verbose = TRUE
)
```

## Arguments

- source:

  Input file path

- dest:

  Output file path

- step1_fn:

  Function to convert source -\> h5seurat (temp)

- step2_fn:

  Function to convert h5seurat (temp) -\> dest

- step1_args:

  Named list of additional args for step1_fn

- step2_args:

  Named list of additional args for step2_fn

- gzip:

  Integer gzip compression level

- verbose:

  Logical
