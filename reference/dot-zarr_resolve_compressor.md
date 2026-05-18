# Resolve a writer-side compressor spec to a zarr-codec list

Accepts a user-supplied compressor argument and returns a list
`list(id, level, ...)` suitable for
[`.zarr_compress()`](https://mianaz.github.io/scConvert/reference/dot-zarr_compress.md)
and for writing into a `.zarray` `compressor` field. `NULL` means
"auto": Zstd when a Zstd codec package (currently zstdlite) is on the
search path; otherwise zlib at the current package compression level. As
of 2026 there is no maintained CRAN Zstd-bytes wrapper, so the auto path
almost always lands on zlib; the resolver still picks Zstd if the user
has installed a provider out of band.

## Usage

``` r
.zarr_resolve_compressor(x = NULL, level = NULL)
```

## Arguments

- x:

  One of `NULL`, a string (`"zstd"`, `"gzip"`, `"zlib"`, `"blosc"`,
  `"none"`), or a list with at least an `id` field.

- level:

  Override the codec level. When `NULL`, a per-codec default is used.

## Value

Either `NULL` (no compression) or a list with `id` and typically
`level`.
