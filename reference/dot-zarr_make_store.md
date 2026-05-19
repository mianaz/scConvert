# Storage abstraction for zarr v2 stores

Returns a list of method closures that hide the difference between a
local-filesystem store and a remote HTTP-backed store. Every `.zarr_*`
read helper that previously took a `store_path` string now also accepts
the list returned here.

## Usage

``` r
.zarr_make_store(x, cache = NA)
```

## Arguments

- x:

  Either a filesystem path, an `s3://`/`gs://`/HTTP URL, or an
  already-constructed store list.

- cache:

  For HTTP stores: `TRUE` to persist fetched objects under
  `tools::R_user_dir("scConvert", "cache")`, `FALSE` for a tempdir that
  is discarded with the session, `NA` for default (TRUE).

## Value

A store list.

## Details

Methods:

- `kind`:

  `"local"` or `"http"`.

- `root`:

  The store root URI (filesystem path or URL).

- `exists(rel)`:

  Logical: does `rel` exist in the store?

- `read_bytes(rel)`:

  Raw vector with the contents of `rel`. Errors on missing keys; callers
  should `store$exists(rel)` first if a missing key is expected.

- `read_json(rel)`:

  Parsed JSON list, or empty list if missing.

- `list_dirs(rel)`:

  Character vector of immediate subgroup names under `rel`
  (zarr-meaningful only; filters `c/`, dotfiles, and `__*`).
