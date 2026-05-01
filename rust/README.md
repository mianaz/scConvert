# `rust/` — scConvert native core (B1 scaffolding)

This directory hosts the Rust migration of the scConvert conversion engine,
as described in Part B of
`/Users/miana/.claude/plans/codex-gave-very-critical-structured-knuth.md`.

**Current status: B1 (contracts only).** No conversion code lives here yet.
The existing C CLI at `src/scconvert` remains the production path. This
workspace exists so the types, error codes, fidelity report schema, route
planner schema, and C ABI can be reviewed and stabilised before any IO
work lands in B2.

## Layout

```
rust/
  Cargo.toml                # workspace manifest
  scconvert-core/           # pure Rust library (data model, schemas,
                            #   route planner, format-agnostic logic)
  scconvert-capi/           # stable C ABI over core; cbindgen-generated
                            #   header `scconvert.h`
  scconvert-cli/            # thin CLI binary that links scconvert-core
  scconvert.h               # (generated during build, committed for R)
```

## Provisional dependency decisions

These defaults come from the Part B plan. Change here; revisit before B2
if something looks wrong.

| Concern | Default | Alternative | Notes |
|---|---|---|---|
| HDF5 binding | `hdf5-metno` (active fork of `hdf5`, vendors libhdf5 via `hdf5-src`) | `hdf5` (original, unmaintained) | Picked for reproducible builds and ongoing maintenance. Not yet added — lands in B2. |
| Zarr backend | `zarrs` | hand-rolled | `zarrs` has Zarr v2+v3 and the only non-Python v3 implementation. Lands in B4. |
| C ABI generation | `cbindgen` | hand-written header | `cbindgen` keeps header and Rust in sync. |
| R ↔ Rust bridge | `.Call` into the C ABI | `extendr` | `.Call` keeps R package builds independent of a Rust toolchain on user machines. |
| Python ↔ Rust | `pyo3` + `maturin` | hand-written CPython | Lands in B7. |
| C CLI retirement | Remove after B4 ships Zarr | Keep forever | ~11k lines of C replaced by the Rust workspace. |

## Build

```bash
cd rust
cargo check                 # validate types and schemas (fast, no codegen)
cargo test -p scconvert-core
cargo build --release       # produces cli binary + capi shared lib
```

The R package never invokes `cargo` during install. If a user wants the
Rust backend they build it explicitly here, same pattern as the C CLI.

## What "done" looks like for B1

- [x] Workspace compiles with `cargo check`.
- [x] `scconvert-core` exports `Format`, `DType`, `SparseLayout`,
      `FidelityReport`, `Route`, `ConvertError`.
- [x] `scconvert-core` has a failing-stub `convert()` returning
      `ConvertError::NotYetImplemented`.
- [x] JSON schemas for `FidelityReport` and `Route` are serde-roundtrippable.
- [x] `scconvert-capi` defines the C entry points (`sc_convert`,
      `sc_plan`, `sc_inspect`, `sc_free`) that all delegate to core and
      currently return error code 1 with `NotYetImplemented`.
- [x] `scconvert-cli` parses `convert`, `plan`, `inspect` subcommands and
      calls core.
- [x] A conformance-fixture manifest (JSON) points at the existing
      `inst/testdata/` files.

B2 lifts `convert()` from stub to real h5ad ↔ h5seurat streaming.
