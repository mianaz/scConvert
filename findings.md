scConvert review findings
========================

1. `readH5MU()` drops neighbor graphs. Inside `R/LoadH5MU.R:381-395` graphs are saved via `seurat_obj[[full_name]] <- graph_obj`; the Seurat setter for `[[<-` only handles assays/dimreducs, so the object does not reach `@graphs` (and may overwrite an assay or error). Graphs should be added through `seurat_obj@graphs[[full_name]] <- graph_obj`, matching the h5ad/zarr loaders.

2. Per-feature metadata from each modality is discarded. At `R/LoadH5MU.R:424-433` the code tries to assign a whole `var` data frame with `seurat_obj[[assay_name]][[names(var_cols)]] <- var_df`, but `[[` only accepts a single column name, so the call errors and is silently dropped. Iterate over columns (or assign to the assay’s `meta.features` slot) to preserve feature annotations.

3. Python fallbacks for h5ad conversion never verify interpreter availability. The helper `system2("python3", "-c import h5py…")` only raises an error when the binary is missing; if it exists but the import fails, the non-zero status is ignored, `python_cmd` is still returned, and each conversion attempt emits “Python failed…” warnings (`R/Convert.R:891-907`, `R/Convert.R:1047-1066`). Check the status/result before accepting the interpreter so the code can immediately fall back to the pure-R path when h5py/pandas are absent.
