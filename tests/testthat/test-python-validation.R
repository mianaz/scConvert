# Test that scConvert outputs are fully compatible with the Python single-cell ecosystem
# Uses reticulate to call scanpy, anndata, squidpy, mudata, and loompy
# Requires: conda env "scverse" with scanpy, anndata, squidpy, mudata, loompy

# Set threading env vars BEFORE Python starts to avoid OMP/numba segfaults on macOS arm64.
# These must be set at the R process level so they are inherited by the Python subprocess.
Sys.setenv(
  NUMBA_THREADING_LAYER = "tbb",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

# Skip entirely if reticulate not available or scverse env not found
skip_if_not_installed("reticulate")

scverse_path <- "/opt/homebrew/Caskroom/miniconda/base/envs/scverse"
skip_if(!dir.exists(scverse_path), "scverse conda env not found")

library(reticulate)
use_condaenv("scverse", required = FALSE)

skip_if(!reticulate::py_module_available("anndata"), "anndata not available in Python")
skip_if(!reticulate::py_module_available("scanpy"), "scanpy not available in Python")

# Set matplotlib to non-interactive backend
py_run_string("
import matplotlib
matplotlib.use('Agg')
")

# Define fix_arrow_strings once (anndata 0.11+ reads strings as ArrowStringArray
# but cannot write them back — this is an anndata bug, not scConvert's)
py_run_string("
import numpy as _np
import pandas as _pd

def _fix_arrow_strings(adata):
    '''Convert ArrowStringArrays to numpy object arrays for anndata write compat.'''
    def _fix_df(df):
        if hasattr(df.index, 'dtype') and 'arrow' in str(df.index.dtype).lower():
            df.index = _pd.Index(_np.array(df.index, dtype=object))
        for col in df.columns:
            if hasattr(df[col], 'dtype') and 'arrow' in str(df[col].dtype).lower():
                df[col] = _pd.array(_np.array(df[col], dtype=object))
    _fix_df(adata.obs)
    _fix_df(adata.var)
    if adata.raw is not None:
        try:
            raw_adata = adata.raw.to_adata()
            _fix_df(raw_adata.obs)
            _fix_df(raw_adata.var)
            adata.raw = raw_adata
        except Exception:
            adata.raw = None
    return adata
")
fix_arrow <- py$`_fix_arrow_strings`

# Helper: run Python code and return result, keeping objects as Python references
py_exec <- function(code) py_run_string(code)

# ── Setup: generate test files ──────────────────────────────────────────────

tmpdir <- tempdir()
library(Seurat)
library(scConvert)

pbmc <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 10, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)

# Direct h5ad
h5ad_direct <- file.path(tmpdir, "test_pbmc_direct.h5ad")
SeuratToH5AD(pbmc, h5ad_direct, overwrite = TRUE)

# Via h5Seurat
h5s_path <- file.path(tmpdir, "test_pbmc.h5Seurat")
scSaveH5Seurat(pbmc, h5s_path, overwrite = TRUE, verbose = FALSE)
h5ad_via <- file.path(tmpdir, "test_pbmc.h5ad")
scConvert(h5s_path, dest = "h5ad", overwrite = TRUE, verbose = FALSE)

# Loom
loom_path <- file.path(tmpdir, "test_pbmc.loom")
scSaveLoom(pbmc, loom_path, overwrite = TRUE, verbose = FALSE)

# ── Helper to read h5ad and inspect via Python ──────────────────────────────

# Use py_run_string to keep adata as a Python object, avoiding
# reticulate's auto-conversion which mangles pandas/numpy types

read_h5ad_py <- function(path, varname = "_adata") {
  py_run_string(sprintf("
import anndata as _ad
%s = _ad.read_h5ad('%s')
", varname, path))
}

# ── Phase 1: anndata reads scConvert h5ad ───────────────────────────────────

test_that("anndata reads DirectSeuratToH5AD output", {
  read_h5ad_py(h5ad_direct)
  py_exec("
_n_obs = _adata.n_obs
_n_vars = _adata.n_vars
_has_X = _adata.X is not None
_obs_cols = list(_adata.obs.columns)
_var_names = list(_adata.var.index[:5])
import numpy as _np
_X_arr = _adata.X.toarray() if hasattr(_adata.X, 'toarray') else _adata.X
_has_nan = bool(_np.isnan(_X_arr).any())
")
  expect_true(py$`_n_obs` > 0)
  expect_true(py$`_n_vars` > 0)
  expect_true(py$`_has_X`)
  expect_true(length(py$`_obs_cols`) > 0)
  expect_true(length(py$`_var_names`) > 0)
  expect_false(py$`_has_nan`)
})

test_that("obsm reductions have correct shape (direct h5ad)", {
  read_h5ad_py(h5ad_direct)
  py_exec("
_obsm_keys = list(_adata.obsm.keys())
_has_pca = 'X_pca' in _adata.obsm
_has_umap = 'X_umap' in _adata.obsm
_pca_shape = tuple(_adata.obsm['X_pca'].shape) if _has_pca else (0, 0)
_umap_shape = tuple(_adata.obsm['X_umap'].shape) if _has_umap else (0, 0)
")
  expect_true(py$`_has_pca`)
  expect_true(py$`_has_umap`)
  expect_equal(py$`_pca_shape`[[1]], py$`_n_obs`)
  expect_equal(py$`_umap_shape`[[1]], py$`_n_obs`)
  expect_equal(py$`_umap_shape`[[2]], 2L)
})

test_that("obsp graphs present (direct h5ad)", {
  read_h5ad_py(h5ad_direct)
  py_exec("_obsp_keys = list(_adata.obsp.keys())")
  expect_true(length(py$`_obsp_keys`) > 0)
})

test_that("anndata reads h5Seurat-converted h5ad", {
  read_h5ad_py(h5ad_via, "_adata2")
  py_exec("
_n_obs2 = _adata2.n_obs
_n_vars2 = _adata2.n_vars
_obs_cols2 = list(_adata2.obs.columns)
")
  expect_true(py$`_n_obs2` > 0)
  expect_true(py$`_n_vars2` > 0)
  expect_true(length(py$`_obs_cols2`) > 0)
})

# ── Phase 2: scanpy operations on scConvert output ─────────────────────────

test_that("scanpy preprocessing runs on scConvert h5ad", {
  read_h5ad_py(h5ad_direct, "_adata_pp")
  py_exec("
import scanpy as _sc
_sc.pp.filter_genes(_adata_pp, min_cells=1)
_pp_nvars = _adata_pp.n_vars
_adata_copy = _adata_pp.copy()
_sc.pp.normalize_total(_adata_copy)
_sc.pp.log1p(_adata_copy)
_pp_ok = True
")
  expect_true(py$`_pp_nvars` > 0)
  expect_true(py$`_pp_ok`)
})

test_that("scanpy neighbors + leiden run on scConvert h5ad", {
  read_h5ad_py(h5ad_direct, "_adata_cl")
  py_exec("
import scanpy as _sc
_sc.pp.neighbors(_adata_cl, use_rep='X_pca', n_neighbors=10)
_sc.tl.leiden(_adata_cl, flavor='igraph')
_has_leiden = 'leiden' in _adata_cl.obs.columns
_obsp_after = list(_adata_cl.obsp.keys())
")
  expect_true(py$`_has_leiden`)
  expect_true(length(py$`_obsp_after`) > 0)
})

test_that("scanpy rank_genes_groups runs on scConvert h5ad", {
  read_h5ad_py(h5ad_direct, "_adata_de")
  py_exec("
import scanpy as _sc
_obs_cols_de = list(_adata_de.obs.columns)
_groupby = 'seurat_clusters' if 'seurat_clusters' in _obs_cols_de else _obs_cols_de[0]
_sc.tl.rank_genes_groups(_adata_de, groupby=_groupby, method='t-test')
_has_rgg = 'rank_genes_groups' in _adata_de.uns
")
  expect_true(py$`_has_rgg`)
})

test_that("scanpy UMAP plot runs on scConvert h5ad", {
  read_h5ad_py(h5ad_direct, "_adata_pl")
  py_exec("
import scanpy as _sc
_obs_cols_pl = list(_adata_pl.obs.columns)
_sc.pl.umap(_adata_pl, color=_obs_cols_pl[0], show=False)
_plot_ok = True
")
  expect_true(py$`_plot_ok`)
})

# ── Phase 3: Python write-back roundtrip ───────────────────────────────────

test_that("anndata can write back scConvert h5ad (roundtrip)", {
  rt_path <- file.path(tmpdir, "test_pbmc_rt.h5ad")
  read_h5ad_py(h5ad_direct, "_adata_rt")
  py_exec(sprintf("
_fix_arrow_strings(_adata_rt)
_adata_rt.write_h5ad('%s')
import anndata as _ad
_adata_rt2 = _ad.read_h5ad('%s')
_rt_n_obs = _adata_rt2.n_obs
_rt_n_vars = _adata_rt2.n_vars
_rt_orig_obs = _adata_rt.n_obs
_rt_orig_vars = _adata_rt.n_vars
", rt_path, rt_path))
  expect_equal(py$`_rt_n_obs`, py$`_rt_orig_obs`)
  expect_equal(py$`_rt_n_vars`, py$`_rt_orig_vars`)
})

# ── Phase 4: loompy reads scConvert loom ───────────────────────────────────

loompy_available <- reticulate::py_module_available("loompy")

test_that("loompy reads scConvert loom file", {
  skip_if(!loompy_available, "loompy not available")

  py_exec(sprintf("
import loompy
_ds = loompy.connect('%s', mode='r')
_loom_shape = _ds.shape
_ra_keys = list(_ds.ra.keys())
_ca_keys = list(_ds.ca.keys())
_has_gene = 'Gene' in _ds.ra
_has_cellid = 'CellID' in _ds.ca
_slice = _ds[0:5, 0:5]
_slice_shape = _slice.shape
_ds.close()
", loom_path))

  expect_true(py$`_loom_shape`[[1]] > 0)
  expect_true(py$`_loom_shape`[[2]] > 0)
  expect_true(length(py$`_ra_keys`) > 0)
  expect_true(length(py$`_ca_keys`) > 0)
  expect_true(py$`_has_gene`)
  expect_true(py$`_has_cellid`)
  expect_equal(py$`_slice_shape`[[1]], 5L)
  expect_equal(py$`_slice_shape`[[2]], 5L)
})

test_that("scanpy reads scConvert loom file", {
  skip_if(!loompy_available, "loompy not available")

  py_exec(sprintf("
import scanpy as _sc
_adata_loom = _sc.read_loom('%s', sparse=True)
_loom_nobs = _adata_loom.n_obs
_loom_nvars = _adata_loom.n_vars
", loom_path))
  expect_true(py$`_loom_nobs` > 0)
  expect_true(py$`_loom_nvars` > 0)
})

# ── Phase 5: Python-native h5ad -> scConvert -> Python roundtrip ───────────

test_that("Python->R->Python roundtrip works", {
  native_path <- file.path(tmpdir, "scanpy_native.h5ad")
  py_exec(sprintf("
import anndata as _ad
import scanpy as _sc
import numpy as _np
import pandas as _pd
from scipy.sparse import random as _sparse_random

_np.random.seed(42)
n_obs, n_vars = 200, 500
X = _sparse_random(n_obs, n_vars, density=0.1, format='csr', dtype=_np.float32)
obs = _pd.DataFrame({
    'cell_type': _pd.Categorical(_np.random.choice(['T', 'B', 'Mono', 'NK'], n_obs)),
    'n_counts': _np.random.randint(500, 5000, n_obs).astype(_np.float64),
}, index=[f'cell_{i}' for i in range(n_obs)])
var = _pd.DataFrame({
    'highly_variable': _np.random.choice([True, False], n_vars),
}, index=[f'gene_{i}' for i in range(n_vars)])
_adata_syn = _ad.AnnData(X=X, obs=obs, var=var)
_sc.pp.normalize_total(_adata_syn)
_sc.pp.log1p(_adata_syn)
_sc.pp.pca(_adata_syn, n_comps=20)
_sc.pp.neighbors(_adata_syn)
_sc.tl.umap(_adata_syn)
_adata_syn.write_h5ad('%s')
_syn_nobs = _adata_syn.n_obs
_syn_nvars = _adata_syn.n_vars
", native_path))

  expect_true(file.exists(native_path))

  # Load in scConvert
  obj <- LoadH5AD(native_path, verbose = FALSE)
  expect_equal(ncol(obj), 200L)
  expect_equal(nrow(obj), 500L)
  expect_true(length(Reductions(obj)) > 0)

  # Write back to h5ad
  rt_path <- file.path(tmpdir, "scanpy_roundtrip.h5ad")
  SeuratToH5AD(obj, rt_path, overwrite = TRUE)

  # Read back in Python and run scanpy
  py_exec(sprintf("
import anndata as _ad
import scanpy as _sc
_adata_pyrt = _ad.read_h5ad('%s')
_pyrt_nobs = _adata_pyrt.n_obs
_pyrt_nvars = _adata_pyrt.n_vars
_sc.pp.neighbors(_adata_pyrt, use_rep='X_pca', n_neighbors=10)
_sc.tl.leiden(_adata_pyrt, flavor='igraph')
_pyrt_has_leiden = 'leiden' in _adata_pyrt.obs.columns
", rt_path))
  expect_equal(py$`_pyrt_nobs`, 200L)
  expect_equal(py$`_pyrt_nvars`, 500L)
  expect_true(py$`_pyrt_has_leiden`)
})

# ── Phase 6: Spatial validation with squidpy ───────────────────────────────

has_stxBrain <- tryCatch({
  requireNamespace("SeuratData", quietly = TRUE) &&
    "stxBrain.SeuratData" %in% rownames(SeuratData::InstalledData()) &&
    SeuratData::InstalledData()["stxBrain.SeuratData", "Installed"]
}, error = function(e) FALSE)

squidpy_available <- reticulate::py_module_available("squidpy")

if (has_stxBrain) {
  library(SeuratData)
  brain <- UpdateSeuratObject(LoadData("stxBrain", type = "anterior1"))
  brain <- NormalizeData(brain, verbose = FALSE)
  spatial_h5ad <- file.path(tmpdir, "test_brain_spatial.h5ad")
  SeuratToH5AD(brain, spatial_h5ad, overwrite = TRUE)
}

test_that("anndata reads spatial h5ad with correct obsm/spatial shape", {
  skip_if(!has_stxBrain, "stxBrain not installed")

  read_h5ad_py(spatial_h5ad, "_adata_sp")
  py_exec("
import numpy as _np
_sp_nobs = _adata_sp.n_obs
_has_spatial = 'spatial' in _adata_sp.obsm
if _has_spatial:
    _sp_coords = _adata_sp.obsm['spatial']
    _sp_coords_shape = tuple(_sp_coords.shape)
    _sp_coords_nan = bool(_np.isnan(_sp_coords).any())
else:
    _sp_coords_shape = (0, 0)
    _sp_coords_nan = True
")
  expect_true(py$`_sp_nobs` > 0)
  expect_true(py$`_has_spatial`)
  expect_equal(py$`_sp_coords_shape`[[1]], py$`_sp_nobs`)
  expect_true(py$`_sp_coords_shape`[[2]] >= 2)
  expect_false(py$`_sp_coords_nan`)
})

test_that("uns/spatial has images and scale factors", {
  skip_if(!has_stxBrain, "stxBrain not installed")

  read_h5ad_py(spatial_h5ad, "_adata_sp2")
  py_exec("
_has_uns_spatial = 'spatial' in _adata_sp2.uns
_lib_ids = []
_has_images = False
_has_sf = False
_img_shape = ()
if _has_uns_spatial:
    _lib_ids = list(_adata_sp2.uns['spatial'].keys())
    if len(_lib_ids) > 0:
        lib = _adata_sp2.uns['spatial'][_lib_ids[0]]
        _has_images = 'images' in lib
        _has_sf = 'scalefactors' in lib
        if _has_images:
            _img_keys = list(lib['images'].keys())
            if len(_img_keys) > 0:
                _img_shape = tuple(lib['images'][_img_keys[0]].shape)
")
  expect_true(py$`_has_uns_spatial`)
  expect_true(length(py$`_lib_ids`) > 0)
  expect_true(py$`_has_images`)
  expect_true(py$`_has_sf`)
  expect_true(length(py$`_img_shape`) == 3)  # H x W x 3
})

test_that("squidpy spatial_neighbors runs on scConvert spatial h5ad", {
  skip_if(!has_stxBrain, "stxBrain not installed")
  skip_if(!squidpy_available, "squidpy not available")

  read_h5ad_py(spatial_h5ad, "_adata_sq")
  py_exec("
import squidpy as _sq
_sq.gr.spatial_neighbors(_adata_sq, coord_type='generic')
_sq_obsp = list(_adata_sq.obsp.keys())
_sq_has_conn = len(_sq_obsp) > 0
")
  expect_true(py$`_sq_has_conn`)
})

test_that("squidpy spatial_autocorr (Moran's I) runs on scConvert spatial h5ad", {
  skip_if(!has_stxBrain, "stxBrain not installed")
  skip_if(!squidpy_available, "squidpy not available")

  read_h5ad_py(spatial_h5ad, "_adata_sq2")
  py_exec("
import squidpy as _sq
_sq.gr.spatial_neighbors(_adata_sq2, coord_type='generic')
_genes = list(_adata_sq2.var_names[:5])
_sq.gr.spatial_autocorr(_adata_sq2, mode='moran', genes=_genes)
_has_moran = 'moranI' in _adata_sq2.uns
")
  expect_true(py$`_has_moran`)
})

# ── Phase 7: h5mu validation with mudata ───────────────────────────────────

has_cbmc <- tryCatch({
  requireNamespace("SeuratData", quietly = TRUE) &&
    "cbmc.SeuratData" %in% rownames(SeuratData::InstalledData())
}, error = function(e) FALSE)

mudata_available <- reticulate::py_module_available("mudata")

if (has_cbmc) {
  library(SeuratData)
  data("cbmc", package = "cbmc.SeuratData")
  cbmc <- UpdateSeuratObject(cbmc)
  overlap <- intersect(rownames(cbmc[["ADT"]]), rownames(cbmc[["RNA"]]))
  if (length(overlap) > 0) {
    adt_keep <- setdiff(rownames(cbmc[["ADT"]]), overlap)
    cbmc[["ADT"]] <- subset(cbmc[["ADT"]], features = adt_keep)
  }
  h5mu_path <- file.path(tmpdir, "test_cbmc.h5mu")
  SaveH5MU(cbmc, h5mu_path, overwrite = TRUE)
}

test_that("mudata reads scConvert h5mu with both modalities", {
  skip_if(!has_cbmc, "cbmc not installed")
  skip_if(!mudata_available, "mudata not available")

  py_exec(sprintf("
import mudata as _md
_mdata = _md.read_h5mu('%s')
_mu_nobs = _mdata.n_obs
_mod_keys = list(_mdata.mod.keys())
_has_rna = 'rna' in _mdata.mod
_has_prot = 'prot' in _mdata.mod
_rna_shape = tuple(_mdata.mod['rna'].shape) if _has_rna else (0, 0)
_prot_shape = tuple(_mdata.mod['prot'].shape) if _has_prot else (0, 0)
_rna_has_X = _mdata.mod['rna'].X is not None if _has_rna else False
", h5mu_path))

  expect_true(py$`_mu_nobs` > 0)
  expect_true(length(py$`_mod_keys`) >= 2)
  expect_true(py$`_has_rna`)
  expect_true(py$`_has_prot`)
  expect_true(py$`_rna_shape`[[1]] > 0 && py$`_rna_shape`[[2]] > 0)
  expect_true(py$`_prot_shape`[[1]] > 0 && py$`_prot_shape`[[2]] > 0)
  expect_true(py$`_rna_has_X`)
})

test_that("scanpy preprocessing runs on h5mu RNA modality", {
  skip_if(!has_cbmc, "cbmc not installed")
  skip_if(!mudata_available, "mudata not available")

  py_exec(sprintf("
import mudata as _md
import scanpy as _sc
_mdata2 = _md.read_h5mu('%s')
_rna2 = _mdata2.mod['rna'].copy()
_sc.pp.normalize_total(_rna2)
_sc.pp.log1p(_rna2)
_mu_pp_ok = True
", h5mu_path))
  expect_true(py$`_mu_pp_ok`)
})

# ── Cleanup ────────────────────────────────────────────────────────────────

test_that("cleanup temp files", {
  files <- list.files(tmpdir, pattern = "^test_|^scanpy_", full.names = TRUE)
  unlink(files)
  expect_true(TRUE)
})
