"""Generate a feature-rich AnnData fixture with whatever anndata version is on the path.

Usage: python make_fixture.py OUTDIR
Writes OUTDIR/fixture_<version>.h5ad (and .zarr when zarr is importable).
"""
import sys, os, warnings, json
warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
# --infer-string on|off controls pandas 3 string inference (StringArray columns/index)
_mode = sys.argv[2] if len(sys.argv) > 2 else "auto"
if _mode in ("on", "off") and hasattr(pd.options, "future") and hasattr(pd.options.future, "infer_string"):
    pd.options.future.infer_string = (_mode == "on")
import scipy.sparse as sp
import anndata as ad

try:
    ver = ad.__version__
except Exception:
    from importlib.metadata import version
    ver = version("anndata")
tag = "ad" + ver.replace(".", "_") + ("_strings" if _mode == "on" else "")
outdir = sys.argv[1]
os.makedirs(outdir, exist_ok=True)

rng = np.random.default_rng(0)
n_obs, n_var = int(os.environ.get("FIX_NOBS", "30")), int(os.environ.get("FIX_NVAR", "12"))
counts = rng.poisson(1.0, size=(n_obs, n_var)).astype(np.float32)
X = sp.csr_matrix(np.log1p(counts))
obs_names = [f"cell_{i}" for i in range(n_obs)]
var_names = [f"gene_{i}" for i in range(n_var)]
var_names[3] = "gene_3_with_underscore"

obs = pd.DataFrame(index=obs_names)
obs["cluster"] = pd.Categorical(rng.choice(["b", "a", "c"], n_obs), categories=["b", "a", "c"])
obs["grade"] = pd.Categorical(rng.choice(["low", "high", "mid"], n_obs), categories=["low", "mid", "high"], ordered=True)
obs["n_counts"] = counts.sum(1).astype(np.float64)
obs["n_genes"] = (counts > 0).sum(1).astype(np.int64)
obs["is_doublet"] = rng.random(n_obs) > 0.8
obs["sample"] = ["s1" if i % 2 else "s2" for i in range(n_obs)]  # plain string/object column
major = tuple(int(x) for x in ver.split(".")[:2])
if major >= (0, 8):
    nullable_int = pd.array(rng.integers(0, 5, n_obs), dtype="Int64")
    nullable_int[0] = pd.NA
    obs["batch_id"] = nullable_int
    nullable_bool = pd.array(rng.random(n_obs) > 0.5, dtype="boolean")
    nullable_bool[1] = pd.NA
    obs["qc_pass"] = nullable_bool
if major >= (0, 11) and _mode == "on":
    try:
        ad.settings.allow_write_nullable_strings = True
    except Exception:
        pass
    strs = pd.array(["x" if i % 3 else "y" for i in range(n_obs)], dtype="string")
    strs[2] = pd.NA
    obs["note"] = strs  # pandas StringArray -> nullable-string-array on new anndata

var = pd.DataFrame(index=var_names)
var["highly_variable"] = rng.random(n_var) > 0.5
var["gene_type"] = pd.Categorical(rng.choice(["coding", "lnc"], n_var))
var["mean"] = counts.mean(0).astype(np.float64)

obsm = {"X_pca": rng.normal(size=(n_obs, 5)).astype(np.float32), "X_umap": rng.normal(size=(n_obs, 2))}
varm = {"PCs": rng.normal(size=(n_var, 5))}
conn = sp.random(n_obs, n_obs, density=0.2, format="csr", random_state=1)
obsp = {"connectivities": conn, "distances": sp.csc_matrix(conn)}
varp = {"gene_corr": sp.random(n_var, n_var, density=0.3, format="csr", random_state=2)}
uns = {
    "title": "fixture",
    "n_pcs": 5,
    "ratio": 0.5,
    "flag": True,
    "arr": np.arange(4).astype(np.int64),
    "nested": {"a": 1, "b": "two", "c": [1.0, 2.0, 3.0]},
    "cluster_colors": np.array(["#ff0000", "#00ff00", "#0000ff"]),
    "df": pd.DataFrame({"x": [1, 2], "y": ["p", "q"]}, index=["r1", "r2"]),
}
if major >= (0, 9):
    uns["none_value"] = None
layers = {"counts": sp.csr_matrix(counts), "dense_counts": counts}

adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm, varm=varm, obsp=obsp, varp=varp, uns=uns, layers=layers)
adata.raw = ad.AnnData(X=sp.csr_matrix(counts), var=var[["gene_type"]].copy())
adata.raw.var.index = var_names

h5 = os.path.join(outdir, f"fixture_{tag}.h5ad")
adata.write_h5ad(h5, compression="gzip")
print("wrote", h5)
try:
    import zarr
    z = os.path.join(outdir, f"fixture_{tag}.zarr")
    adata.write_zarr(z)
    print("wrote", z)
except Exception as e:
    print("zarr skipped:", type(e).__name__, e)

# manifest of what we expect, for R-side assertions
manifest = {
    "version": ver,
    "obs_columns": list(obs.columns),
    "cluster_levels": ["b", "a", "c"],
    "grade_levels": ["low", "mid", "high"],
    "n_obs": n_obs, "n_var": n_var,
    "var_names": var_names,
}
with open(os.path.join(outdir, f"fixture_{tag}.json"), "w") as fh:
    json.dump(manifest, fh)
