"""compare_h5ad.py A.h5ad B.h5ad  -> reports differences between two h5ad files as read by this anndata."""
import sys, warnings
warnings.simplefilter("ignore")
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

a = ad.read_h5ad(sys.argv[1]); b = ad.read_h5ad(sys.argv[2])
problems = []
def mat_eq(x, y):
    if sp.issparse(x) != sp.issparse(y): return False
    if sp.issparse(x): return (x != y).nnz == 0 and x.shape == y.shape
    return np.array_equal(np.asarray(x), np.asarray(y))
def df_eq(name, x, y):
    if list(x.index) != list(y.index): problems.append(f"{name}: index differs")
    if list(x.columns) != list(y.columns): problems.append(f"{name}: columns {list(x.columns)} vs {list(y.columns)}"); return
    for c in x.columns:
        xc, yc = x[c], y[c]
        if str(xc.dtype) != str(yc.dtype): problems.append(f"{name}.{c}: dtype {xc.dtype} vs {yc.dtype}")
        if isinstance(xc.dtype, pd.CategoricalDtype) and isinstance(yc.dtype, pd.CategoricalDtype):
            if list(xc.cat.categories) != list(yc.cat.categories): problems.append(f"{name}.{c}: categories differ")
            if xc.cat.ordered != yc.cat.ordered: problems.append(f"{name}.{c}: ordered {xc.cat.ordered} vs {yc.cat.ordered}")
        try:
            same = xc.astype(object).where(~xc.isna(), None).tolist() == yc.astype(object).where(~yc.isna(), None).tolist()
        except Exception as e:
            same = False
        if not same: problems.append(f"{name}.{c}: values differ")
if a.shape != b.shape: problems.append(f"shape {a.shape} vs {b.shape}")
if not mat_eq(a.X, b.X): problems.append("X differs")
df_eq("obs", a.obs, b.obs); df_eq("var", a.var, b.var)
for k in set(a.layers) | set(b.layers):
    if k not in a.layers or k not in b.layers or not mat_eq(a.layers[k], b.layers[k]): problems.append(f"layers.{k} differs")
for slot in ["obsm", "varm", "obsp", "varp"]:
    xa, xb = getattr(a, slot), getattr(b, slot)
    for k in set(xa) | set(xb):
        if k not in xa or k not in xb or not mat_eq(xa[k], xb[k]): problems.append(f"{slot}.{k} differs")
def uns_eq(x, y, path):
    if isinstance(x, dict) and isinstance(y, dict):
        for k in set(x) | set(y):
            if k not in x or k not in y: problems.append(f"uns{path}.{k} missing on one side"); continue
            uns_eq(x[k], y[k], path + "." + k)
    elif isinstance(x, pd.DataFrame) and isinstance(y, pd.DataFrame):
        df_eq("uns" + path, x, y)
    elif isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
        if not np.array_equal(np.asarray(x), np.asarray(y)): problems.append(f"uns{path} array differs")
    else:
        if not (x == y or (x is None and y is None)): problems.append(f"uns{path}: {x!r} vs {y!r}")
uns_eq(dict(a.uns), dict(b.uns), "")
if (a.raw is None) != (b.raw is None): problems.append("raw presence differs")
elif a.raw is not None:
    if not mat_eq(a.raw.X, b.raw.X): problems.append("raw.X differs")
    df_eq("raw.var", a.raw.var, b.raw.var)
print(f"anndata {ad.__version__}: {sys.argv[1]} vs {sys.argv[2]}: " + ("IDENTICAL" if not problems else "DIFFERENCES:"))
for p in problems: print("   -", p)
