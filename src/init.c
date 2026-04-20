/*
 * init.c — R native routine registration for scConvert
 *
 * Registers .Call entry points so R can invoke compiled C code directly.
 * Built into the package shared library (scConvert.so / scConvert.dll)
 * via Makevars; NOT part of the standalone CLI binary (which uses main.c).
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Defined in sc_rcall.c */
extern SEXP C_read_h5ad(SEXP, SEXP);
extern SEXP C_read_h5ad_matrix(SEXP);
extern SEXP C_read_h5ad_obs(SEXP);
extern SEXP C_read_h5ad_obsm(SEXP);
extern SEXP C_read_h5ad_obsp(SEXP);
extern SEXP C_write_h5seurat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_write_h5ad(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_read_h5ad",        (DL_FUNC) &C_read_h5ad,        2},
    {"C_read_h5ad_matrix", (DL_FUNC) &C_read_h5ad_matrix, 1},
    {"C_read_h5ad_obs",    (DL_FUNC) &C_read_h5ad_obs,    1},
    {"C_read_h5ad_obsm",   (DL_FUNC) &C_read_h5ad_obsm,   1},
    {"C_read_h5ad_obsp",   (DL_FUNC) &C_read_h5ad_obsp,   1},
    {"C_write_h5seurat",   (DL_FUNC) &C_write_h5seurat,   7},
    {"C_write_h5ad",       (DL_FUNC) &C_write_h5ad,       8},
    {NULL, NULL, 0}
};

void R_init_scConvert(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
