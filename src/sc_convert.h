/*
 * sc_convert.h — Core header for scConvert C CLI
 *
 * Streaming on-disk conversion between h5ad, h5seurat, and h5mu formats.
 * No full dataset materialization in RAM — reads and writes in chunks.
 */

#ifndef SC_CONVERT_H
#define SC_CONVERT_H

#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

/* ── Constants ──────────────────────────────────────────────────────────────── */

#define SC_CHUNK_SIZE       (1 << 20)   /* 1M elements per streaming chunk */
#define SC_GZIP_LEVEL       1           /* level 1: 2-3x faster than 4, ~10% larger */
#define SC_MAX_NAME_LEN     256
#define SC_MAX_COLUMNS      1024
#define SC_MAX_MODALITIES   32
#define SC_VERSION_STRING   "1.0.0"

/* ── Error handling ─────────────────────────────────────────────────────────── */

#define SC_OK       0
#define SC_ERR     -1
#define SC_ERR_IO  -2
#define SC_ERR_HDF -3
#define SC_ERR_ARG -4

#define SC_CHECK_HDF(expr, msg) do { \
    if ((expr) < 0) { \
        fprintf(stderr, "HDF5 error: %s\n", msg); \
        return SC_ERR_HDF; \
    } \
} while(0)

#define SC_CHECK_HDF_GOTO(expr, msg, label) do { \
    if ((expr) < 0) { \
        fprintf(stderr, "HDF5 error: %s\n", msg); \
        goto label; \
    } \
} while(0)

/* ── Data structures ────────────────────────────────────────────────────────── */

/* Sparse CSR matrix descriptor (no data owned — just metadata) */
typedef struct {
    hsize_t nnz;           /* number of non-zeros */
    hsize_t n_rows;        /* number of rows (cells for AnnData) */
    hsize_t n_cols;        /* number of columns (genes for AnnData) */
} sc_csr_info_t;

/* Column type in obs/var dataframes */
typedef enum {
    SC_COL_NUMERIC,
    SC_COL_INTEGER,
    SC_COL_STRING,
    SC_COL_CATEGORICAL,
    SC_COL_BOOLEAN,
    SC_COL_UNKNOWN
} sc_col_type_t;

/* Conversion direction */
typedef enum {
    SC_H5AD_TO_H5SEURAT,
    SC_H5SEURAT_TO_H5AD,
    SC_H5MU_TO_H5SEURAT,
    SC_H5SEURAT_TO_H5MU,
    SC_H5MU_TO_H5AD,
    SC_H5AD_TO_H5MU,
    SC_DIRECTION_UNKNOWN
} sc_direction_t;

/* Conversion options */
typedef struct {
    const char   *input_path;
    const char   *output_path;
    sc_direction_t direction;
    const char   *assay_name;      /* default: "RNA" */
    int           gzip_level;      /* default: SC_GZIP_LEVEL */
    int           verbose;
    int           overwrite;
} sc_opts_t;

/* ── Function declarations ──────────────────────────────────────────────────── */

/* Main conversion entry points (sc_h5ad.c) */
int sc_h5ad_to_h5seurat(const sc_opts_t *opts);
int sc_h5seurat_to_h5ad(const sc_opts_t *opts);

/* h5mu conversion entry points (sc_h5mu.c) */
int sc_h5mu_to_h5seurat(const sc_opts_t *opts);
int sc_h5seurat_to_h5mu(const sc_opts_t *opts);
int sc_h5mu_to_h5ad(const sc_opts_t *opts);
int sc_h5ad_to_h5mu(const sc_opts_t *opts);

/* Modality name mapping (sc_modality.c) */
const char *sc_modality_to_assay(const char *modality);
const char *sc_assay_to_modality(const char *assay);

/* Sparse matrix streaming (sc_sparse.c) */
int sc_read_csr_info(hid_t grp, sc_csr_info_t *info);
int sc_stream_csr_transpose(hid_t src_grp, hid_t dst_grp, int gzip_level);
int sc_stream_csr_copy(hid_t src_grp, hid_t dst_grp, int gzip_level);
int sc_write_csr_from_csc(hid_t src_grp, hid_t dst_grp,
                           hsize_t n_rows, hsize_t n_cols, int gzip_level);

/* DataFrame streaming (sc_dataframe.c) */
int sc_stream_obs_h5ad_to_h5seurat(hid_t src, hid_t dst, const char *assay);
int sc_stream_obs_h5seurat_to_h5ad(hid_t src, hid_t dst, const char *assay);
int sc_stream_var_h5ad_to_h5seurat(hid_t src, hid_t dst, const char *assay);
int sc_stream_var_h5seurat_to_h5ad(hid_t src, hid_t dst, const char *assay);
int sc_stream_df_group(hid_t src_grp, hid_t dst_grp, int gzip_level);

/* Metadata / group transfer (sc_groups.c) */
int sc_stream_obsm(hid_t src, hid_t dst, sc_direction_t dir);
int sc_stream_obsp(hid_t src, hid_t dst, sc_direction_t dir,
                    const char *assay);
int sc_stream_layers(hid_t src, hid_t dst, sc_direction_t dir, int gzip_level);
int sc_stream_uns(hid_t src, hid_t dst, sc_direction_t dir);
int sc_ensure_empty_groups(hid_t file, sc_direction_t dir);
int sc_copy_dataset_chunked(hid_t src_dset, hid_t dst_grp,
                             const char *name, int gzip_level);
int sc_copy_2d_transposed(hid_t src_dset, hid_t dst_grp,
                            const char *name, int gzip_level);
int sc_copy_group_recursive(hid_t src, hid_t dst, int gzip_level);
int sc_copy_group_h5ocopy(hid_t src_loc, const char *src_name,
                           hid_t dst_loc, const char *dst_name);

/* Utility (sc_util.c) */
int sc_set_str_attr(hid_t loc, const char *name, const char *value);
int sc_set_str_array_attr(hid_t loc, const char *name,
                           const char **values, hsize_t n);
int sc_set_int_array_attr(hid_t loc, const char *name,
                           const int64_t *values, hsize_t n);
int sc_get_str_attr(hid_t loc, const char *name, char *buf, size_t buflen);
int sc_get_str_array_attr(hid_t loc, const char *name,
                           char ***values, hsize_t *n);
void sc_free_str_array(char **values, hsize_t n);
hid_t sc_create_vlen_str_type(void);
int sc_has_group(hid_t loc, const char *name);
int sc_has_dataset(hid_t loc, const char *name);
hid_t sc_create_or_open_group(hid_t loc, const char *name);
int sc_get_encoding_type(hid_t loc, char *buf, size_t buflen);
int sc_copy_group_attrs(hid_t src, hid_t dst);

#endif /* SC_CONVERT_H */
