/*
 * sc_loom.c — Loom ↔ h5seurat conversion (+ chained loom ↔ h5ad/h5mu)
 *
 * Loom files are HDF5-based single-cell containers with a different layout
 * than h5ad/h5seurat.  Key structural differences:
 *
 *   - /matrix is DENSE (genes × cells), stored as a 2D HDF5 dataset
 *   - /layers/ has additional dense matrices (same orientation)
 *   - /row_attrs/ holds feature metadata (Gene names, etc.)
 *   - /col_attrs/ holds cell metadata (CellID, clusters, etc.)
 *   - /col_graphs/{name}/ stores cell-cell graphs in COO format (a, b, w)
 *   - /row_graphs/ stores gene-gene graphs in the same COO format
 *
 * All conversions route through h5seurat as the hub format:
 *   loom ↔ h5seurat  — direct (this file)
 *   loom ↔ h5ad      — via temp h5seurat, chaining existing converters
 *   loom ↔ h5mu      — via temp h5seurat, chaining existing converters
 */

#include "sc_convert.h"
#include <ctype.h>
#include <unistd.h>

/* ── Forward declarations ──────────────────────────────────────────────────── */

static int sc_dense_to_csc(hid_t src_dset, hid_t dst_grp,
                            hsize_t n_rows, hsize_t n_cols,
                            int gzip_level);
static int sc_csc_to_dense(hid_t src_grp, hid_t dst_loc,
                            const char *dst_name,
                            hsize_t n_rows, hsize_t n_cols,
                            int gzip_level);
static int sc_coo_to_csc(hid_t src_grp, hid_t dst_grp,
                           hsize_t n_dim, int gzip_level);
static int sc_csc_to_coo(hid_t src_grp, hid_t dst_grp,
                           int gzip_level);

/* ── Helper: convert via temporary h5seurat ────────────────────────────────── */

static int sc_convert_via_temp_h5seurat(
    const sc_opts_t *opts,
    int (*to_h5seurat)(const sc_opts_t *),
    int (*from_h5seurat)(const sc_opts_t *)
) {
    char tmp[1024];
    snprintf(tmp, sizeof(tmp), "%s.tmp.h5seurat", opts->output_path);

    /* Step 1: source → temp h5seurat */
    sc_opts_t step1 = *opts;
    step1.output_path = tmp;
    int rc = to_h5seurat(&step1);
    if (rc != SC_OK) {
        unlink(tmp);
        return rc;
    }

    /* Step 2: temp h5seurat → final output */
    sc_opts_t step2 = *opts;
    step2.input_path = tmp;
    rc = from_h5seurat(&step2);
    unlink(tmp);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Dense ↔ Sparse conversion helpers
 * ══════════════════════════════════════════════════════════════════════════════ */

/*
 * sc_dense_to_csc — Convert a dense [n_rows × n_cols] HDF5 dataset to CSC.
 *
 * The loom /matrix is [n_genes × n_cells]. h5seurat stores as CSC with
 * column = cell, row = gene.  We read one column (one cell) at a time
 * via hyperslab selection to limit memory usage.
 *
 * Two-pass approach:
 *   Pass 1: count non-zeros per column → build indptr
 *   Pass 2: fill data[] and indices[] arrays
 *
 * For very large matrices we stream in column batches to avoid O(nnz) memory.
 * However, for the indptr/indices/data arrays, we do need O(nnz) memory for
 * the final write.  For truly enormous matrices, the caller could chunk
 * further, but in practice single-cell matrices fit in RAM as sparse.
 */
static int sc_dense_to_csc(hid_t src_dset, hid_t dst_grp,
                            hsize_t n_rows, hsize_t n_cols,
                            int gzip_level)
{
    int rc = SC_OK;
    double *col_buf = NULL;
    int64_t *nnz_per_col = NULL;
    int64_t *indptr = NULL;
    double *data = NULL;
    int32_t *indices = NULL;

    if (n_rows == 0 || n_cols == 0) {
        /* Empty matrix: write trivial CSC */
        hsize_t zero = 0;
        hsize_t ptr_len = n_cols + 1;

        /* Empty data */
        {
            hid_t sp = H5Screate_simple(1, &zero, NULL);
            hid_t ds = H5Dcreate2(dst_grp, "data", H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(ds);
            H5Sclose(sp);
        }
        /* Empty indices */
        {
            hid_t sp = H5Screate_simple(1, &zero, NULL);
            hid_t ds = H5Dcreate2(dst_grp, "indices", H5T_NATIVE_INT32, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(ds);
            H5Sclose(sp);
        }
        /* All-zero indptr */
        {
            hid_t sp = H5Screate_simple(1, &ptr_len, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk = ptr_len;
            H5Pset_chunk(dcpl, 1, &chunk);
            hid_t ds = H5Dcreate2(dst_grp, "indptr", H5T_NATIVE_INT64, sp,
                                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
            int64_t *zeros = (int64_t *)calloc(ptr_len, sizeof(int64_t));
            if (zeros) {
                H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, zeros);
                free(zeros);
            }
            H5Dclose(ds);
            H5Pclose(dcpl);
            H5Sclose(sp);
        }
        return SC_OK;
    }

    /* Allocate column buffer */
    col_buf = (double *)malloc(n_rows * sizeof(double));
    if (!col_buf) return SC_ERR;

    /* Allocate nnz_per_col for pass 1 */
    nnz_per_col = (int64_t *)calloc(n_cols, sizeof(int64_t));
    if (!nnz_per_col) { free(col_buf); return SC_ERR; }

    hid_t fspace = H5Dget_space(src_dset);

    /* ── Pass 1: Count non-zeros per column ────────────────────────────── */
    hsize_t total_nnz = 0;
    for (hsize_t j = 0; j < n_cols; j++) {
        hsize_t start[2] = {0, j};
        hsize_t count[2] = {n_rows, 1};
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);

        hid_t mspace = H5Screate_simple(1, &n_rows, NULL);
        if (H5Dread(src_dset, H5T_NATIVE_DOUBLE, mspace, fspace,
                     H5P_DEFAULT, col_buf) < 0) {
            H5Sclose(mspace);
            rc = SC_ERR_HDF;
            goto cleanup;
        }
        H5Sclose(mspace);

        int64_t count_nz = 0;
        for (hsize_t i = 0; i < n_rows; i++) {
            if (col_buf[i] != 0.0)
                count_nz++;
        }
        nnz_per_col[j] = count_nz;
        total_nnz += (hsize_t)count_nz;
    }

    /* Build indptr */
    {
        hsize_t ptr_len = n_cols + 1;
        indptr = (int64_t *)malloc(ptr_len * sizeof(int64_t));
        if (!indptr) { rc = SC_ERR; goto cleanup; }

        indptr[0] = 0;
        for (hsize_t j = 0; j < n_cols; j++)
            indptr[j + 1] = indptr[j] + nnz_per_col[j];
    }

    /* Allocate data and indices for pass 2 */
    if (total_nnz > 0) {
        data = (double *)malloc(total_nnz * sizeof(double));
        indices = (int32_t *)malloc(total_nnz * sizeof(int32_t));
        if (!data || !indices) { rc = SC_ERR; goto cleanup; }
    }

    /* ── Pass 2: Fill data and indices ─────────────────────────────────── */
    {
        hsize_t pos = 0;
        for (hsize_t j = 0; j < n_cols; j++) {
            hsize_t start[2] = {0, j};
            hsize_t count[2] = {n_rows, 1};
            H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL,
                                 count, NULL);

            hid_t mspace = H5Screate_simple(1, &n_rows, NULL);
            if (H5Dread(src_dset, H5T_NATIVE_DOUBLE, mspace, fspace,
                         H5P_DEFAULT, col_buf) < 0) {
                H5Sclose(mspace);
                rc = SC_ERR_HDF;
                goto cleanup;
            }
            H5Sclose(mspace);

            for (hsize_t i = 0; i < n_rows; i++) {
                if (col_buf[i] != 0.0) {
                    data[pos] = col_buf[i];
                    indices[pos] = (int32_t)i;
                    pos++;
                }
            }
        }
    }

    H5Sclose(fspace);
    fspace = H5I_INVALID_HID;

    /* ── Write CSC arrays to destination group ─────────────────────────── */

    /* data */
    {
        hsize_t dim = total_nnz > 0 ? total_nnz : 0;
        if (dim == 0) {
            hid_t sp = H5Screate_simple(1, &dim, NULL);
            hid_t ds = H5Dcreate2(dst_grp, "data", H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(ds);
            H5Sclose(sp);
        } else {
            hid_t sp = H5Screate_simple(1, &dim, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk = dim < (hsize_t)SC_CHUNK_SIZE ? dim : (hsize_t)SC_CHUNK_SIZE;
            H5Pset_chunk(dcpl, 1, &chunk);
            if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);
            hid_t ds = H5Dcreate2(dst_grp, "data", H5T_NATIVE_DOUBLE, sp,
                                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
            H5Dclose(ds);
            H5Pclose(dcpl);
            H5Sclose(sp);
        }
    }

    /* indices */
    {
        hsize_t dim = total_nnz > 0 ? total_nnz : 0;
        if (dim == 0) {
            hid_t sp = H5Screate_simple(1, &dim, NULL);
            hid_t ds = H5Dcreate2(dst_grp, "indices", H5T_NATIVE_INT32, sp,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dclose(ds);
            H5Sclose(sp);
        } else {
            hid_t sp = H5Screate_simple(1, &dim, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk = dim < (hsize_t)SC_CHUNK_SIZE ? dim : (hsize_t)SC_CHUNK_SIZE;
            H5Pset_chunk(dcpl, 1, &chunk);
            if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);
            hid_t ds = H5Dcreate2(dst_grp, "indices", H5T_NATIVE_INT32, sp,
                                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
            H5Dclose(ds);
            H5Pclose(dcpl);
            H5Sclose(sp);
        }
    }

    /* indptr */
    {
        hsize_t ptr_len = n_cols + 1;
        hid_t sp = H5Screate_simple(1, &ptr_len, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = ptr_len < (hsize_t)SC_CHUNK_SIZE ? ptr_len : (hsize_t)SC_CHUNK_SIZE;
        H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);
        /* Write as int64 for indptr (matches h5seurat convention) */
        hid_t ds = H5Dcreate2(dst_grp, "indptr", H5T_NATIVE_INT64, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
        H5Dclose(ds);
        H5Pclose(dcpl);
        H5Sclose(sp);
    }

    rc = SC_OK;

cleanup:
    if (fspace != H5I_INVALID_HID) H5Sclose(fspace);
    free(col_buf);
    free(nnz_per_col);
    free(indptr);
    free(data);
    free(indices);
    return rc;
}

/*
 * sc_csc_to_dense — Convert CSC sparse group to dense [n_rows × n_cols] dataset.
 *
 * Reads indptr, indices, data arrays from src_grp and writes a dense 2D
 * dataset to dst_loc with name dst_name.
 */
static int sc_csc_to_dense(hid_t src_grp, hid_t dst_loc,
                            const char *dst_name,
                            hsize_t n_rows, hsize_t n_cols,
                            int gzip_level)
{
    if (n_rows == 0 || n_cols == 0) {
        /* Write empty 2D dataset */
        hsize_t dims[2] = {n_rows, n_cols};
        hid_t sp = H5Screate_simple(2, dims, NULL);
        hid_t ds = H5Dcreate2(dst_loc, dst_name, H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(ds);
        H5Sclose(sp);
        return SC_OK;
    }

    /* Read indptr */
    hid_t ptr_dset = H5Dopen2(src_grp, "indptr", H5P_DEFAULT);
    hid_t ptr_space = H5Dget_space(ptr_dset);
    hsize_t ptr_len;
    H5Sget_simple_extent_dims(ptr_space, &ptr_len, NULL);
    H5Sclose(ptr_space);

    int64_t *indptr = (int64_t *)malloc(ptr_len * sizeof(int64_t));
    if (!indptr) { H5Dclose(ptr_dset); return SC_ERR; }

    /* Read as int64 — h5seurat may store as int32 or int64, so use native conversion */
    hid_t ptr_type = H5Dget_type(ptr_dset);
    if (H5Tget_size(ptr_type) <= 4) {
        int32_t *tmp32 = (int32_t *)malloc(ptr_len * sizeof(int32_t));
        if (!tmp32) { H5Tclose(ptr_type); H5Dclose(ptr_dset); free(indptr); return SC_ERR; }
        H5Dread(ptr_dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp32);
        for (hsize_t i = 0; i < ptr_len; i++) indptr[i] = tmp32[i];
        free(tmp32);
    } else {
        H5Dread(ptr_dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
    }
    H5Tclose(ptr_type);
    H5Dclose(ptr_dset);

    hsize_t nnz = (hsize_t)indptr[ptr_len - 1];

    /* Read indices */
    int32_t *indices = NULL;
    if (nnz > 0) {
        hid_t idx_dset = H5Dopen2(src_grp, "indices", H5P_DEFAULT);
        indices = (int32_t *)malloc(nnz * sizeof(int32_t));
        if (!indices) { free(indptr); H5Dclose(idx_dset); return SC_ERR; }
        H5Dread(idx_dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
        H5Dclose(idx_dset);
    }

    /* Read data */
    double *vals = NULL;
    if (nnz > 0) {
        hid_t dat_dset = H5Dopen2(src_grp, "data", H5P_DEFAULT);
        vals = (double *)malloc(nnz * sizeof(double));
        if (!vals) { free(indptr); free(indices); H5Dclose(dat_dset); return SC_ERR; }
        H5Dread(dat_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
        H5Dclose(dat_dset);
    }

    /* Build dense matrix — write column by column to limit memory.
     * We allocate one column buffer and write via hyperslab. */
    hsize_t dims[2] = {n_rows, n_cols};
    hid_t dst_space = H5Screate_simple(2, dims, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    /* Use row-based chunking for decent compression and access patterns */
    {
        hsize_t chunk_rows = n_rows < 1024 ? n_rows : 1024;
        hsize_t chunk_cols = n_cols < 1024 ? n_cols : 1024;
        hsize_t chunks[2] = {chunk_rows, chunk_cols};
        H5Pset_chunk(dcpl, 2, chunks);
    }
    if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);

    hid_t dst_dset = H5Dcreate2(dst_loc, dst_name, H5T_NATIVE_DOUBLE,
                                 dst_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Initialize to zero — write a zero buffer for each column */
    double *col_buf = (double *)calloc(n_rows, sizeof(double));
    if (!col_buf) {
        free(indptr); free(indices); free(vals);
        H5Dclose(dst_dset); H5Pclose(dcpl); H5Sclose(dst_space);
        return SC_ERR;
    }

    hid_t fspace = H5Dget_space(dst_dset);

    for (hsize_t j = 0; j < n_cols; j++) {
        /* Zero the buffer */
        memset(col_buf, 0, n_rows * sizeof(double));

        /* Fill in non-zero entries from CSC */
        int64_t start_k = indptr[j];
        int64_t end_k = indptr[j + 1];
        for (int64_t k = start_k; k < end_k; k++) {
            int32_t row = indices[k];
            if (row >= 0 && (hsize_t)row < n_rows)
                col_buf[row] = vals[k];
        }

        /* Write column via hyperslab */
        hsize_t hs_start[2] = {0, j};
        hsize_t hs_count[2] = {n_rows, 1};
        H5Sselect_hyperslab(fspace, H5S_SELECT_SET, hs_start, NULL,
                             hs_count, NULL);
        hid_t mspace = H5Screate_simple(1, &n_rows, NULL);
        H5Dwrite(dst_dset, H5T_NATIVE_DOUBLE, mspace, fspace,
                  H5P_DEFAULT, col_buf);
        H5Sclose(mspace);
    }

    H5Sclose(fspace);
    H5Dclose(dst_dset);
    H5Pclose(dcpl);
    H5Sclose(dst_space);

    free(col_buf);
    free(indptr);
    free(indices);
    free(vals);
    return SC_OK;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  COO ↔ CSC conversion helpers (for loom col_graphs ↔ h5seurat graphs)
 * ══════════════════════════════════════════════════════════════════════════════ */

/*
 * Comparison function for COO entries: sort by column (b), then row (a).
 */
typedef struct {
    int32_t row;
    int32_t col;
    double  val;
} coo_entry_t;

static int coo_cmp(const void *x, const void *y) {
    const coo_entry_t *a = (const coo_entry_t *)x;
    const coo_entry_t *b = (const coo_entry_t *)y;
    if (a->col != b->col) return (a->col < b->col) ? -1 : 1;
    if (a->row != b->row) return (a->row < b->row) ? -1 : 1;
    return 0;
}

/*
 * sc_coo_to_csc — Convert loom COO graph (a, b, w) to h5seurat CSC.
 *
 * Loom COO: a = row indices, b = column indices, w = weights (all 0-based).
 * n_dim = number of cells (square matrix: n_dim × n_dim).
 *
 * Writes data, indices, indptr + dims attr to dst_grp.
 */
static int sc_coo_to_csc(hid_t src_grp, hid_t dst_grp,
                           hsize_t n_dim, int gzip_level)
{
    /* Read a (row indices) */
    hid_t a_dset = H5Dopen2(src_grp, "a", H5P_DEFAULT);
    if (a_dset < 0) return SC_ERR_HDF;
    hid_t a_space = H5Dget_space(a_dset);
    hsize_t nnz;
    H5Sget_simple_extent_dims(a_space, &nnz, NULL);
    H5Sclose(a_space);

    if (nnz == 0) {
        H5Dclose(a_dset);
        /* Write empty CSC */
        hsize_t zero = 0;
        hsize_t ptr_len = n_dim + 1;

        hid_t sp = H5Screate_simple(1, &zero, NULL);
        hid_t ds = H5Dcreate2(dst_grp, "data", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(ds);
        ds = H5Dcreate2(dst_grp, "indices", H5T_NATIVE_INT32, sp,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(ds);
        H5Sclose(sp);

        int64_t *zeros = (int64_t *)calloc(ptr_len, sizeof(int64_t));
        hid_t psp = H5Screate_simple(1, &ptr_len, NULL);
        ds = H5Dcreate2(dst_grp, "indptr", H5T_NATIVE_INT64, psp,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (zeros) {
            H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, zeros);
            free(zeros);
        }
        H5Dclose(ds);
        H5Sclose(psp);
        return SC_OK;
    }

    int32_t *a_arr = (int32_t *)malloc(nnz * sizeof(int32_t));
    if (!a_arr) { H5Dclose(a_dset); return SC_ERR; }
    H5Dread(a_dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_arr);
    H5Dclose(a_dset);

    /* Read b (column indices) */
    hid_t b_dset = H5Dopen2(src_grp, "b", H5P_DEFAULT);
    int32_t *b_arr = (int32_t *)malloc(nnz * sizeof(int32_t));
    if (!b_arr) { free(a_arr); H5Dclose(b_dset); return SC_ERR; }
    H5Dread(b_dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, b_arr);
    H5Dclose(b_dset);

    /* Read w (weights) */
    double *w_arr = (double *)malloc(nnz * sizeof(double));
    if (!w_arr) { free(a_arr); free(b_arr); return SC_ERR; }
    if (sc_has_dataset(src_grp, "w")) {
        hid_t w_dset = H5Dopen2(src_grp, "w", H5P_DEFAULT);
        H5Dread(w_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, w_arr);
        H5Dclose(w_dset);
    } else {
        /* No weights — default to 1.0 */
        for (hsize_t k = 0; k < nnz; k++) w_arr[k] = 1.0;
    }

    /* Build COO entry array and sort by column */
    coo_entry_t *entries = (coo_entry_t *)malloc(nnz * sizeof(coo_entry_t));
    if (!entries) {
        free(a_arr); free(b_arr); free(w_arr);
        return SC_ERR;
    }
    for (hsize_t k = 0; k < nnz; k++) {
        entries[k].row = a_arr[k];
        entries[k].col = b_arr[k];
        entries[k].val = w_arr[k];
    }
    free(a_arr);
    free(b_arr);
    free(w_arr);

    qsort(entries, (size_t)nnz, sizeof(coo_entry_t), coo_cmp);

    /* Build CSC arrays */
    hsize_t ptr_len = n_dim + 1;
    int64_t *indptr = (int64_t *)calloc(ptr_len, sizeof(int64_t));
    if (!indptr) { free(entries); return SC_ERR; }
    int32_t *indices = (int32_t *)malloc(nnz * sizeof(int32_t));
    if (!indices) { free(entries); free(indptr); return SC_ERR; }
    double *data = (double *)malloc(nnz * sizeof(double));
    if (!data) { free(entries); free(indptr); free(indices); return SC_ERR; }

    /* Count entries per column */
    for (hsize_t k = 0; k < nnz; k++) {
        int32_t col = entries[k].col;
        if (col >= 0 && (hsize_t)col < n_dim)
            indptr[col + 1]++;
    }
    /* Cumulative sum */
    for (hsize_t j = 1; j <= n_dim; j++)
        indptr[j] += indptr[j - 1];

    /* Fill data and indices (entries are already sorted by column) */
    for (hsize_t k = 0; k < nnz; k++) {
        indices[k] = entries[k].row;
        data[k] = entries[k].val;
    }
    free(entries);

    /* Write CSC arrays */
    int rc = SC_OK;

    /* data */
    {
        hid_t sp = H5Screate_simple(1, &nnz, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = nnz < (hsize_t)SC_CHUNK_SIZE ? nnz : (hsize_t)SC_CHUNK_SIZE;
        H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);
        hid_t ds = H5Dcreate2(dst_grp, "data", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(ds);
        H5Pclose(dcpl);
        H5Sclose(sp);
    }

    /* indices */
    {
        hid_t sp = H5Screate_simple(1, &nnz, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = nnz < (hsize_t)SC_CHUNK_SIZE ? nnz : (hsize_t)SC_CHUNK_SIZE;
        H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);
        hid_t ds = H5Dcreate2(dst_grp, "indices", H5T_NATIVE_INT32, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
        H5Dclose(ds);
        H5Pclose(dcpl);
        H5Sclose(sp);
    }

    /* indptr */
    {
        hid_t sp = H5Screate_simple(1, &ptr_len, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = ptr_len < (hsize_t)SC_CHUNK_SIZE ? ptr_len : (hsize_t)SC_CHUNK_SIZE;
        H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);
        hid_t ds = H5Dcreate2(dst_grp, "indptr", H5T_NATIVE_INT64, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
        H5Dclose(ds);
        H5Pclose(dcpl);
        H5Sclose(sp);
    }

    free(indptr);
    free(indices);
    free(data);
    return rc;
}

/*
 * sc_csc_to_coo — Convert h5seurat CSC graph to loom COO (a, b, w).
 *
 * Reads data, indices, indptr from src_grp.
 * Writes a (row indices), b (column indices), w (values) to dst_grp.
 */
static int sc_csc_to_coo(hid_t src_grp, hid_t dst_grp, int gzip_level)
{
    /* Read indptr */
    hid_t ptr_dset = H5Dopen2(src_grp, "indptr", H5P_DEFAULT);
    if (ptr_dset < 0) return SC_ERR_HDF;
    hid_t ptr_space = H5Dget_space(ptr_dset);
    hsize_t ptr_len;
    H5Sget_simple_extent_dims(ptr_space, &ptr_len, NULL);
    H5Sclose(ptr_space);

    int64_t *indptr = (int64_t *)malloc(ptr_len * sizeof(int64_t));
    if (!indptr) { H5Dclose(ptr_dset); return SC_ERR; }
    hid_t ptr_type = H5Dget_type(ptr_dset);
    if (H5Tget_size(ptr_type) <= 4) {
        int32_t *tmp32 = (int32_t *)malloc(ptr_len * sizeof(int32_t));
        if (!tmp32) { H5Tclose(ptr_type); H5Dclose(ptr_dset); free(indptr); return SC_ERR; }
        H5Dread(ptr_dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp32);
        for (hsize_t i = 0; i < ptr_len; i++) indptr[i] = tmp32[i];
        free(tmp32);
    } else {
        H5Dread(ptr_dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
    }
    H5Tclose(ptr_type);
    H5Dclose(ptr_dset);

    hsize_t n_cols = ptr_len - 1;
    hsize_t nnz = (hsize_t)indptr[n_cols];

    if (nnz == 0) {
        free(indptr);
        /* Write empty COO arrays */
        hsize_t zero = 0;
        hid_t sp = H5Screate_simple(1, &zero, NULL);
        hid_t ds;
        ds = H5Dcreate2(dst_grp, "a", H5T_NATIVE_INT32, sp,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(ds);
        ds = H5Dcreate2(dst_grp, "b", H5T_NATIVE_INT32, sp,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(ds);
        ds = H5Dcreate2(dst_grp, "w", H5T_NATIVE_DOUBLE, sp,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dclose(ds);
        H5Sclose(sp);
        return SC_OK;
    }

    /* Read indices (row indices in CSC) */
    hid_t idx_dset = H5Dopen2(src_grp, "indices", H5P_DEFAULT);
    int32_t *indices = (int32_t *)malloc(nnz * sizeof(int32_t));
    if (!indices) { free(indptr); H5Dclose(idx_dset); return SC_ERR; }
    H5Dread(idx_dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
    H5Dclose(idx_dset);

    /* Read data (values) */
    hid_t dat_dset = H5Dopen2(src_grp, "data", H5P_DEFAULT);
    double *data = (double *)malloc(nnz * sizeof(double));
    if (!data) { free(indptr); free(indices); H5Dclose(dat_dset); return SC_ERR; }
    H5Dread(dat_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dat_dset);

    /* Expand CSC to COO */
    int32_t *a_arr = (int32_t *)malloc(nnz * sizeof(int32_t));  /* row indices */
    if (!a_arr) { free(indptr); free(indices); free(data); return SC_ERR; }
    int32_t *b_arr = (int32_t *)malloc(nnz * sizeof(int32_t));  /* col indices */
    if (!b_arr) { free(indptr); free(indices); free(data); free(a_arr); return SC_ERR; }

    for (hsize_t j = 0; j < n_cols; j++) {
        for (int64_t k = indptr[j]; k < indptr[j + 1]; k++) {
            a_arr[k] = indices[k];
            b_arr[k] = (int32_t)j;
        }
    }
    free(indptr);
    free(indices);

    /* Write COO arrays */
    hid_t sp = H5Screate_simple(1, &nnz, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk = nnz < (hsize_t)SC_CHUNK_SIZE ? nnz : (hsize_t)SC_CHUNK_SIZE;
    H5Pset_chunk(dcpl, 1, &chunk);
    if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);

    /* a (row indices) */
    {
        hid_t ds = H5Dcreate2(dst_grp, "a", H5T_NATIVE_INT32, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_arr);
        H5Dclose(ds);
    }

    /* b (column indices) */
    {
        hid_t ds = H5Dcreate2(dst_grp, "b", H5T_NATIVE_INT32, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, b_arr);
        H5Dclose(ds);
    }

    /* w (values) */
    {
        hid_t ds = H5Dcreate2(dst_grp, "w", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(ds);
    }

    H5Pclose(dcpl);
    H5Sclose(sp);

    free(a_arr);
    free(b_arr);
    free(data);
    return SC_OK;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Loom → h5seurat
 * ══════════════════════════════════════════════════════════════════════════════ */

int sc_loom_to_h5seurat(const sc_opts_t *opts) {
    const char *assay = opts->assay_name ? opts->assay_name : "RNA";
    char assay_buf[SC_MAX_NAME_LEN] = {0};
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (loom → h5seurat)\n",
                opts->input_path, opts->output_path);

    /* Open loom file (read-only) */
    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        SC_MSG("Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    /* Create h5seurat file */
    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dst < 0) {
        SC_MSG("Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Check for SEURAT_ASSAY attribute on root (written by scConvert) */
    if (sc_get_str_attr(src, "SEURAT_ASSAY", assay_buf, sizeof(assay_buf)) == SC_OK
        && assay_buf[0] != '\0') {
        assay = assay_buf;
    }

    /* Set h5seurat root attributes */
    sc_set_str_attr(dst, "version", "5.3.0");
    sc_set_str_attr(dst, "project", "scConvert");
    sc_set_str_attr(dst, "active.assay", assay);

    /* Read matrix dimensions: /matrix is [n_genes × n_cells] */
    hsize_t n_genes = 0, n_cells = 0;
    if (sc_has_dataset(src, "matrix")) {
        hid_t mat_dset = H5Dopen2(src, "matrix", H5P_DEFAULT);
        hid_t mat_space = H5Dget_space(mat_dset);
        int ndims = H5Sget_simple_extent_ndims(mat_space);
        if (ndims == 2) {
            hsize_t mat_dims[2];
            H5Sget_simple_extent_dims(mat_space, mat_dims, NULL);
            n_genes = mat_dims[0];
            n_cells = mat_dims[1];
        }
        H5Sclose(mat_space);
        H5Dclose(mat_dset);
    }

    if (opts->verbose)
        SC_MSG("  Matrix dimensions: %llu genes × %llu cells\n",
                (unsigned long long)n_genes, (unsigned long long)n_cells);

    /* ── 1. Transfer /matrix → assays/{assay}/layers/data (dense → CSC) ── */
    if (opts->verbose) SC_MSG("  [1/6] Converting matrix to sparse...\n");
    if (n_genes > 0 && n_cells > 0) {
        hid_t assays = sc_create_or_open_group(dst, "assays");
        hid_t assay_grp = sc_create_or_open_group(assays, assay);

        /* Set assay key */
        {
            char key_buf[SC_MAX_NAME_LEN];
            size_t len = strlen(assay);
            if (len > sizeof(key_buf) - 2) len = sizeof(key_buf) - 2;
            for (size_t k = 0; k < len; k++)
                key_buf[k] = (char)tolower((unsigned char)assay[k]);
            key_buf[len] = '_';
            key_buf[len + 1] = '\0';
            sc_set_str_attr(assay_grp, "key", key_buf);
        }
        sc_set_str_attr(assay_grp, "s4class", "SeuratObject::Assay5");

        hid_t layers_grp = sc_create_or_open_group(assay_grp, "layers");

        hid_t mat_dset = H5Dopen2(src, "matrix", H5P_DEFAULT);

        /* Convert dense [n_genes × n_cells] → CSC (genes×cells) */
        hid_t dst_data = sc_create_or_open_group(layers_grp, "data");
        rc = sc_dense_to_csc(mat_dset, dst_data, n_genes, n_cells, gzip);

        /* Set dims = [n_genes, n_cells] (h5seurat convention) */
        int64_t dims[2] = {(int64_t)n_genes, (int64_t)n_cells};
        sc_set_int_array_attr(dst_data, "dims", dims, 2);

        H5Gclose(dst_data);
        H5Dclose(mat_dset);

        /* Also write to counts (so Seurat has both layers) */
        if (rc == SC_OK) {
            mat_dset = H5Dopen2(src, "matrix", H5P_DEFAULT);
            hid_t dst_counts = sc_create_or_open_group(layers_grp, "counts");
            rc = sc_dense_to_csc(mat_dset, dst_counts, n_genes, n_cells, gzip);
            sc_set_int_array_attr(dst_counts, "dims", dims, 2);
            H5Gclose(dst_counts);
            H5Dclose(mat_dset);
        }

        H5Gclose(layers_grp);
        H5Gclose(assay_grp);
        H5Gclose(assays);

        if (rc != SC_OK) goto cleanup;
    }

    /* ── 2. Transfer /row_attrs/Gene → assays/{assay}/features ─────────── */
    if (opts->verbose) SC_MSG("  [2/6] Transferring features...\n");
    {
        hid_t assays = sc_create_or_open_group(dst, "assays");
        hid_t assay_grp = sc_create_or_open_group(assays, assay);

        if (sc_has_group(src, "row_attrs")) {
            hid_t ra = H5Gopen2(src, "row_attrs", H5P_DEFAULT);

            /* Try Gene, then gene, then _index, then first string dataset */
            const char *feat_names[] = {"Gene", "gene", "Gene_Name",
                                         "gene_name", "_index", "var_names"};
            int found = 0;
            for (int i = 0; i < 6 && !found; i++) {
                if (sc_has_dataset(ra, feat_names[i])) {
                    hid_t src_dset = H5Dopen2(ra, feat_names[i], H5P_DEFAULT);
                    sc_copy_dataset_chunked(src_dset, assay_grp, "features", gzip);
                    H5Dclose(src_dset);
                    found = 1;
                }
            }

            /* If not found by name, try first string dataset in row_attrs */
            if (!found) {
                hsize_t n_members;
                H5Gget_num_objs(ra, &n_members);
                for (hsize_t i = 0; i < n_members && !found; i++) {
                    char name[SC_MAX_NAME_LEN];
                    H5Gget_objname_by_idx(ra, i, name, sizeof(name));
                    int obj_type = H5Gget_objtype_by_idx(ra, i);
                    if (obj_type == H5G_DATASET) {
                        hid_t ds = H5Dopen2(ra, name, H5P_DEFAULT);
                        hid_t dtype = H5Dget_type(ds);
                        if (H5Tget_class(dtype) == H5T_STRING) {
                            sc_copy_dataset_chunked(ds, assay_grp, "features", gzip);
                            found = 1;
                        }
                        H5Tclose(dtype);
                        H5Dclose(ds);
                    }
                }
            }
            H5Gclose(ra);
        }

        H5Gclose(assay_grp);
        H5Gclose(assays);
    }

    /* ── 3. Transfer /col_attrs/CellID → cell.names ────────────────────── */
    if (opts->verbose) SC_MSG("  [3/6] Transferring cell names...\n");
    if (sc_has_group(src, "col_attrs")) {
        hid_t ca = H5Gopen2(src, "col_attrs", H5P_DEFAULT);

        const char *cell_names[] = {"CellID", "cellid", "cell_id",
                                     "barcode", "barcodes", "obs_names"};
        int found = 0;
        for (int i = 0; i < 6 && !found; i++) {
            if (sc_has_dataset(ca, cell_names[i])) {
                hid_t src_dset = H5Dopen2(ca, cell_names[i], H5P_DEFAULT);
                sc_copy_dataset_chunked(src_dset, dst, "cell.names", gzip);
                H5Dclose(src_dset);
                found = 1;
            }
        }

        /* Fallback: first string dataset */
        if (!found) {
            hsize_t n_members;
            H5Gget_num_objs(ca, &n_members);
            for (hsize_t i = 0; i < n_members && !found; i++) {
                char name[SC_MAX_NAME_LEN];
                H5Gget_objname_by_idx(ca, i, name, sizeof(name));
                int obj_type = H5Gget_objtype_by_idx(ca, i);
                if (obj_type == H5G_DATASET) {
                    hid_t ds = H5Dopen2(ca, name, H5P_DEFAULT);
                    hid_t dtype = H5Dget_type(ds);
                    if (H5Tget_class(dtype) == H5T_STRING) {
                        sc_copy_dataset_chunked(ds, dst, "cell.names", gzip);
                        found = 1;
                    }
                    H5Tclose(dtype);
                    H5Dclose(ds);
                }
            }
        }
        H5Gclose(ca);
    }

    /* ── 4. Transfer /col_attrs → meta.data ────────────────────────────── */
    if (opts->verbose) SC_MSG("  [4/6] Transferring cell metadata...\n");
    if (sc_has_group(src, "col_attrs")) {
        hid_t ca = H5Gopen2(src, "col_attrs", H5P_DEFAULT);
        hid_t meta = sc_create_or_open_group(dst, "meta.data");

        /* Collect column names for colnames attribute */
        char *col_names[SC_MAX_COLUMNS];
        int n_cols = 0;

        hsize_t n_members;
        H5Gget_num_objs(ca, &n_members);

        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(ca, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(ca, i);

            /* Skip CellID (already stored as cell.names) */
            if (strcmp(name, "CellID") == 0 || strcmp(name, "cellid") == 0 ||
                strcmp(name, "cell_id") == 0 || strcmp(name, "barcode") == 0 ||
                strcmp(name, "barcodes") == 0 || strcmp(name, "obs_names") == 0)
                continue;

            if (obj_type == H5G_DATASET) {
                hid_t src_dset = H5Dopen2(ca, name, H5P_DEFAULT);
                hid_t src_space = H5Dget_space(src_dset);
                int ndims = H5Sget_simple_extent_ndims(src_space);
                H5Sclose(src_space);

                if (ndims <= 1) {
                    /* 1D dataset: copy directly as metadata column */
                    sc_copy_dataset_chunked(src_dset, meta, name, gzip);
                    if (n_cols < SC_MAX_COLUMNS)
                        col_names[n_cols++] = strdup(name);
                }
                /* Skip 2D+ datasets in col_attrs (embeddings handled separately) */
                H5Dclose(src_dset);
            }
        }

        /* Set colnames attribute */
        if (n_cols > 0) {
            sc_set_str_array_attr(meta, "colnames",
                                  (const char **)col_names, (hsize_t)n_cols);
            for (int i = 0; i < n_cols; i++) free(col_names[i]);
        }

        H5Gclose(meta);
        H5Gclose(ca);
    }

    /* ── 5. Transfer /col_graphs → graphs ──────────────────────────────── */
    if (opts->verbose) SC_MSG("  [5/6] Transferring graphs...\n");
    if (sc_has_group(src, "col_graphs")) {
        hid_t cg = H5Gopen2(src, "col_graphs", H5P_DEFAULT);
        hid_t graphs = sc_create_or_open_group(dst, "graphs");

        hsize_t n_members;
        H5Gget_num_objs(cg, &n_members);

        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(cg, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(cg, i);

            if (obj_type == H5G_GROUP) {
                hid_t src_graph = H5Gopen2(cg, name, H5P_DEFAULT);

                /* Check for COO arrays (a, b) */
                if (sc_has_dataset(src_graph, "a") &&
                    sc_has_dataset(src_graph, "b")) {
                    /* Prefix graph name with assay */
                    char dst_name[SC_MAX_NAME_LEN];
                    snprintf(dst_name, sizeof(dst_name), "%s_%s", assay, name);

                    hid_t dst_graph = sc_create_or_open_group(graphs, dst_name);
                    rc = sc_coo_to_csc(src_graph, dst_graph, n_cells, gzip);

                    /* Set dims and assay.used attributes */
                    int64_t gdims[2] = {(int64_t)n_cells, (int64_t)n_cells};
                    sc_set_int_array_attr(dst_graph, "dims", gdims, 2);
                    sc_set_str_attr(dst_graph, "assay.used", assay);

                    H5Gclose(dst_graph);
                    if (rc != SC_OK) {
                        H5Gclose(src_graph);
                        H5Gclose(graphs);
                        H5Gclose(cg);
                        goto cleanup;
                    }
                }
                H5Gclose(src_graph);
            }
        }

        H5Gclose(graphs);
        H5Gclose(cg);
    }

    /* ── 6. Transfer /reductions → reductions (scConvert extension) ────── */
    if (opts->verbose) SC_MSG("  [6/6] Transferring reductions...\n");
    if (sc_has_group(src, "reductions")) {
        /* scConvert loom extension: /reductions/{name}/embeddings */
        int copy_rc = sc_copy_group_h5ocopy(src, "reductions", dst, "reductions");
        if (copy_rc != SC_OK) {
            hid_t src_red = H5Gopen2(src, "reductions", H5P_DEFAULT);
            hid_t dst_red = sc_create_or_open_group(dst, "reductions");
            sc_copy_group_recursive(src_red, dst_red, gzip);
            H5Gclose(dst_red);
            H5Gclose(src_red);
        }
    }

    /* ── 7. Transfer /layers → additional expression layers ────────────── */
    if (sc_has_group(src, "layers")) {
        hid_t src_layers = H5Gopen2(src, "layers", H5P_DEFAULT);
        hsize_t n_members;
        H5Gget_num_objs(src_layers, &n_members);

        if (n_members > 0 && opts->verbose)
            SC_MSG("  [+] Transferring %llu additional layers...\n",
                    (unsigned long long)n_members);

        hid_t assays = sc_create_or_open_group(dst, "assays");
        hid_t assay_grp = sc_create_or_open_group(assays, assay);
        hid_t layers_grp = sc_create_or_open_group(assay_grp, "layers");

        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_layers, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(src_layers, i);

            if (obj_type == H5G_DATASET) {
                /* Skip if this layer name already exists (e.g. "data" or "counts") */
                if (sc_has_group(layers_grp, name) || sc_has_dataset(layers_grp, name))
                    continue;

                hid_t src_dset = H5Dopen2(src_layers, name, H5P_DEFAULT);
                hid_t src_space = H5Dget_space(src_dset);
                int ndims = H5Sget_simple_extent_ndims(src_space);

                if (ndims == 2) {
                    hsize_t layer_dims[2];
                    H5Sget_simple_extent_dims(src_space, layer_dims, NULL);
                    H5Sclose(src_space);

                    if (layer_dims[0] == n_genes && layer_dims[1] == n_cells) {
                        /* Dense [n_genes × n_cells] → sparse CSC */
                        hid_t dst_layer = sc_create_or_open_group(layers_grp, name);
                        rc = sc_dense_to_csc(src_dset, dst_layer,
                                              n_genes, n_cells, gzip);
                        int64_t dims[2] = {(int64_t)n_genes, (int64_t)n_cells};
                        sc_set_int_array_attr(dst_layer, "dims", dims, 2);
                        H5Gclose(dst_layer);
                    }
                } else {
                    H5Sclose(src_space);
                }
                H5Dclose(src_dset);

                if (rc != SC_OK) break;
            }
        }

        H5Gclose(layers_grp);
        H5Gclose(assay_grp);
        H5Gclose(assays);
        H5Gclose(src_layers);

        if (rc != SC_OK) goto cleanup;
    }

    /* ── 8. Create required empty groups for h5seurat ──────────────────── */
    {
        const char *required[] = {"tools", "commands", "images",
                                   "neighbors", "misc", "graphs",
                                   "reductions"};
        for (int i = 0; i < 7; i++) {
            if (!sc_has_group(dst, required[i])) {
                hid_t g = H5Gcreate2(dst, required[i], H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
                H5Gclose(g);
            }
        }
    }

    /* ── 9. Create active.ident ────────────────────────────────────────── */
    {
        hid_t ai = sc_create_or_open_group(dst, "active.ident");

        /* levels: ["scConvert"] */
        if (!sc_has_dataset(ai, "levels")) {
            hsize_t one = 1;
            hid_t space = H5Screate_simple(1, &one, NULL);
            hid_t tid = sc_create_vlen_str_type();
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(dcpl, 1, &one);
            hid_t ds = H5Dcreate2(ai, "levels", tid, space,
                                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
            const char *level = "scConvert";
            H5Dwrite(ds, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &level);
            H5Dclose(ds);
            H5Pclose(dcpl);
            H5Tclose(tid);
            H5Sclose(space);
        }

        /* values: n_cells × 1 (all = 1, 1-based) */
        if (!sc_has_dataset(ai, "values") && n_cells > 0) {
            int *vals = (int *)malloc(n_cells * sizeof(int));
            if (!vals) { H5Gclose(ai); goto cleanup; }
            for (hsize_t c = 0; c < n_cells; c++) vals[c] = 1;

            hid_t space = H5Screate_simple(1, &n_cells, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            hsize_t chunk = n_cells;
            H5Pset_chunk(dcpl, 1, &chunk);
            hid_t ds = H5Dcreate2(ai, "values", H5T_NATIVE_INT, space,
                                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
            H5Dclose(ds);
            H5Pclose(dcpl);
            H5Sclose(space);
            free(vals);
        }

        H5Gclose(ai);
    }

    if (opts->verbose && rc == SC_OK)
        SC_MSG("[scConvert] Done.\n");

cleanup:
    H5Fclose(dst);
    H5Fclose(src);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  h5seurat → Loom
 * ══════════════════════════════════════════════════════════════════════════════ */

int sc_h5seurat_to_loom(const sc_opts_t *opts) {
    const char *assay = opts->assay_name ? opts->assay_name : "RNA";
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5seurat → loom)\n",
                opts->input_path, opts->output_path);

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        SC_MSG("Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dst < 0) {
        SC_MSG("Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Read active.assay from source if available */
    {
        char active_buf[SC_MAX_NAME_LEN] = {0};
        if (sc_get_str_attr(src, "active.assay", active_buf, sizeof(active_buf)) == SC_OK
            && active_buf[0] != '\0'
            && strcmp(opts->assay_name, "RNA") == 0) {
            /* Only override if user didn't explicitly set an assay */
            assay = strdup(active_buf);
        }
    }

    /* Get matrix dimensions */
    hsize_t n_genes = 0, n_cells = 0;

    /* Try to get n_cells from cell.names */
    if (sc_has_dataset(src, "cell.names")) {
        hid_t cn = H5Dopen2(src, "cell.names", H5P_DEFAULT);
        hid_t cn_sp = H5Dget_space(cn);
        H5Sget_simple_extent_dims(cn_sp, &n_cells, NULL);
        H5Sclose(cn_sp);
        H5Dclose(cn);
    }

    /* Try to get n_genes from features */
    {
        char feat_path[512];
        snprintf(feat_path, sizeof(feat_path), "assays/%s/features", assay);
        if (sc_has_dataset(src, feat_path)) {
            hid_t fn = H5Dopen2(src, feat_path, H5P_DEFAULT);
            hid_t fn_sp = H5Dget_space(fn);
            H5Sget_simple_extent_dims(fn_sp, &n_genes, NULL);
            H5Sclose(fn_sp);
            H5Dclose(fn);
        }
    }

    if (opts->verbose)
        SC_MSG("  Matrix dimensions: %llu genes × %llu cells\n",
                (unsigned long long)n_genes, (unsigned long long)n_cells);

    /* ── 1. Transfer expression matrix → /matrix (CSC → dense) ─────────── */
    if (opts->verbose) SC_MSG("  [1/7] Converting sparse to dense matrix...\n");
    {
        /* Find sparse data group — try v5 path first, then v4 */
        char data_path[512];
        snprintf(data_path, sizeof(data_path), "assays/%s/layers/data", assay);
        if (!sc_has_group(src, data_path)) {
            snprintf(data_path, sizeof(data_path), "assays/%s/data", assay);
        }
        if (!sc_has_group(src, data_path)) {
            /* Try counts if no data */
            snprintf(data_path, sizeof(data_path), "assays/%s/layers/counts", assay);
            if (!sc_has_group(src, data_path)) {
                snprintf(data_path, sizeof(data_path), "assays/%s/counts", assay);
            }
        }

        if (sc_has_group(src, data_path)) {
            hid_t src_data = H5Gopen2(src, data_path, H5P_DEFAULT);

            /* Get dims from the sparse group if we don't have them yet */
            if (n_genes == 0 || n_cells == 0) {
                sc_csr_info_t info = {0};
                sc_read_csr_info(src_data, &info);
                /* h5seurat dims: [n_genes, n_cells] */
                n_genes = info.n_rows;
                n_cells = info.n_cols;
            }

            rc = sc_csc_to_dense(src_data, dst, "matrix",
                                  n_genes, n_cells, gzip);
            H5Gclose(src_data);
            if (rc != SC_OK) goto cleanup_loom;
        } else if (sc_has_dataset(src, data_path)) {
            /* Dense data layer — copy directly */
            hid_t src_dset = H5Dopen2(src, data_path, H5P_DEFAULT);
            sc_copy_dataset_chunked(src_dset, dst, "matrix", gzip);
            H5Dclose(src_dset);
        } else {
            /* No matrix found — write empty */
            if (n_genes > 0 && n_cells > 0) {
                hsize_t dims[2] = {n_genes, n_cells};
                hid_t sp = H5Screate_simple(2, dims, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                hsize_t chunks[2] = {n_genes < 1024 ? n_genes : 1024,
                                      n_cells < 1024 ? n_cells : 1024};
                H5Pset_chunk(dcpl, 2, chunks);
                if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
                hid_t ds = H5Dcreate2(dst, "matrix", H5T_NATIVE_DOUBLE, sp,
                                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
                /* Zero-fill */
                double *zeros = (double *)calloc(n_cells, sizeof(double));
                if (zeros) {
                    for (hsize_t i = 0; i < n_genes; i++) {
                        hsize_t start[2] = {i, 0};
                        hsize_t count[2] = {1, n_cells};
                        hid_t fsp = H5Dget_space(ds);
                        H5Sselect_hyperslab(fsp, H5S_SELECT_SET, start, NULL, count, NULL);
                        hid_t msp = H5Screate_simple(1, &n_cells, NULL);
                        H5Dwrite(ds, H5T_NATIVE_DOUBLE, msp, fsp, H5P_DEFAULT, zeros);
                        H5Sclose(msp);
                        H5Sclose(fsp);
                    }
                    free(zeros);
                }
                H5Dclose(ds);
                H5Pclose(dcpl);
                H5Sclose(sp);
            }
        }
    }

    /* ── 2. Transfer cell.names → /col_attrs/CellID ────────────────────── */
    if (opts->verbose) SC_MSG("  [2/7] Writing cell names...\n");
    {
        hid_t col_attrs = sc_create_or_open_group(dst, "col_attrs");

        if (sc_has_dataset(src, "cell.names")) {
            hid_t cn = H5Dopen2(src, "cell.names", H5P_DEFAULT);
            sc_copy_dataset_chunked(cn, col_attrs, "CellID", gzip);
            H5Dclose(cn);
        }

        H5Gclose(col_attrs);
    }

    /* ── 3. Transfer features → /row_attrs/Gene ────────────────────────── */
    if (opts->verbose) SC_MSG("  [3/7] Writing feature names...\n");
    {
        hid_t row_attrs = sc_create_or_open_group(dst, "row_attrs");

        char feat_path[512];
        snprintf(feat_path, sizeof(feat_path), "assays/%s/features", assay);
        if (sc_has_dataset(src, feat_path)) {
            hid_t fn = H5Dopen2(src, feat_path, H5P_DEFAULT);
            sc_copy_dataset_chunked(fn, row_attrs, "Gene", gzip);
            H5Dclose(fn);
        }

        H5Gclose(row_attrs);
    }

    /* ── 4. Transfer meta.data → /col_attrs (each column as dataset) ──── */
    if (opts->verbose) SC_MSG("  [4/7] Writing cell metadata...\n");
    if (sc_has_group(src, "meta.data")) {
        hid_t meta = H5Gopen2(src, "meta.data", H5P_DEFAULT);
        hid_t col_attrs = sc_create_or_open_group(dst, "col_attrs");

        hsize_t n_members;
        H5Gget_num_objs(meta, &n_members);

        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(meta, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(meta, i);

            /* Skip _index (that's CellID, already handled) */
            if (strcmp(name, "_index") == 0) continue;

            if (obj_type == H5G_DATASET) {
                /* Direct copy of numeric/string datasets */
                if (!sc_has_dataset(col_attrs, name)) {
                    hid_t src_dset = H5Dopen2(meta, name, H5P_DEFAULT);
                    sc_copy_dataset_chunked(src_dset, col_attrs, name, gzip);
                    H5Dclose(src_dset);
                }
            } else if (obj_type == H5G_GROUP) {
                /* Factor group (levels/values) — expand to string array for loom.
                 * Loom doesn't support factor encoding, so we store the
                 * string representation. */
                hid_t fac = H5Gopen2(meta, name, H5P_DEFAULT);

                if (sc_has_dataset(fac, "levels") && sc_has_dataset(fac, "values")) {
                    /* Read levels */
                    hid_t lev_dset = H5Dopen2(fac, "levels", H5P_DEFAULT);
                    hid_t lev_space = H5Dget_space(lev_dset);
                    hsize_t n_levels;
                    H5Sget_simple_extent_dims(lev_space, &n_levels, NULL);
                    H5Sclose(lev_space);

                    hid_t str_type = sc_create_vlen_str_type();
                    char **levels = (char **)calloc(n_levels, sizeof(char *));
                    if (!levels) {
                        H5Tclose(str_type); H5Dclose(lev_dset); H5Gclose(fac);
                        continue;
                    }
                    H5Dread(lev_dset, str_type, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, levels);
                    H5Dclose(lev_dset);

                    /* Read values (1-based) */
                    hid_t val_dset = H5Dopen2(fac, "values", H5P_DEFAULT);
                    hid_t val_space = H5Dget_space(val_dset);
                    hsize_t n_vals;
                    H5Sget_simple_extent_dims(val_space, &n_vals, NULL);
                    H5Sclose(val_space);

                    int *vals = (int *)malloc(n_vals * sizeof(int));
                    if (!vals) {
                        hid_t rs = H5Screate_simple(1, &n_levels, NULL);
                        H5Treclaim(str_type, rs, H5P_DEFAULT, levels);
                        H5Sclose(rs);
                        free(levels); H5Tclose(str_type);
                        H5Dclose(val_dset); H5Gclose(fac);
                        continue;
                    }
                    H5Dread(val_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, vals);
                    H5Dclose(val_dset);

                    /* Expand to string array */
                    char **str_arr = (char **)calloc(n_vals, sizeof(char *));
                    if (!str_arr) {
                        free(vals);
                        hid_t rs = H5Screate_simple(1, &n_levels, NULL);
                        H5Treclaim(str_type, rs, H5P_DEFAULT, levels);
                        H5Sclose(rs);
                        free(levels); H5Tclose(str_type); H5Gclose(fac);
                        continue;
                    }
                    for (hsize_t j = 0; j < n_vals; j++) {
                        int idx = vals[j] - 1;  /* 1-based → 0-based */
                        if (idx >= 0 && (hsize_t)idx < n_levels && levels[idx])
                            str_arr[j] = levels[idx];
                        else
                            str_arr[j] = "";
                    }

                    /* Write expanded string array to col_attrs */
                    hid_t sp = H5Screate_simple(1, &n_vals, NULL);
                    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                    hsize_t chunk = n_vals;
                    H5Pset_chunk(dcpl, 1, &chunk);
                    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);

                    hid_t ds = H5Dcreate2(col_attrs, name, str_type, sp,
                                           H5P_DEFAULT, dcpl, H5P_DEFAULT);
                    H5Dwrite(ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, str_arr);
                    H5Dclose(ds);
                    H5Pclose(dcpl);
                    H5Sclose(sp);

                    free(str_arr);
                    free(vals);

                    /* Free HDF5-allocated level strings */
                    hid_t reclaim_space = H5Screate_simple(1, &n_levels, NULL);
                    H5Treclaim(str_type, reclaim_space, H5P_DEFAULT, levels);
                    H5Sclose(reclaim_space);
                    free(levels);
                    H5Tclose(str_type);
                }

                H5Gclose(fac);
            }
        }

        H5Gclose(col_attrs);
        H5Gclose(meta);
    }

    /* ── 5. Transfer graphs → /col_graphs (CSC → COO) ─────────────────── */
    if (opts->verbose) SC_MSG("  [5/7] Transferring graphs...\n");
    if (sc_has_group(src, "graphs")) {
        hid_t graphs = H5Gopen2(src, "graphs", H5P_DEFAULT);
        hid_t col_graphs = sc_create_or_open_group(dst, "col_graphs");

        hsize_t n_members;
        H5Gget_num_objs(graphs, &n_members);

        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(graphs, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(graphs, i);

            if (obj_type == H5G_GROUP) {
                hid_t src_graph = H5Gopen2(graphs, name, H5P_DEFAULT);

                /* Check it's a sparse graph (has data/indices/indptr) */
                if (sc_has_dataset(src_graph, "data") &&
                    sc_has_dataset(src_graph, "indices") &&
                    sc_has_dataset(src_graph, "indptr")) {

                    /* Strip assay prefix for loom name
                     * e.g. "RNA_snn" → "snn" */
                    const char *loom_name = name;
                    size_t prefix_len = strlen(assay);
                    if (strncmp(name, assay, prefix_len) == 0 &&
                        name[prefix_len] == '_') {
                        loom_name = name + prefix_len + 1;
                    }

                    hid_t dst_graph = sc_create_or_open_group(col_graphs,
                                                                loom_name);
                    rc = sc_csc_to_coo(src_graph, dst_graph, gzip);
                    H5Gclose(dst_graph);

                    if (rc != SC_OK) {
                        H5Gclose(src_graph);
                        H5Gclose(col_graphs);
                        H5Gclose(graphs);
                        goto cleanup_loom;
                    }
                }
                H5Gclose(src_graph);
            }
        }

        H5Gclose(col_graphs);
        H5Gclose(graphs);
    }

    /* ── 6. Transfer reductions → /reductions (scConvert extension) ────── */
    if (opts->verbose) SC_MSG("  [6/7] Transferring reductions...\n");
    if (sc_has_group(src, "reductions")) {
        int copy_rc = sc_copy_group_h5ocopy(src, "reductions", dst, "reductions");
        if (copy_rc != SC_OK) {
            hid_t src_red = H5Gopen2(src, "reductions", H5P_DEFAULT);
            hid_t dst_red = sc_create_or_open_group(dst, "reductions");
            sc_copy_group_recursive(src_red, dst_red, gzip);
            H5Gclose(dst_red);
            H5Gclose(src_red);
        }
    }

    /* ── 7. Transfer additional layers → /layers ───────────────────────── */
    if (opts->verbose) SC_MSG("  [7/7] Transferring layers...\n");
    {
        hid_t layers_dst = sc_create_or_open_group(dst, "layers");

        /* Check v5 path: assays/{assay}/layers/ */
        char layers_path[512];
        snprintf(layers_path, sizeof(layers_path), "assays/%s/layers", assay);

        if (sc_has_group(src, layers_path)) {
            hid_t src_layers = H5Gopen2(src, layers_path, H5P_DEFAULT);
            hsize_t n_members;
            H5Gget_num_objs(src_layers, &n_members);

            for (hsize_t i = 0; i < n_members; i++) {
                char name[SC_MAX_NAME_LEN];
                H5Gget_objname_by_idx(src_layers, i, name, sizeof(name));
                int obj_type = H5Gget_objtype_by_idx(src_layers, i);

                /* Skip "data" (already written as /matrix) and "counts"
                 * (also already written as /matrix) */
                if (strcmp(name, "data") == 0 || strcmp(name, "counts") == 0)
                    continue;

                if (obj_type == H5G_GROUP) {
                    /* Sparse layer → dense */
                    hid_t src_layer = H5Gopen2(src_layers, name, H5P_DEFAULT);
                    if (sc_has_dataset(src_layer, "data") &&
                        sc_has_dataset(src_layer, "indices") &&
                        sc_has_dataset(src_layer, "indptr")) {
                        rc = sc_csc_to_dense(src_layer, layers_dst, name,
                                              n_genes, n_cells, gzip);
                    }
                    H5Gclose(src_layer);
                } else if (obj_type == H5G_DATASET) {
                    /* Dense layer — copy directly */
                    hid_t src_dset = H5Dopen2(src_layers, name, H5P_DEFAULT);
                    sc_copy_dataset_chunked(src_dset, layers_dst, name, gzip);
                    H5Dclose(src_dset);
                }

                if (rc != SC_OK) break;
            }
            H5Gclose(src_layers);
        }

        H5Gclose(layers_dst);
        if (rc != SC_OK) goto cleanup_loom;
    }

    /* ── 8. Ensure required empty groups ───────────────────────────────── */
    {
        const char *required[] = {"col_attrs", "row_attrs", "col_graphs",
                                   "row_graphs", "layers"};
        for (int i = 0; i < 5; i++) {
            if (!sc_has_group(dst, required[i])) {
                hid_t g = H5Gcreate2(dst, required[i], H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
                H5Gclose(g);
            }
        }
    }

    /* ── 9. Set loom root attributes ───────────────────────────────────── */
    sc_set_str_attr(dst, "LOOM_SPEC_VERSION", "3.0.0");
    sc_set_str_attr(dst, "SEURAT_ASSAY", assay);

    if (opts->verbose && rc == SC_OK)
        SC_MSG("[scConvert] Done.\n");

cleanup_loom:
    H5Fclose(dst);
    H5Fclose(src);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Composite conversions via temp h5seurat
 * ══════════════════════════════════════════════════════════════════════════════ */

int sc_loom_to_h5ad(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (loom → h5ad via h5seurat)\n",
                opts->input_path, opts->output_path);

    return sc_convert_via_temp_h5seurat(opts,
                                         sc_loom_to_h5seurat,
                                         sc_h5seurat_to_h5ad);
}

int sc_h5ad_to_loom(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5ad → loom via h5seurat)\n",
                opts->input_path, opts->output_path);

    return sc_convert_via_temp_h5seurat(opts,
                                         sc_h5ad_to_h5seurat,
                                         sc_h5seurat_to_loom);
}

int sc_loom_to_h5mu(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (loom → h5mu via h5seurat)\n",
                opts->input_path, opts->output_path);

    return sc_convert_via_temp_h5seurat(opts,
                                         sc_loom_to_h5seurat,
                                         sc_h5seurat_to_h5mu);
}

int sc_h5mu_to_loom(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5mu → loom via h5seurat)\n",
                opts->input_path, opts->output_path);

    return sc_convert_via_temp_h5seurat(opts,
                                         sc_h5mu_to_h5seurat,
                                         sc_h5seurat_to_loom);
}
