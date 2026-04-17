/*
 * sc_rcall.c — R .Call() entry points for fast h5ad reading
 *
 * Provides C-level readers that return R objects (SEXP) directly,
 * bypassing the R/hdf5r overhead for performance-critical paths.
 *
 * All functions use PROTECT/UNPROTECT and close HDF5 handles before returning.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <ctype.h>
#include "sc_convert.h"

/* ── Helpers ─────────────────────────────────────────────────────────────────── */

/*
 * Read a 1D vlen-string dataset into an R character vector (STRSXP).
 * Caller must PROTECT the result. Returns R_NilValue on failure.
 */
static SEXP read_string_dataset(hid_t loc, const char *name) {
    if (!sc_has_dataset(loc, name))
        return R_NilValue;

    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    if (dset < 0) return R_NilValue;

    hid_t space = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims != 1) {
        H5Sclose(space);
        H5Dclose(dset);
        return R_NilValue;
    }

    hsize_t n;
    H5Sget_simple_extent_dims(space, &n, NULL);
    H5Sclose(space);

    /* Use native file type for reading (avoids ASCII↔UTF8 conversion failures) */
    hid_t ftype = H5Dget_type(dset);
    hid_t memtype;
    if (H5Tis_variable_str(ftype) > 0) {
        memtype = H5Tcopy(ftype);
    } else {
        memtype = sc_create_vlen_str_type();
    }
    H5Tclose(ftype);

    char **raw = (char **)calloc((size_t)n, sizeof(char *));
    if (!raw) {
        H5Tclose(memtype);
        H5Dclose(dset);
        return R_NilValue;
    }

    if (H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw) < 0) {
        free(raw);
        H5Tclose(memtype);
        H5Dclose(dset);
        return R_NilValue;
    }

    SEXP result = PROTECT(allocVector(STRSXP, (R_xlen_t)n));
    for (hsize_t i = 0; i < n; i++) {
        if (raw[i]) {
            SET_STRING_ELT(result, (R_xlen_t)i,
                           mkCharCE(raw[i], CE_UTF8));
            H5free_memory(raw[i]);
        } else {
            SET_STRING_ELT(result, (R_xlen_t)i, NA_STRING);
        }
    }
    free(raw);

    H5Tclose(memtype);
    H5Dclose(dset);
    UNPROTECT(1);
    return result;
}

/*
 * Read a 1D numeric dataset as REALSXP.
 * Caller must PROTECT the result.
 */
static SEXP read_numeric_dataset(hid_t loc, const char *name, hsize_t expected_n) {
    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    if (dset < 0) return R_NilValue;

    hid_t space = H5Dget_space(dset);
    hsize_t n;
    H5Sget_simple_extent_dims(space, &n, NULL);
    H5Sclose(space);

    if (expected_n > 0 && n != expected_n) {
        H5Dclose(dset);
        return R_NilValue;
    }

    SEXP result = PROTECT(allocVector(REALSXP, (R_xlen_t)n));
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            REAL(result));
    H5Dclose(dset);
    UNPROTECT(1);
    return result;
}

/*
 * Read a 1D integer dataset as INTSXP.
 * Caller must PROTECT the result.
 */
static SEXP read_int_dataset(hid_t loc, const char *name, hsize_t expected_n) {
    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    if (dset < 0) return R_NilValue;

    hid_t space = H5Dget_space(dset);
    hsize_t n;
    H5Sget_simple_extent_dims(space, &n, NULL);
    H5Sclose(space);

    if (expected_n > 0 && n != expected_n) {
        H5Dclose(dset);
        return R_NilValue;
    }

    SEXP result = PROTECT(allocVector(INTSXP, (R_xlen_t)n));
    H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            INTEGER(result));
    H5Dclose(dset);
    UNPROTECT(1);
    return result;
}

/*
 * Sort row indices within each column of a CSC matrix.
 * dgCMatrix requires @i to be sorted within each @p segment.
 * Simultaneously sorts values in x to stay aligned with i.
 */
static void sort_csc_indices(int *i_arr, double *x_arr, int *p_arr, int ncol) {
    for (int col = 0; col < ncol; col++) {
        int start = p_arr[col];
        int end = p_arr[col + 1];
        int len = end - start;
        if (len <= 1) continue;

        /* Simple insertion sort — columns are typically short */
        for (int a = start + 1; a < end; a++) {
            int key_i = i_arr[a];
            double key_x = x_arr[a];
            int b = a - 1;
            while (b >= start && i_arr[b] > key_i) {
                i_arr[b + 1] = i_arr[b];
                x_arr[b + 1] = x_arr[b];
                b--;
            }
            i_arr[b + 1] = key_i;
            x_arr[b + 1] = key_x;
        }
    }
}

/* ── 1. C_read_h5ad_matrix ──────────────────────────────────────────────────── */

SEXP C_read_h5ad_matrix(SEXP path_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));

    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        Rf_warning("Cannot open HDF5 file: %s", path);
        return R_NilValue;
    }

    int nprotect = 0;

    /* Read obs/_index (cell barcodes) and var/_index (gene names) */
    SEXP colnames = R_NilValue;  /* cells = columns in genes x cells */
    SEXP rownames = R_NilValue;  /* genes = rows in genes x cells */

    if (sc_has_group(file, "obs")) {
        hid_t obs = H5Gopen2(file, "obs", H5P_DEFAULT);
        colnames = PROTECT(read_string_dataset(obs, "_index"));
        nprotect++;
        H5Gclose(obs);
    }
    if (sc_has_group(file, "var")) {
        hid_t var = H5Gopen2(file, "var", H5P_DEFAULT);
        rownames = PROTECT(read_string_dataset(var, "_index"));
        nprotect++;
        H5Gclose(var);
    }

    /* Determine if X is a group (sparse) or dataset (dense) */
    int x_is_group = sc_has_group(file, "X");
    int x_is_dataset = sc_has_dataset(file, "X");

    if (!x_is_group && !x_is_dataset) {
        Rf_warning("No X matrix found in %s", path);
        H5Fclose(file);
        UNPROTECT(nprotect);
        return R_NilValue;
    }

    /* ── Dense X ────────────────────────────────────────────────────────────── */
    if (x_is_dataset) {
        hid_t dset = H5Dopen2(file, "X", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        int ndims = H5Sget_simple_extent_ndims(space);

        if (ndims != 2) {
            H5Sclose(space);
            H5Dclose(dset);
            H5Fclose(file);
            UNPROTECT(nprotect);
            Rf_warning("X is a dataset but not 2D");
            return R_NilValue;
        }

        hsize_t dims[2];
        H5Sget_simple_extent_dims(space, dims, NULL);
        H5Sclose(space);

        /* h5ad dense X: [n_obs, n_vars] in C row-major (HDF5 storage)
         * R needs genes x cells in column-major
         * HDF5 reads row-major into flat buffer: buf[i*n_vars + j] = X[cell_i, gene_j]
         * R column-major matrix(nrow=n_genes, ncol=n_cells): mat[gene + n_genes*cell]
         * So we transpose: mat[j * n_obs + i] = buf[i * n_vars + j] */
        hsize_t n_obs = dims[0];
        hsize_t n_vars = dims[1];
        hsize_t total = n_obs * n_vars;

        double *buf = (double *)malloc(total * sizeof(double));
        if (!buf) {
            H5Dclose(dset);
            H5Fclose(file);
            UNPROTECT(nprotect);
            return R_NilValue;
        }
        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dclose(dset);

        /* Allocate R matrix (genes x cells) and transpose */
        SEXP mat = PROTECT(allocMatrix(REALSXP, (int)n_vars, (int)n_obs));
        nprotect++;
        double *mat_ptr = REAL(mat);

        for (hsize_t i = 0; i < n_obs; i++) {
            for (hsize_t j = 0; j < n_vars; j++) {
                mat_ptr[j + n_vars * i] = buf[i * n_vars + j];
            }
        }
        free(buf);

        /* Build result list */
        SEXP result = PROTECT(allocVector(VECSXP, 5));
        nprotect++;
        SEXP names = PROTECT(allocVector(STRSXP, 5));
        nprotect++;

        SET_VECTOR_ELT(result, 0, mat);
        SET_STRING_ELT(names, 0, mkChar("matrix"));

        SEXP nrow_s = PROTECT(ScalarInteger((int)n_vars));
        nprotect++;
        SET_VECTOR_ELT(result, 1, nrow_s);
        SET_STRING_ELT(names, 1, mkChar("nrow"));

        SEXP ncol_s = PROTECT(ScalarInteger((int)n_obs));
        nprotect++;
        SET_VECTOR_ELT(result, 2, ncol_s);
        SET_STRING_ELT(names, 2, mkChar("ncol"));

        SET_VECTOR_ELT(result, 3, (rownames != R_NilValue) ? rownames : R_NilValue);
        SET_STRING_ELT(names, 3, mkChar("rownames"));

        SET_VECTOR_ELT(result, 4, (colnames != R_NilValue) ? colnames : R_NilValue);
        SET_STRING_ELT(names, 4, mkChar("colnames"));

        /* Mark as dense */
        SEXP dense_flag = PROTECT(ScalarLogical(TRUE));
        nprotect++;
        setAttrib(result, install("dense"), dense_flag);

        setAttrib(result, R_NamesSymbol, names);
        H5Fclose(file);
        UNPROTECT(nprotect);
        return result;
    }

    /* ── Sparse X (group with data/indices/indptr) ──────────────────────────── */
    hid_t x_grp = H5Gopen2(file, "X", H5P_DEFAULT);
    if (x_grp < 0) {
        H5Fclose(file);
        UNPROTECT(nprotect);
        return R_NilValue;
    }

    /* Get encoding type: csr_matrix or csc_matrix */
    char enc[64] = "";
    sc_get_encoding_type(x_grp, enc, sizeof(enc));
    int is_csr = (strcmp(enc, "csr_matrix") == 0);
    /* If no encoding-type attr, assume CSR (standard h5ad) */
    if (enc[0] == '\0') is_csr = 1;

    /* Read shape attribute */
    int64_t shape_i64[2] = {0, 0};
    int has_shape = 0;
    if (H5Aexists(x_grp, "shape") > 0) {
        hid_t attr = H5Aopen(x_grp, "shape", H5P_DEFAULT);
        hid_t aspace = H5Aget_space(attr);
        hsize_t adims;
        H5Sget_simple_extent_dims(aspace, &adims, NULL);
        if (adims == 2) {
            hid_t atype = H5Aget_type(attr);
            if (H5Tget_size(atype) <= 4) {
                int32_t s32[2];
                H5Aread(attr, H5T_NATIVE_INT32, s32);
                shape_i64[0] = s32[0];
                shape_i64[1] = s32[1];
            } else {
                H5Aread(attr, H5T_NATIVE_INT64, shape_i64);
            }
            H5Tclose(atype);
            has_shape = 1;
        }
        H5Sclose(aspace);
        H5Aclose(attr);
    }

    /* Read nnz from data dataset */
    hsize_t nnz = 0;
    {
        hid_t data_dset = H5Dopen2(x_grp, "data", H5P_DEFAULT);
        if (data_dset < 0) {
            H5Gclose(x_grp);
            H5Fclose(file);
            UNPROTECT(nprotect);
            return R_NilValue;
        }
        hid_t dspace = H5Dget_space(data_dset);
        H5Sget_simple_extent_dims(dspace, &nnz, NULL);
        H5Sclose(dspace);
        H5Dclose(data_dset);
    }

    /* Read indptr length */
    hsize_t indptr_len = 0;
    {
        hid_t ptr_dset = H5Dopen2(x_grp, "indptr", H5P_DEFAULT);
        if (ptr_dset < 0) {
            H5Gclose(x_grp);
            H5Fclose(file);
            UNPROTECT(nprotect);
            return R_NilValue;
        }
        hid_t pspace = H5Dget_space(ptr_dset);
        H5Sget_simple_extent_dims(pspace, &indptr_len, NULL);
        H5Sclose(pspace);
        H5Dclose(ptr_dset);
    }

    /*
     * For dgCMatrix (CSC of genes x cells):
     *   @i   = row indices (gene indices), 0-based, length nnz
     *   @p   = column pointers (cell pointers), length ncol+1
     *   @x   = values, length nnz
     *   @Dim = c(nrow=n_genes, ncol=n_cells)
     *
     * h5ad CSR stores cells x genes:
     *   indptr[n_cells+1], indices[nnz] = gene indices, data[nnz]
     *   CSR of (cells x genes) == CSC of (genes x cells)
     *   So: indptr -> @p, indices -> @i directly (zero-copy reinterpretation)
     *
     * h5ad CSC stores cells x genes:
     *   indptr[n_genes+1], indices[nnz] = cell indices, data[nnz]
     *   This is CSC of (cells x genes). We need CSC of (genes x cells).
     *   That requires actual transposition.
     *   For CSC, treat as: indptr -> column ptrs over genes, indices -> cell row indices
     *   So this is already genes-as-columns. We need genes-as-rows.
     *   Actually, CSC of (cells x genes) means columns are genes, rows are cells.
     *   We want CSC of (genes x cells): columns are cells, rows are genes.
     *   CSC of (cells x genes) == CSR of (genes x cells).
     *   So indptr[n_genes+1], indices = cell indices. We'd need to convert CSR to CSC
     *   for the transposed matrix. That's expensive, so we do the simple approach:
     *   just read indices as-is and swap the dimension interpretation.
     *
     * Actually, for CSC h5ad (rare): the CSC of A (cells x genes) is equivalent to
     * CSR of A^T (genes x cells). To get dgCMatrix (CSC of genes x cells), we need
     * to actually perform a CSR->CSC conversion. This is the expensive path.
     * In practice, h5ad almost always uses CSR, so we optimize that path.
     */

    int n_genes, n_cells;

    if (is_csr) {
        /* CSR of (cells x genes): indptr has n_cells+1 entries */
        if (has_shape) {
            n_cells = (int)shape_i64[0];
            n_genes = (int)shape_i64[1];
        } else {
            n_cells = (int)(indptr_len - 1);
            n_genes = (rownames != R_NilValue) ? LENGTH(rownames) : 0;
        }

        /* Zero-copy reinterpretation: CSR(cells x genes) = CSC(genes x cells) */
        /* indptr -> @p (column pointers for cells) */
        /* indices -> @i (row indices = gene indices) */

        SEXP r_x = PROTECT(allocVector(REALSXP, (R_xlen_t)nnz));
        nprotect++;
        SEXP r_i = PROTECT(allocVector(INTSXP, (R_xlen_t)nnz));
        nprotect++;
        SEXP r_p = PROTECT(allocVector(INTSXP, (R_xlen_t)indptr_len));
        nprotect++;

        /* Read data -> x */
        {
            hid_t dset = H5Dopen2(x_grp, "data", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    REAL(r_x));
            H5Dclose(dset);
        }

        /* Read indices -> i */
        {
            hid_t dset = H5Dopen2(x_grp, "indices", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    INTEGER(r_i));
            H5Dclose(dset);
        }

        /* Read indptr -> p */
        {
            hid_t dset = H5Dopen2(x_grp, "indptr", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    INTEGER(r_p));
            H5Dclose(dset);
        }

        /* Sort row indices within each column (dgCMatrix requires sorted @i) */
        sort_csc_indices(INTEGER(r_i), REAL(r_x), INTEGER(r_p), n_cells);

        /* Build result list */
        SEXP result = PROTECT(allocVector(VECSXP, 7));
        nprotect++;
        SEXP names = PROTECT(allocVector(STRSXP, 7));
        nprotect++;

        SET_VECTOR_ELT(result, 0, r_i);
        SET_STRING_ELT(names, 0, mkChar("i"));

        SET_VECTOR_ELT(result, 1, r_p);
        SET_STRING_ELT(names, 1, mkChar("p"));

        SET_VECTOR_ELT(result, 2, r_x);
        SET_STRING_ELT(names, 2, mkChar("x"));

        SEXP nrow_s = PROTECT(ScalarInteger(n_genes));
        nprotect++;
        SET_VECTOR_ELT(result, 3, nrow_s);
        SET_STRING_ELT(names, 3, mkChar("nrow"));

        SEXP ncol_s = PROTECT(ScalarInteger(n_cells));
        nprotect++;
        SET_VECTOR_ELT(result, 4, ncol_s);
        SET_STRING_ELT(names, 4, mkChar("ncol"));

        SET_VECTOR_ELT(result, 5, (rownames != R_NilValue) ? rownames : R_NilValue);
        SET_STRING_ELT(names, 5, mkChar("rownames"));

        SET_VECTOR_ELT(result, 6, (colnames != R_NilValue) ? colnames : R_NilValue);
        SET_STRING_ELT(names, 6, mkChar("colnames"));

        setAttrib(result, R_NamesSymbol, names);

        H5Gclose(x_grp);
        H5Fclose(file);
        UNPROTECT(nprotect);
        return result;

    } else {
        /* CSC of (cells x genes): columns are genes, rows are cells
         * indptr has n_genes+1 entries, indices are cell indices
         * We need CSC of (genes x cells): columns are cells, rows are genes
         * This requires actual CSR-to-CSC conversion of the transposed matrix.
         */
        if (has_shape) {
            n_cells = (int)shape_i64[0];
            n_genes = (int)shape_i64[1];
        } else {
            /* indptr has n_genes+1 entries for CSC of cells x genes */
            n_genes = (int)(indptr_len - 1);
            n_cells = (colnames != R_NilValue) ? LENGTH(colnames) : 0;
        }

        /* Read the source arrays */
        int *src_indptr = (int *)malloc(indptr_len * sizeof(int));
        int *src_indices = (int *)malloc((size_t)nnz * sizeof(int));
        double *src_data = (double *)malloc((size_t)nnz * sizeof(double));

        if (!src_indptr || !src_indices || !src_data) {
            free(src_indptr);
            free(src_indices);
            free(src_data);
            H5Gclose(x_grp);
            H5Fclose(file);
            UNPROTECT(nprotect);
            return R_NilValue;
        }

        {
            hid_t dset = H5Dopen2(x_grp, "indptr", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_indptr);
            H5Dclose(dset);
        }
        {
            hid_t dset = H5Dopen2(x_grp, "indices", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_indices);
            H5Dclose(dset);
        }
        {
            hid_t dset = H5Dopen2(x_grp, "data", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_data);
            H5Dclose(dset);
        }

        H5Gclose(x_grp);
        H5Fclose(file);

        /* CSC of (cells x genes) → CSC of (genes x cells)
         * Source: col j (gene j) has entries src_indptr[j]..src_indptr[j+1]-1
         *         each entry has row index = cell index in src_indices
         * Target: col c (cell c) has entries with row index = gene index
         *
         * This is a standard CSC transpose: count entries per target column,
         * then scatter.
         */

        /* Count entries per cell (target column) */
        int *col_count = (int *)calloc((size_t)n_cells, sizeof(int));
        if (!col_count) {
            free(src_indptr);
            free(src_indices);
            free(src_data);
            UNPROTECT(nprotect);
            return R_NilValue;
        }
        for (hsize_t k = 0; k < nnz; k++) {
            int cell_idx = src_indices[k];
            if (cell_idx >= 0 && cell_idx < n_cells)
                col_count[cell_idx]++;
        }

        /* Build target indptr (p) */
        SEXP r_p = PROTECT(allocVector(INTSXP, (R_xlen_t)(n_cells + 1)));
        nprotect++;
        int *p = INTEGER(r_p);
        p[0] = 0;
        for (int c = 0; c < n_cells; c++)
            p[c + 1] = p[c] + col_count[c];

        /* Scatter into target arrays */
        SEXP r_i = PROTECT(allocVector(INTSXP, (R_xlen_t)nnz));
        nprotect++;
        SEXP r_x = PROTECT(allocVector(REALSXP, (R_xlen_t)nnz));
        nprotect++;
        int *dst_i = INTEGER(r_i);
        double *dst_x = REAL(r_x);

        /* Use col_count as running offset */
        memset(col_count, 0, (size_t)n_cells * sizeof(int));

        for (int gene = 0; gene < n_genes; gene++) {
            for (int k = src_indptr[gene]; k < src_indptr[gene + 1]; k++) {
                int cell = src_indices[k];
                if (cell < 0 || cell >= n_cells) continue;
                int pos = p[cell] + col_count[cell];
                dst_i[pos] = gene;
                dst_x[pos] = src_data[k];
                col_count[cell]++;
            }
        }

        free(col_count);
        free(src_indptr);
        free(src_indices);
        free(src_data);

        /* Sort row indices within each column */
        sort_csc_indices(dst_i, dst_x, p, n_cells);

        /* Build result list */
        SEXP result = PROTECT(allocVector(VECSXP, 7));
        nprotect++;
        SEXP rnames = PROTECT(allocVector(STRSXP, 7));
        nprotect++;

        SET_VECTOR_ELT(result, 0, r_i);
        SET_STRING_ELT(rnames, 0, mkChar("i"));

        SET_VECTOR_ELT(result, 1, r_p);
        SET_STRING_ELT(rnames, 1, mkChar("p"));

        SET_VECTOR_ELT(result, 2, r_x);
        SET_STRING_ELT(rnames, 2, mkChar("x"));

        SEXP nrow_s = PROTECT(ScalarInteger(n_genes));
        nprotect++;
        SET_VECTOR_ELT(result, 3, nrow_s);
        SET_STRING_ELT(rnames, 3, mkChar("nrow"));

        SEXP ncol_s = PROTECT(ScalarInteger(n_cells));
        nprotect++;
        SET_VECTOR_ELT(result, 4, ncol_s);
        SET_STRING_ELT(rnames, 4, mkChar("ncol"));

        SET_VECTOR_ELT(result, 5, (rownames != R_NilValue) ? rownames : R_NilValue);
        SET_STRING_ELT(rnames, 5, mkChar("rownames"));

        SET_VECTOR_ELT(result, 6, (colnames != R_NilValue) ? colnames : R_NilValue);
        SET_STRING_ELT(rnames, 6, mkChar("colnames"));

        setAttrib(result, R_NamesSymbol, rnames);
        UNPROTECT(nprotect);
        return result;
    }
}

/* ── 2. C_read_h5ad_obs ─────────────────────────────────────────────────────── */

SEXP C_read_h5ad_obs(SEXP path_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));

    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        Rf_warning("Cannot open HDF5 file: %s", path);
        return R_NilValue;
    }

    if (!sc_has_group(file, "obs")) {
        H5Fclose(file);
        return R_NilValue;
    }

    hid_t obs = H5Gopen2(file, "obs", H5P_DEFAULT);
    if (obs < 0) {
        H5Fclose(file);
        return R_NilValue;
    }

    int nprotect = 0;

    /* Get column-order to know which columns to read */
    char **col_order = NULL;
    hsize_t n_cols = 0;
    sc_get_str_array_attr(obs, "column-order", &col_order, &n_cols);

    /* Check for legacy __categories group */
    int has_legacy_cats = sc_has_group(obs, "__categories");
    hid_t legacy_cat_grp = -1;
    if (has_legacy_cats)
        legacy_cat_grp = H5Gopen2(obs, "__categories", H5P_DEFAULT);

    /* Total elements = _index + n_cols */
    R_xlen_t n_elements = (R_xlen_t)(1 + n_cols);
    SEXP result = PROTECT(allocVector(VECSXP, n_elements));
    nprotect++;
    SEXP names = PROTECT(allocVector(STRSXP, n_elements));
    nprotect++;

    /* Read _index first */
    SEXP idx = PROTECT(read_string_dataset(obs, "_index"));
    nprotect++;
    SET_VECTOR_ELT(result, 0, (idx != R_NilValue) ? idx : R_NilValue);
    SET_STRING_ELT(names, 0, mkChar("_index"));

    /* Read each column */
    for (hsize_t ci = 0; ci < n_cols; ci++) {
        const char *col_name = col_order[ci];
        R_xlen_t ri = (R_xlen_t)(ci + 1);
        SET_STRING_ELT(names, ri, mkCharCE(col_name, CE_UTF8));

        /* Check if it's a categorical (group with codes/categories) */
        if (sc_has_group(obs, col_name)) {
            hid_t cat_grp = H5Gopen2(obs, col_name, H5P_DEFAULT);

            if (sc_has_dataset(cat_grp, "codes") &&
                sc_has_dataset(cat_grp, "categories")) {

                /* Read categories (levels) */
                hid_t cat_dset = H5Dopen2(cat_grp, "categories", H5P_DEFAULT);
                hid_t cat_space = H5Dget_space(cat_dset);
                hsize_t n_cats;
                H5Sget_simple_extent_dims(cat_space, &n_cats, NULL);
                H5Sclose(cat_space);

                hid_t cat_ftype = H5Dget_type(cat_dset);
                hid_t str_type = (H5Tis_variable_str(cat_ftype) > 0) ?
                    H5Tcopy(cat_ftype) : sc_create_vlen_str_type();
                H5Tclose(cat_ftype);
                char **cats = (char **)calloc((size_t)n_cats, sizeof(char *));
                if (cats) {
                    H5Dread(cat_dset, str_type, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, cats);
                }
                H5Dclose(cat_dset);

                /* Read codes (0-based) */
                hid_t code_dset = H5Dopen2(cat_grp, "codes", H5P_DEFAULT);
                hid_t code_space = H5Dget_space(code_dset);
                hsize_t n_codes;
                H5Sget_simple_extent_dims(code_space, &n_codes, NULL);
                H5Sclose(code_space);

                /* Build R factor: 1-based integer with levels attribute */
                SEXP factor = PROTECT(allocVector(INTSXP, (R_xlen_t)n_codes));
                nprotect++;

                /* Read as int8 or int32 depending on source type */
                hid_t code_type = H5Dget_type(code_dset);
                size_t code_sz = H5Tget_size(code_type);
                H5Tclose(code_type);

                if (code_sz == 1) {
                    /* int8 codes — read into temp buffer, convert */
                    int8_t *tmp = (int8_t *)malloc((size_t)n_codes);
                    if (tmp) {
                        H5Dread(code_dset, H5T_NATIVE_INT8, H5S_ALL,
                                H5S_ALL, H5P_DEFAULT, tmp);
                        for (hsize_t k = 0; k < n_codes; k++) {
                            INTEGER(factor)[k] = (tmp[k] < 0) ? NA_INTEGER
                                                               : (int)(tmp[k] + 1);
                        }
                        free(tmp);
                    }
                } else {
                    int *codes_buf = (int *)malloc((size_t)n_codes * sizeof(int));
                    if (codes_buf) {
                        H5Dread(code_dset, H5T_NATIVE_INT, H5S_ALL,
                                H5S_ALL, H5P_DEFAULT, codes_buf);
                        for (hsize_t k = 0; k < n_codes; k++) {
                            INTEGER(factor)[k] = (codes_buf[k] < 0)
                                                     ? NA_INTEGER
                                                     : codes_buf[k] + 1;
                        }
                        free(codes_buf);
                    }
                }
                H5Dclose(code_dset);

                /* Set levels */
                if (cats) {
                    SEXP levels = PROTECT(allocVector(STRSXP, (R_xlen_t)n_cats));
                    nprotect++;
                    for (hsize_t k = 0; k < n_cats; k++) {
                        if (cats[k]) {
                            SET_STRING_ELT(levels, (R_xlen_t)k,
                                           mkCharCE(cats[k], CE_UTF8));
                            H5free_memory(cats[k]);
                        } else {
                            SET_STRING_ELT(levels, (R_xlen_t)k, NA_STRING);
                        }
                    }
                    free(cats);
                    setAttrib(factor, install("levels"), levels);
                }

                /* Set class = "factor" */
                SEXP cls = PROTECT(allocVector(STRSXP, 1));
                nprotect++;
                SET_STRING_ELT(cls, 0, mkChar("factor"));
                setAttrib(factor, R_ClassSymbol, cls);

                SET_VECTOR_ELT(result, ri, factor);
                H5Tclose(str_type);

            } else {
                /* Unknown group type — skip */
                SET_VECTOR_ELT(result, ri, R_NilValue);
            }
            H5Gclose(cat_grp);
            continue;
        }

        /* Plain dataset */
        if (!sc_has_dataset(obs, col_name)) {
            SET_VECTOR_ELT(result, ri, R_NilValue);
            continue;
        }

        hid_t dset = H5Dopen2(obs, col_name, H5P_DEFAULT);
        hid_t dtype = H5Dget_type(dset);
        H5T_class_t cls = H5Tget_class(dtype);

        hid_t dspace = H5Dget_space(dset);
        hsize_t dlen;
        H5Sget_simple_extent_dims(dspace, &dlen, NULL);
        H5Sclose(dspace);

        if (cls == H5T_STRING) {
            /* String column — use native file type to avoid charset conversion issues */
            hid_t memtype = (H5Tis_variable_str(dtype) > 0) ?
                H5Tcopy(dtype) : sc_create_vlen_str_type();
            char **strs = (char **)calloc((size_t)dlen, sizeof(char *));
            if (strs) {
                H5Dread(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);

                SEXP svec = PROTECT(allocVector(STRSXP, (R_xlen_t)dlen));
                nprotect++;
                for (hsize_t k = 0; k < dlen; k++) {
                    if (strs[k]) {
                        SET_STRING_ELT(svec, (R_xlen_t)k,
                                       mkCharCE(strs[k], CE_UTF8));
                        H5free_memory(strs[k]);
                    } else {
                        SET_STRING_ELT(svec, (R_xlen_t)k, NA_STRING);
                    }
                }
                free(strs);

                /* Check legacy __categories for this column */
                if (has_legacy_cats && legacy_cat_grp >= 0 &&
                    sc_has_dataset(legacy_cat_grp, col_name)) {
                    /* This is a legacy categorical stored as string dataset
                     * with categories in __categories/<col_name> */
                    SEXP cat_strs = PROTECT(
                        read_string_dataset(legacy_cat_grp, col_name));
                    nprotect++;
                    if (cat_strs != R_NilValue) {
                        /* Convert string column to factor */
                        R_xlen_t n_levs = XLENGTH(cat_strs);
                        SEXP factor = PROTECT(allocVector(INTSXP, (R_xlen_t)dlen));
                        nprotect++;

                        /* Map each string to its level index */
                        for (R_xlen_t k = 0; k < (R_xlen_t)dlen; k++) {
                            const char *val = CHAR(STRING_ELT(svec, k));
                            int found = 0;
                            for (R_xlen_t lv = 0; lv < n_levs; lv++) {
                                if (strcmp(val, CHAR(STRING_ELT(cat_strs, lv))) == 0) {
                                    INTEGER(factor)[k] = (int)(lv + 1);
                                    found = 1;
                                    break;
                                }
                            }
                            if (!found) INTEGER(factor)[k] = NA_INTEGER;
                        }
                        setAttrib(factor, install("levels"), cat_strs);
                        SEXP fcls = PROTECT(allocVector(STRSXP, 1));
                        nprotect++;
                        SET_STRING_ELT(fcls, 0, mkChar("factor"));
                        setAttrib(factor, R_ClassSymbol, fcls);
                        SET_VECTOR_ELT(result, ri, factor);
                    } else {
                        SET_VECTOR_ELT(result, ri, svec);
                    }
                } else {
                    SET_VECTOR_ELT(result, ri, svec);
                }
            }
            H5Tclose(memtype);

        } else if (cls == H5T_FLOAT) {
            /* Numeric column */
            SEXP nvec = PROTECT(allocVector(REALSXP, (R_xlen_t)dlen));
            nprotect++;
            H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, REAL(nvec));
            SET_VECTOR_ELT(result, ri, nvec);

        } else if (cls == H5T_INTEGER) {
            /* Could be integer or boolean — check size and signedness */
            size_t sz = H5Tget_size(dtype);
            H5T_sign_t sign = H5Tget_sign(dtype);

            /* Check for boolean encoding attribute */
            char enc_buf[64] = "";
            sc_get_str_attr(dset, "encoding-type", enc_buf, sizeof(enc_buf));
            int is_bool = (strcmp(enc_buf, "bool") == 0);

            /* Also treat 1-byte as potential boolean if no encoding attr */
            if (!is_bool && sz == 1 && sign == H5T_SGN_NONE) {
                /* Unsigned 1-byte — likely boolean */
                is_bool = 1;
            }

            if (is_bool) {
                /* Boolean column */
                SEXP bvec = PROTECT(allocVector(LGLSXP, (R_xlen_t)dlen));
                nprotect++;
                if (sz == 1) {
                    uint8_t *tmp = (uint8_t *)malloc((size_t)dlen);
                    if (tmp) {
                        H5Dread(dset, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, tmp);
                        for (hsize_t k = 0; k < dlen; k++)
                            LOGICAL(bvec)[k] = (int)tmp[k];
                        free(tmp);
                    }
                } else {
                    int *tmp = (int *)malloc((size_t)dlen * sizeof(int));
                    if (tmp) {
                        H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, tmp);
                        for (hsize_t k = 0; k < dlen; k++)
                            LOGICAL(bvec)[k] = tmp[k] ? 1 : 0;
                        free(tmp);
                    }
                }
                SET_VECTOR_ELT(result, ri, bvec);
            } else {
                /* Integer column */
                SEXP ivec = PROTECT(allocVector(INTSXP, (R_xlen_t)dlen));
                nprotect++;
                H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, INTEGER(ivec));
                SET_VECTOR_ELT(result, ri, ivec);
            }

        } else if (cls == H5T_ENUM) {
            /* HDF5 boolean encoded as enum (anndata >= 0.8) */
            SEXP bvec = PROTECT(allocVector(LGLSXP, (R_xlen_t)dlen));
            nprotect++;
            int8_t *tmp = (int8_t *)malloc((size_t)dlen);
            if (tmp) {
                H5Dread(dset, H5T_NATIVE_INT8, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, tmp);
                for (hsize_t k = 0; k < dlen; k++)
                    LOGICAL(bvec)[k] = (int)tmp[k];
                free(tmp);
            }
            SET_VECTOR_ELT(result, ri, bvec);

        } else {
            /* Unknown type — skip */
            SET_VECTOR_ELT(result, ri, R_NilValue);
        }

        H5Tclose(dtype);
        H5Dclose(dset);
    }

    setAttrib(result, R_NamesSymbol, names);

    if (legacy_cat_grp >= 0)
        H5Gclose(legacy_cat_grp);

    sc_free_str_array(col_order, n_cols);
    H5Gclose(obs);
    H5Fclose(file);
    UNPROTECT(nprotect);
    return result;
}

/* ── 3. C_read_h5ad_obsm ────────────────────────────────────────────────────── */

SEXP C_read_h5ad_obsm(SEXP path_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));

    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        Rf_warning("Cannot open HDF5 file: %s", path);
        return R_NilValue;
    }

    if (!sc_has_group(file, "obsm")) {
        H5Fclose(file);
        return R_NilValue;
    }

    hid_t obsm = H5Gopen2(file, "obsm", H5P_DEFAULT);
    if (obsm < 0) {
        H5Fclose(file);
        return R_NilValue;
    }

    int nprotect = 0;

    /* Read obs/_index for rownames */
    SEXP obs_index = R_NilValue;
    if (sc_has_group(file, "obs")) {
        hid_t obs_grp = H5Gopen2(file, "obs", H5P_DEFAULT);
        obs_index = PROTECT(read_string_dataset(obs_grp, "_index"));
        nprotect++;
        H5Gclose(obs_grp);
    }

    /* Count eligible members (skip "spatial") */
    hsize_t n_members;
    H5Gget_num_objs(obsm, &n_members);

    /* First pass: count non-spatial 2D datasets */
    int n_emb = 0;
    for (hsize_t i = 0; i < n_members; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(obsm, i, name, sizeof(name));
        if (strcmp(name, "spatial") == 0) continue;
        int obj_type = H5Gget_objtype_by_idx(obsm, i);
        if (obj_type == H5G_DATASET) n_emb++;
    }

    SEXP result = PROTECT(allocVector(VECSXP, n_emb));
    nprotect++;
    SEXP rnames = PROTECT(allocVector(STRSXP, n_emb));
    nprotect++;

    int ei = 0;
    for (hsize_t i = 0; i < n_members && ei < n_emb; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(obsm, i, name, sizeof(name));
        if (strcmp(name, "spatial") == 0) continue;
        int obj_type = H5Gget_objtype_by_idx(obsm, i);
        if (obj_type != H5G_DATASET) continue;

        hid_t dset = H5Dopen2(obsm, name, H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        int ndims = H5Sget_simple_extent_ndims(space);

        if (ndims != 2) {
            H5Sclose(space);
            H5Dclose(dset);
            SET_VECTOR_ELT(result, ei, R_NilValue);
            SET_STRING_ELT(rnames, ei, mkCharCE(name, CE_UTF8));
            ei++;
            continue;
        }

        hsize_t dims[2];
        H5Sget_simple_extent_dims(space, dims, NULL);
        H5Sclose(space);

        /* h5ad obsm: HDF5 row-major [n_cells, n_dims]
         * R matrix: column-major, nrow=n_cells, ncol=n_dims
         * Reading row-major data into R column-major requires transpose.
         *
         * HDF5 stores buf[cell * n_dims + dim].
         * R column-major: mat[cell + n_cells * dim].
         */
        hsize_t n_cells = dims[0];
        hsize_t n_dims = dims[1];
        hsize_t total = n_cells * n_dims;

        double *buf = (double *)malloc(total * sizeof(double));
        if (!buf) {
            H5Dclose(dset);
            SET_VECTOR_ELT(result, ei, R_NilValue);
            SET_STRING_ELT(rnames, ei, mkCharCE(name, CE_UTF8));
            ei++;
            continue;
        }

        H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dclose(dset);

        SEXP mat = PROTECT(allocMatrix(REALSXP, (int)n_cells, (int)n_dims));
        nprotect++;
        double *mat_ptr = REAL(mat);

        /* Transpose: buf[cell * n_dims + dim] -> mat[cell + n_cells * dim] */
        for (hsize_t c = 0; c < n_cells; c++) {
            for (hsize_t d = 0; d < n_dims; d++) {
                mat_ptr[c + n_cells * d] = buf[c * n_dims + d];
            }
        }
        free(buf);

        /* Set dimnames: rownames = obs/_index, colnames = key_1, key_2, ... */
        SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
        nprotect++;

        /* Rownames from obs/_index */
        SET_VECTOR_ELT(dimnames, 0,
                        (obs_index != R_NilValue) ? obs_index : R_NilValue);

        /* Colnames: strip X_ prefix if present, then name_1, name_2, ... */
        const char *key = name;
        if (strncmp(name, "X_", 2) == 0)
            key = name + 2;

        SEXP cnames = PROTECT(allocVector(STRSXP, (R_xlen_t)n_dims));
        nprotect++;
        for (hsize_t d = 0; d < n_dims; d++) {
            char colname[SC_MAX_NAME_LEN];
            snprintf(colname, sizeof(colname), "%s_%d", key, (int)(d + 1));
            SET_STRING_ELT(cnames, (R_xlen_t)d, mkChar(colname));
        }
        SET_VECTOR_ELT(dimnames, 1, cnames);
        setAttrib(mat, R_DimNamesSymbol, dimnames);

        SET_VECTOR_ELT(result, ei, mat);
        SET_STRING_ELT(rnames, ei, mkCharCE(name, CE_UTF8));
        ei++;
    }

    setAttrib(result, R_NamesSymbol, rnames);

    H5Gclose(obsm);
    H5Fclose(file);
    UNPROTECT(nprotect);
    return result;
}

/* ── 4. C_read_h5ad_obsp ────────────────────────────────────────────────────── */

SEXP C_read_h5ad_obsp(SEXP path_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));

    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        Rf_warning("Cannot open HDF5 file: %s", path);
        return R_NilValue;
    }

    if (!sc_has_group(file, "obsp")) {
        H5Fclose(file);
        return R_NilValue;
    }

    hid_t obsp = H5Gopen2(file, "obsp", H5P_DEFAULT);
    if (obsp < 0) {
        H5Fclose(file);
        return R_NilValue;
    }

    int nprotect = 0;

    hsize_t n_members;
    H5Gget_num_objs(obsp, &n_members);

    /* Count sparse graph groups */
    int n_graphs = 0;
    for (hsize_t i = 0; i < n_members; i++) {
        int obj_type = H5Gget_objtype_by_idx(obsp, i);
        if (obj_type == H5G_GROUP) n_graphs++;
    }

    SEXP result = PROTECT(allocVector(VECSXP, n_graphs));
    nprotect++;
    SEXP rnames = PROTECT(allocVector(STRSXP, n_graphs));
    nprotect++;

    int gi = 0;
    for (hsize_t i = 0; i < n_members && gi < n_graphs; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(obsp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(obsp, i);

        if (obj_type != H5G_GROUP) continue;

        SET_STRING_ELT(rnames, gi, mkCharCE(name, CE_UTF8));

        hid_t grp = H5Gopen2(obsp, name, H5P_DEFAULT);
        if (grp < 0 || !sc_has_dataset(grp, "data") ||
            !sc_has_dataset(grp, "indices") || !sc_has_dataset(grp, "indptr")) {
            if (grp >= 0) H5Gclose(grp);
            SET_VECTOR_ELT(result, gi, R_NilValue);
            gi++;
            continue;
        }

        /* Get encoding type to determine CSR vs CSC */
        char enc[64] = "";
        sc_get_encoding_type(grp, enc, sizeof(enc));
        int is_csc = (strcmp(enc, "csc_matrix") == 0);

        /* Read shape */
        int64_t shape[2] = {0, 0};
        if (H5Aexists(grp, "shape") > 0) {
            hid_t attr = H5Aopen(grp, "shape", H5P_DEFAULT);
            hid_t aspace = H5Aget_space(attr);
            hsize_t adims;
            H5Sget_simple_extent_dims(aspace, &adims, NULL);
            if (adims == 2) {
                hid_t atype = H5Aget_type(attr);
                if (H5Tget_size(atype) <= 4) {
                    int32_t s32[2];
                    H5Aread(attr, H5T_NATIVE_INT32, s32);
                    shape[0] = s32[0]; shape[1] = s32[1];
                } else {
                    H5Aread(attr, H5T_NATIVE_INT64, shape);
                }
                H5Tclose(atype);
            }
            H5Sclose(aspace);
            H5Aclose(attr);
        }

        /* Read nnz */
        hsize_t nnz = 0;
        {
            hid_t dset = H5Dopen2(grp, "data", H5P_DEFAULT);
            hid_t dspace = H5Dget_space(dset);
            H5Sget_simple_extent_dims(dspace, &nnz, NULL);
            H5Sclose(dspace);
            H5Dclose(dset);
        }

        /* Read indptr length */
        hsize_t indptr_len = 0;
        {
            hid_t dset = H5Dopen2(grp, "indptr", H5P_DEFAULT);
            hid_t dspace = H5Dget_space(dset);
            H5Sget_simple_extent_dims(dspace, &indptr_len, NULL);
            H5Sclose(dspace);
            H5Dclose(dset);
        }

        /* For obsp, graphs are typically square (n_cells x n_cells).
         * Return raw i/p/x/nrow/ncol and let R construct the dgCMatrix.
         * For CSC: indptr = column pointers, indices = row indices → direct.
         * For CSR: indptr = row pointers, indices = col indices →
         *   reinterpret as CSC of transposed.
         */
        int nrow_val, ncol_val;
        if (is_csc) {
            /* CSC: indptr has ncol+1 entries */
            nrow_val = (int)shape[0];
            ncol_val = (int)(indptr_len - 1);
        } else {
            /* CSR: reinterpret as CSC of transposed */
            nrow_val = (int)shape[1]; /* columns become rows */
            ncol_val = (int)(indptr_len - 1); /* rows become columns */
        }

        /* Read arrays */
        SEXP r_x = PROTECT(allocVector(REALSXP, (R_xlen_t)nnz));
        nprotect++;
        SEXP r_i = PROTECT(allocVector(INTSXP, (R_xlen_t)nnz));
        nprotect++;
        SEXP r_p = PROTECT(allocVector(INTSXP, (R_xlen_t)indptr_len));
        nprotect++;

        {
            hid_t dset = H5Dopen2(grp, "data", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, REAL(r_x));
            H5Dclose(dset);
        }
        {
            hid_t dset = H5Dopen2(grp, "indices", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, INTEGER(r_i));
            H5Dclose(dset);
        }
        {
            hid_t dset = H5Dopen2(grp, "indptr", H5P_DEFAULT);
            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, INTEGER(r_p));
            H5Dclose(dset);
        }

        H5Gclose(grp);

        /* Sort indices within columns for dgCMatrix compliance */
        sort_csc_indices(INTEGER(r_i), REAL(r_x), INTEGER(r_p),
                         (int)(indptr_len - 1));

        /* Build sub-list: i, p, x, nrow, ncol */
        SEXP entry = PROTECT(allocVector(VECSXP, 5));
        nprotect++;
        SEXP enames = PROTECT(allocVector(STRSXP, 5));
        nprotect++;

        SET_VECTOR_ELT(entry, 0, r_i);
        SET_STRING_ELT(enames, 0, mkChar("i"));
        SET_VECTOR_ELT(entry, 1, r_p);
        SET_STRING_ELT(enames, 1, mkChar("p"));
        SET_VECTOR_ELT(entry, 2, r_x);
        SET_STRING_ELT(enames, 2, mkChar("x"));

        SEXP nrow_s = PROTECT(ScalarInteger(nrow_val));
        nprotect++;
        SET_VECTOR_ELT(entry, 3, nrow_s);
        SET_STRING_ELT(enames, 3, mkChar("nrow"));

        SEXP ncol_s = PROTECT(ScalarInteger(ncol_val));
        nprotect++;
        SET_VECTOR_ELT(entry, 4, ncol_s);
        SET_STRING_ELT(enames, 4, mkChar("ncol"));

        setAttrib(entry, R_NamesSymbol, enames);
        SET_VECTOR_ELT(result, gi, entry);
        gi++;
    }

    setAttrib(result, R_NamesSymbol, rnames);

    H5Gclose(obsp);
    H5Fclose(file);
    UNPROTECT(nprotect);
    return result;
}

/* ── 5. C_read_h5ad ─────────────────────────────────────────────────────────── */

SEXP C_read_h5ad(SEXP path_sexp, SEXP components_sexp) {
    R_xlen_t n_comp = XLENGTH(components_sexp);
    int nprotect = 0;

    SEXP result = PROTECT(allocVector(VECSXP, n_comp));
    nprotect++;
    SEXP names = PROTECT(allocVector(STRSXP, n_comp));
    nprotect++;

    for (R_xlen_t i = 0; i < n_comp; i++) {
        const char *comp = CHAR(STRING_ELT(components_sexp, i));
        SET_STRING_ELT(names, i, mkChar(comp));

        SEXP val = R_NilValue;

        if (strcmp(comp, "X") == 0) {
            val = PROTECT(C_read_h5ad_matrix(path_sexp));
            nprotect++;
        } else if (strcmp(comp, "obs") == 0) {
            val = PROTECT(C_read_h5ad_obs(path_sexp));
            nprotect++;
        } else if (strcmp(comp, "var") == 0) {
            /* Read var the same way as obs — reuse obs reader on /var */
            /* But var has different structure. Read _index + column-order columns. */
            const char *path = CHAR(STRING_ELT(path_sexp, 0));
            hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
            if (file >= 0 && sc_has_group(file, "var")) {
                hid_t var_grp = H5Gopen2(file, "var", H5P_DEFAULT);

                /* Get column-order */
                char **col_order = NULL;
                hsize_t n_cols = 0;
                sc_get_str_array_attr(var_grp, "column-order", &col_order, &n_cols);

                R_xlen_t n_elements = (R_xlen_t)(1 + n_cols);
                val = PROTECT(allocVector(VECSXP, n_elements));
                nprotect++;
                SEXP vnames = PROTECT(allocVector(STRSXP, n_elements));
                nprotect++;

                /* Read _index */
                SEXP vidx = PROTECT(read_string_dataset(var_grp, "_index"));
                nprotect++;
                SET_VECTOR_ELT(val, 0, (vidx != R_NilValue) ? vidx : R_NilValue);
                SET_STRING_ELT(vnames, 0, mkChar("_index"));

                /* Read each column as string/numeric/integer */
                for (hsize_t ci = 0; ci < n_cols; ci++) {
                    const char *col_name = col_order[ci];
                    R_xlen_t ri = (R_xlen_t)(ci + 1);
                    SET_STRING_ELT(vnames, ri, mkCharCE(col_name, CE_UTF8));

                    if (sc_has_dataset(var_grp, col_name)) {
                        hid_t dset = H5Dopen2(var_grp, col_name, H5P_DEFAULT);
                        hid_t dtype = H5Dget_type(dset);
                        H5T_class_t dcls = H5Tget_class(dtype);
                        hid_t dspace = H5Dget_space(dset);
                        hsize_t dlen;
                        H5Sget_simple_extent_dims(dspace, &dlen, NULL);
                        H5Sclose(dspace);

                        if (dcls == H5T_STRING) {
                            hid_t mt = (H5Tis_variable_str(dtype) > 0) ?
                                H5Tcopy(dtype) : sc_create_vlen_str_type();
                            char **strs = (char **)calloc((size_t)dlen, sizeof(char *));
                            if (strs) {
                                H5Dread(dset, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);
                                SEXP sv = PROTECT(allocVector(STRSXP, (R_xlen_t)dlen));
                                nprotect++;
                                for (hsize_t k = 0; k < dlen; k++) {
                                    if (strs[k]) {
                                        SET_STRING_ELT(sv, (R_xlen_t)k,
                                                       mkCharCE(strs[k], CE_UTF8));
                                        H5free_memory(strs[k]);
                                    } else {
                                        SET_STRING_ELT(sv, (R_xlen_t)k, NA_STRING);
                                    }
                                }
                                free(strs);
                                SET_VECTOR_ELT(val, ri, sv);
                            }
                            H5Tclose(mt);
                        } else if (dcls == H5T_FLOAT) {
                            SEXP nv = PROTECT(allocVector(REALSXP, (R_xlen_t)dlen));
                            nprotect++;
                            H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                    H5P_DEFAULT, REAL(nv));
                            SET_VECTOR_ELT(val, ri, nv);
                        } else if (dcls == H5T_INTEGER) {
                            SEXP iv = PROTECT(allocVector(INTSXP, (R_xlen_t)dlen));
                            nprotect++;
                            H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                                    H5P_DEFAULT, INTEGER(iv));
                            SET_VECTOR_ELT(val, ri, iv);
                        } else {
                            SET_VECTOR_ELT(val, ri, R_NilValue);
                        }

                        H5Tclose(dtype);
                        H5Dclose(dset);
                    } else if (sc_has_group(var_grp, col_name)) {
                        /* Categorical in var (rare but possible) */
                        hid_t cat_grp = H5Gopen2(var_grp, col_name, H5P_DEFAULT);
                        if (sc_has_dataset(cat_grp, "codes") &&
                            sc_has_dataset(cat_grp, "categories")) {
                            /* Simplified factor construction */
                            SEXP cats_vec = PROTECT(
                                read_string_dataset(cat_grp, "categories"));
                            nprotect++;
                            SEXP codes_vec = PROTECT(
                                read_int_dataset(cat_grp, "codes", 0));
                            nprotect++;

                            if (cats_vec != R_NilValue && codes_vec != R_NilValue) {
                                R_xlen_t nc = XLENGTH(codes_vec);
                                SEXP factor = PROTECT(allocVector(INTSXP, nc));
                                nprotect++;
                                for (R_xlen_t k = 0; k < nc; k++) {
                                    int code = INTEGER(codes_vec)[k];
                                    INTEGER(factor)[k] = (code < 0) ? NA_INTEGER
                                                                    : code + 1;
                                }
                                setAttrib(factor, install("levels"), cats_vec);
                                SEXP fcls = PROTECT(allocVector(STRSXP, 1));
                                nprotect++;
                                SET_STRING_ELT(fcls, 0, mkChar("factor"));
                                setAttrib(factor, R_ClassSymbol, fcls);
                                SET_VECTOR_ELT(val, ri, factor);
                            }
                        }
                        H5Gclose(cat_grp);
                    } else {
                        SET_VECTOR_ELT(val, ri, R_NilValue);
                    }
                }

                setAttrib(val, R_NamesSymbol, vnames);
                sc_free_str_array(col_order, n_cols);
                H5Gclose(var_grp);
            }
            if (file >= 0) H5Fclose(file);

        } else if (strcmp(comp, "obsm") == 0) {
            val = PROTECT(C_read_h5ad_obsm(path_sexp));
            nprotect++;
        } else if (strcmp(comp, "obsp") == 0) {
            val = PROTECT(C_read_h5ad_obsp(path_sexp));
            nprotect++;
        }
        /* else: unknown component, leave as R_NilValue */

        SET_VECTOR_ELT(result, i, val);
    }

    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(nprotect);
    return result;
}

/* ── C_write_h5seurat ────────────────────────────────────────────────────────
 *
 * Write a complete h5seurat file directly from R objects, bypassing hdf5r.
 * Returns TRUE on success, FALSE on error.
 */

static hid_t vlen_str_type(void) {
    hid_t t = H5Tcopy(H5T_C_S1);
    H5Tset_size(t, H5T_VARIABLE);
    return t;
}

/* Write a 1D chunked+compressed dataset from a raw buffer */
static hid_t write_1d_dataset(hid_t loc, const char *name, hid_t dtype,
                               hsize_t n, const void *buf, int gzip)
{
    hsize_t chunk = n < (1 << 20) ? n : (1 << 20);
    if (chunk == 0) chunk = 1;

    hid_t space = H5Screate_simple(1, &n, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 1, &chunk);
    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);

    hid_t dset = H5Dcreate2(loc, name, dtype, space, H5P_DEFAULT, dcpl,
                             H5P_DEFAULT);
    if (dset >= 0 && buf)
        H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    H5Pclose(dcpl);
    H5Sclose(space);
    return dset; /* caller closes */
}

/* Write an STRSXP vector as a 1D vlen-string dataset */
static herr_t write_string_dataset(hid_t loc, const char *name,
                                    SEXP vec, int gzip)
{
    R_xlen_t n = XLENGTH(vec);
    hsize_t hn = (hsize_t)n;
    hsize_t chunk = hn < (1 << 20) ? hn : (1 << 20);
    if (chunk == 0) chunk = 1;

    const char **ptrs = (const char **)malloc((size_t)n * sizeof(char *));
    if (!ptrs) return -1;
    for (R_xlen_t i = 0; i < n; i++)
        ptrs[i] = CHAR(STRING_ELT(vec, i));

    hid_t strtype = vlen_str_type();
    hid_t space = H5Screate_simple(1, &hn, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 1, &chunk);
    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);

    hid_t dset = H5Dcreate2(loc, name, strtype, space, H5P_DEFAULT, dcpl,
                             H5P_DEFAULT);
    herr_t status = -1;
    if (dset >= 0) {
        status = H5Dwrite(dset, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs);
        H5Dclose(dset);
    }

    H5Pclose(dcpl);
    H5Sclose(space);
    H5Tclose(strtype);
    free(ptrs);
    return status;
}

/* Write C string array as a 1D vlen-string dataset */
static herr_t write_cstring_dataset(hid_t loc, const char *name,
                                     const char **strs, hsize_t n, int gzip)
{
    hsize_t chunk = n < (1 << 20) ? n : (1 << 20);
    if (chunk == 0) chunk = 1;

    hid_t strtype = vlen_str_type();
    hid_t space = H5Screate_simple(1, &n, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 1, &chunk);
    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);

    hid_t dset = H5Dcreate2(loc, name, strtype, space, H5P_DEFAULT, dcpl,
                             H5P_DEFAULT);
    herr_t status = -1;
    if (dset >= 0) {
        status = H5Dwrite(dset, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);
        H5Dclose(dset);
    }
    H5Pclose(dcpl);
    H5Sclose(space);
    H5Tclose(strtype);
    return status;
}

/* Write a string attribute using ASCII vlen type (matching hdf5r output) */
static void set_str_attr_ascii(hid_t loc, const char *attr_name,
                                const char *value)
{
    hid_t strtype = vlen_str_type();
    hid_t space = H5Screate(H5S_SCALAR);
    if (H5Aexists(loc, attr_name) > 0)
        H5Adelete(loc, attr_name);
    hid_t attr = H5Acreate2(loc, attr_name, strtype, space, H5P_DEFAULT,
                             H5P_DEFAULT);
    if (attr >= 0) {
        const char *ptr = value;
        H5Awrite(attr, strtype, &ptr);
        H5Aclose(attr);
    }
    H5Sclose(space);
    H5Tclose(strtype);
}

/* Write a string-array attribute using ASCII vlen type */
static void set_str_array_attr_ascii(hid_t loc, const char *attr_name,
                                      const char **values, hsize_t n)
{
    hid_t strtype = vlen_str_type();
    hid_t space = H5Screate_simple(1, &n, NULL);
    if (H5Aexists(loc, attr_name) > 0)
        H5Adelete(loc, attr_name);
    hid_t attr = H5Acreate2(loc, attr_name, strtype, space, H5P_DEFAULT,
                             H5P_DEFAULT);
    if (attr >= 0) {
        H5Awrite(attr, strtype, values);
        H5Aclose(attr);
    }
    H5Sclose(space);
    H5Tclose(strtype);
}

/* Write an integer scalar attribute */
static void set_int_attr(hid_t loc, const char *attr_name, int value) {
    hid_t space = H5Screate(H5S_SCALAR);
    if (H5Aexists(loc, attr_name) > 0)
        H5Adelete(loc, attr_name);
    hid_t attr = H5Acreate2(loc, attr_name, H5T_NATIVE_INT, space,
                             H5P_DEFAULT, H5P_DEFAULT);
    if (attr >= 0) {
        H5Awrite(attr, H5T_NATIVE_INT, &value);
        H5Aclose(attr);
    }
    H5Sclose(space);
}


SEXP C_write_h5seurat(SEXP path_sexp, SEXP mat_sexp, SEXP meta_sexp,
                       SEXP reductions_sexp, SEXP graphs_sexp,
                       SEXP assay_sexp, SEXP gzip_sexp)
{
    const char *path     = CHAR(STRING_ELT(path_sexp, 0));
    const char *assay    = CHAR(STRING_ELT(assay_sexp, 0));
    int gzip             = INTEGER(gzip_sexp)[0];
    int nprotect         = 0;

    /* Extract sparse matrix components */
    SEXP mat_i        = VECTOR_ELT(mat_sexp, 0);  /* row indices, 0-based */
    SEXP mat_p        = VECTOR_ELT(mat_sexp, 1);  /* col pointers */
    SEXP mat_x        = VECTOR_ELT(mat_sexp, 2);  /* values */
    SEXP mat_dim      = VECTOR_ELT(mat_sexp, 3);  /* c(n_genes, n_cells) */
    SEXP gene_names   = VECTOR_ELT(mat_sexp, 4);  /* rownames */
    SEXP cell_barcodes = VECTOR_ELT(mat_sexp, 5); /* colnames */

    int n_genes = INTEGER(mat_dim)[0];
    int n_cells = INTEGER(mat_dim)[1];
    R_xlen_t nnz = XLENGTH(mat_x);

    /* Check for optional variable.features */
    SEXP var_features = R_NilValue;
    SEXP mat_names = getAttrib(mat_sexp, R_NamesSymbol);
    if (mat_names != R_NilValue) {
        for (R_xlen_t k = 0; k < XLENGTH(mat_names); k++) {
            if (strcmp(CHAR(STRING_ELT(mat_names, k)), "variable.features") == 0) {
                var_features = VECTOR_ELT(mat_sexp, k);
                break;
            }
        }
    }

    /* Build lowercase assay key (e.g. "RNA" -> "rna_") */
    char assay_key[64];
    size_t alen = strlen(assay);
    if (alen > 58) alen = 58;
    for (size_t ci = 0; ci < alen; ci++)
        assay_key[ci] = (char)tolower((unsigned char)assay[ci]);
    assay_key[alen] = '_';
    assay_key[alen + 1] = '\0';

    /* ── Create file ─────────────────────────────────────────────────────────── */
    hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        Rf_warning("C_write_h5seurat: cannot create file %s", path);
        return ScalarLogical(FALSE);
    }

    /* Root attributes */
    set_str_attr_ascii(file, "version", "5.0.0");
    set_str_attr_ascii(file, "project", "scConvert");
    set_str_attr_ascii(file, "active.assay", assay);

    /* ── cell.names ──────────────────────────────────────────────────────────── */
    write_string_dataset(file, "cell.names", cell_barcodes, gzip);

    /* ── active.ident (factor: all cells in one level "scConvert") ───────────── */
    {
        hid_t ident_grp = H5Gcreate2(file, "active.ident", H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);
        /* values: 1-based integer array, all 1s */
        int *ident_vals = (int *)calloc((size_t)n_cells, sizeof(int));
        if (ident_vals) {
            for (int ci = 0; ci < n_cells; ci++) ident_vals[ci] = 1;
            hid_t dset = write_1d_dataset(ident_grp, "values",
                                           H5T_NATIVE_INT,
                                           (hsize_t)n_cells, ident_vals, gzip);
            if (dset >= 0) H5Dclose(dset);
            free(ident_vals);
        }
        /* levels: single string "scConvert" */
        const char *lvl = "scConvert";
        write_cstring_dataset(ident_grp, "levels", &lvl, 1, gzip);
        H5Gclose(ident_grp);
    }

    /* ── meta.data ───────────────────────────────────────────────────────────── */
    {
        hid_t meta_grp = H5Gcreate2(file, "meta.data", H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);

        /* _index = cell barcodes */
        write_string_dataset(meta_grp, "_index", cell_barcodes, gzip);

        int n_meta_cols = (int)XLENGTH(meta_sexp);
        SEXP meta_names = getAttrib(meta_sexp, R_NamesSymbol);

        /* Gather column names and logical column names */
        const char **col_names = (const char **)malloc(
            (size_t)n_meta_cols * sizeof(char *));
        const char **logical_names = (const char **)malloc(
            (size_t)n_meta_cols * sizeof(char *));
        int n_logical = 0;

        for (int ci = 0; ci < n_meta_cols; ci++) {
            col_names[ci] = CHAR(STRING_ELT(meta_names, ci));
            SEXP col = VECTOR_ELT(meta_sexp, ci);

            if (TYPEOF(col) == LGLSXP)
                logical_names[n_logical++] = col_names[ci];
        }

        /* meta.data group attributes */
        set_str_attr_ascii(meta_grp, "_index", "_index");
        if (n_meta_cols > 0)
            set_str_array_attr_ascii(meta_grp, "colnames",
                                      col_names, (hsize_t)n_meta_cols);
        if (n_logical > 0)
            set_str_array_attr_ascii(meta_grp, "logicals",
                                      logical_names, (hsize_t)n_logical);

        /* Write each column */
        for (int ci = 0; ci < n_meta_cols; ci++) {
            const char *cname = col_names[ci];
            SEXP col = VECTOR_ELT(meta_sexp, ci);

            switch (TYPEOF(col)) {
            case REALSXP: {
                hid_t dset = write_1d_dataset(meta_grp, cname,
                                               H5T_NATIVE_DOUBLE,
                                               (hsize_t)XLENGTH(col),
                                               REAL(col), gzip);
                if (dset >= 0) H5Dclose(dset);
                break;
            }
            case INTSXP: {
                SEXP levels = getAttrib(col, R_LevelsSymbol);
                if (levels != R_NilValue) {
                    /* Factor: write as subgroup with values + levels */
                    hid_t fac_grp = H5Gcreate2(meta_grp, cname, H5P_DEFAULT,
                                                 H5P_DEFAULT, H5P_DEFAULT);
                    hid_t dset = write_1d_dataset(fac_grp, "values",
                                                   H5T_NATIVE_INT,
                                                   (hsize_t)XLENGTH(col),
                                                   INTEGER(col), gzip);
                    if (dset >= 0) H5Dclose(dset);
                    write_string_dataset(fac_grp, "levels", levels, gzip);
                    H5Gclose(fac_grp);
                } else {
                    /* Plain integer */
                    hid_t dset = write_1d_dataset(meta_grp, cname,
                                                   H5T_NATIVE_INT,
                                                   (hsize_t)XLENGTH(col),
                                                   INTEGER(col), gzip);
                    if (dset >= 0) H5Dclose(dset);
                }
                break;
            }
            case LGLSXP: {
                /* Logical -> integer (TRUE=1, FALSE=0, NA=INT_MIN) */
                R_xlen_t len = XLENGTH(col);
                int *buf = (int *)malloc((size_t)len * sizeof(int));
                if (buf) {
                    int *src = LOGICAL(col);
                    for (R_xlen_t li = 0; li < len; li++)
                        buf[li] = (src[li] == NA_LOGICAL) ? INT_MIN : src[li];
                    hid_t dset = write_1d_dataset(meta_grp, cname,
                                                   H5T_NATIVE_INT,
                                                   (hsize_t)len, buf, gzip);
                    if (dset >= 0) H5Dclose(dset);
                    free(buf);
                }
                break;
            }
            case STRSXP:
                write_string_dataset(meta_grp, cname, col, gzip);
                break;
            default:
                /* Skip unsupported types */
                break;
            }
        }

        free(col_names);
        free(logical_names);
        H5Gclose(meta_grp);
    }

    /* ── assays/{assay}/ ─────────────────────────────────────────────────────── */
    {
        hid_t assays_grp = H5Gcreate2(file, "assays", H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
        hid_t assay_grp = H5Gcreate2(assays_grp, assay, H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);

        set_str_attr_ascii(assay_grp, "s4class", "Assay5");
        set_str_attr_ascii(assay_grp, "key", assay_key);

        /* features (gene names) */
        write_string_dataset(assay_grp, "features", gene_names, gzip);

        /* variable.features (if present) */
        if (var_features != R_NilValue && TYPEOF(var_features) == STRSXP)
            write_string_dataset(assay_grp, "variable.features",
                                 var_features, gzip);

        /* layers/counts/ — CSC sparse: data, indices, indptr */
        hid_t layers_grp = H5Gcreate2(assay_grp, "layers", H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
        hid_t counts_grp = H5Gcreate2(layers_grp, "counts", H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);

        /* data (non-zero values) */
        {
            hid_t dset = write_1d_dataset(counts_grp, "data",
                                           H5T_NATIVE_DOUBLE,
                                           (hsize_t)nnz, REAL(mat_x), gzip);
            if (dset >= 0) H5Dclose(dset);
        }
        /* indices (row indices, int32) */
        {
            hid_t dset = write_1d_dataset(counts_grp, "indices",
                                           H5T_NATIVE_INT,
                                           (hsize_t)XLENGTH(mat_i),
                                           INTEGER(mat_i), gzip);
            if (dset >= 0) H5Dclose(dset);
        }
        /* indptr (column pointers, int32) */
        {
            hid_t dset = write_1d_dataset(counts_grp, "indptr",
                                           H5T_NATIVE_INT,
                                           (hsize_t)XLENGTH(mat_p),
                                           INTEGER(mat_p), gzip);
            if (dset >= 0) H5Dclose(dset);
        }

        /* dims attribute on counts group: c(n_genes, n_cells) */
        {
            int dims[2] = { n_genes, n_cells };
            hsize_t two = 2;
            hid_t space = H5Screate_simple(1, &two, NULL);
            hid_t attr = H5Acreate2(counts_grp, "dims", H5T_NATIVE_INT,
                                     space, H5P_DEFAULT, H5P_DEFAULT);
            if (attr >= 0) {
                H5Awrite(attr, H5T_NATIVE_INT, dims);
                H5Aclose(attr);
            }
            H5Sclose(space);
        }

        H5Gclose(counts_grp);
        H5Gclose(layers_grp);
        H5Gclose(assay_grp);
        H5Gclose(assays_grp);
    }

    /* ── reductions/ ─────────────────────────────────────────────────────────── */
    {
        hid_t red_grp = H5Gcreate2(file, "reductions", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
        int n_reductions = (int)XLENGTH(reductions_sexp);
        SEXP red_names = getAttrib(reductions_sexp, R_NamesSymbol);

        for (int ri = 0; ri < n_reductions; ri++) {
            const char *red_name = CHAR(STRING_ELT(red_names, ri));
            SEXP embedding = VECTOR_ELT(reductions_sexp, ri);
            SEXP embed_dim = getAttrib(embedding, R_DimSymbol);
            int embed_nrow = INTEGER(embed_dim)[0]; /* n_cells */
            int embed_ncol = INTEGER(embed_dim)[1]; /* n_dims */

            hid_t this_red = H5Gcreate2(red_grp, red_name, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);

            set_str_attr_ascii(this_red, "active.assay", assay);
            set_int_attr(this_red, "global", 0);

            /* key attribute from the matrix */
            SEXP key_attr = getAttrib(embedding, install("key"));
            if (key_attr != R_NilValue && TYPEOF(key_attr) == STRSXP)
                set_str_attr_ascii(this_red, "key",
                                    CHAR(STRING_ELT(key_attr, 0)));

            /* cell.embeddings: R column-major [n_cells, n_dims] buffer
             * maps to HDF5 row-major [n_dims, n_cells] — swap dims so
             * hdf5r reads it back correctly as [n_cells, n_dims] in R */
            {
                hsize_t dims2[2] = { (hsize_t)embed_ncol,
                                      (hsize_t)embed_nrow };
                hsize_t chunks2[2];
                chunks2[0] = dims2[0];  /* n_dims — always small */
                chunks2[1] = dims2[1] < (1 << 16) ? dims2[1] : (1 << 16);

                hid_t space = H5Screate_simple(2, dims2, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl, 2, chunks2);
                if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);

                hid_t dset = H5Dcreate2(this_red, "cell.embeddings",
                                          H5T_NATIVE_DOUBLE, space,
                                          H5P_DEFAULT, dcpl, H5P_DEFAULT);
                if (dset >= 0) {
                    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, REAL(embedding));
                    H5Dclose(dset);
                }
                H5Pclose(dcpl);
                H5Sclose(space);
            }

            H5Gclose(this_red);
        }

        H5Gclose(red_grp);
    }

    /* ── graphs/ ─────────────────────────────────────────────────────────────── */
    {
        hid_t graphs_grp = H5Gcreate2(file, "graphs", H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
        int n_graphs = (int)XLENGTH(graphs_sexp);
        SEXP graph_names = getAttrib(graphs_sexp, R_NamesSymbol);

        for (int gi = 0; gi < n_graphs; gi++) {
            const char *graph_name = CHAR(STRING_ELT(graph_names, gi));
            SEXP graph = VECTOR_ELT(graphs_sexp, gi);

            SEXP g_i   = VECTOR_ELT(graph, 0);
            SEXP g_p   = VECTOR_ELT(graph, 1);
            SEXP g_x   = VECTOR_ELT(graph, 2);
            SEXP g_dim = VECTOR_ELT(graph, 3);

            hid_t this_graph = H5Gcreate2(graphs_grp, graph_name, H5P_DEFAULT,
                                            H5P_DEFAULT, H5P_DEFAULT);

            /* data */
            {
                hid_t dset = write_1d_dataset(this_graph, "data",
                                               H5T_NATIVE_DOUBLE,
                                               (hsize_t)XLENGTH(g_x),
                                               REAL(g_x), gzip);
                if (dset >= 0) H5Dclose(dset);
            }
            /* indices */
            {
                hid_t dset = write_1d_dataset(this_graph, "indices",
                                               H5T_NATIVE_INT,
                                               (hsize_t)XLENGTH(g_i),
                                               INTEGER(g_i), gzip);
                if (dset >= 0) H5Dclose(dset);
            }
            /* indptr */
            {
                hid_t dset = write_1d_dataset(this_graph, "indptr",
                                               H5T_NATIVE_INT,
                                               (hsize_t)XLENGTH(g_p),
                                               INTEGER(g_p), gzip);
                if (dset >= 0) H5Dclose(dset);
            }

            /* dims attribute */
            {
                int gdims[2] = { INTEGER(g_dim)[0], INTEGER(g_dim)[1] };
                hsize_t two = 2;
                hid_t space = H5Screate_simple(1, &two, NULL);
                hid_t attr = H5Acreate2(this_graph, "dims", H5T_NATIVE_INT,
                                         space, H5P_DEFAULT, H5P_DEFAULT);
                if (attr >= 0) {
                    H5Awrite(attr, H5T_NATIVE_INT, gdims);
                    H5Aclose(attr);
                }
                H5Sclose(space);
            }

            /* assay.used attribute */
            set_str_attr_ascii(this_graph, "assay.used", assay);

            H5Gclose(this_graph);
        }

        H5Gclose(graphs_grp);
    }

    /* ── Empty groups: neighbors, commands, misc, tools, images ──────────────── */
    {
        const char *empty_groups[] = {
            "neighbors", "commands", "misc", "tools", "images"
        };
        for (int ei = 0; ei < 5; ei++) {
            hid_t grp = H5Gcreate2(file, empty_groups[ei], H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
            if (grp >= 0) H5Gclose(grp);
        }
    }

    H5Fclose(file);
    (void)nprotect; /* no R allocations needed in this function */
    return ScalarLogical(TRUE);
}

/* ── 5. C_write_h5ad ─────────────────────────────────────────────────────────── */

/*
 * Write an h5ad file directly from R SEXP objects.
 *
 * Key insight: R's dgCMatrix stores genes x cells in CSC format.
 * h5ad stores cells x genes in CSR format.
 * CSC of (genes x cells) == CSR of (cells x genes) — zero-copy reinterpretation.
 * Just write @p as indptr, @i as indices, @x as data, shape = [n_cells, n_genes].
 */

/* Helper: set encoding-type + encoding-version string attrs on a string dataset */
static void set_string_encoding_attrs(hid_t loc, const char *dset_name) {
    hid_t dset = H5Dopen2(loc, dset_name, H5P_DEFAULT);
    if (dset >= 0) {
        set_str_attr_ascii(dset, "encoding-type", "string-array");
        set_str_attr_ascii(dset, "encoding-version", "0.2.0");
        H5Dclose(dset);
    }
}

/* Helper: write a boolean attribute (H5T_NATIVE_HBOOL) */
static void set_bool_attr(hid_t loc, const char *attr_name, int value) {
    hid_t space = H5Screate(H5S_SCALAR);
    if (H5Aexists(loc, attr_name) > 0)
        H5Adelete(loc, attr_name);
    hid_t attr = H5Acreate2(loc, attr_name, H5T_NATIVE_HBOOL, space,
                             H5P_DEFAULT, H5P_DEFAULT);
    if (attr >= 0) {
        hbool_t bval = (hbool_t)value;
        H5Awrite(attr, H5T_NATIVE_HBOOL, &bval);
        H5Aclose(attr);
    }
    H5Sclose(space);
}

SEXP C_write_h5ad(SEXP path_sexp, SEXP mat_sexp, SEXP meta_sexp,
                   SEXP reductions_sexp, SEXP graphs_sexp,
                   SEXP assay_sexp, SEXP gzip_sexp)
{
    const char *path  = CHAR(STRING_ELT(path_sexp, 0));
    int gzip          = INTEGER(gzip_sexp)[0];

    /* Extract sparse matrix components from mat_sexp list */
    SEXP mat_i         = VECTOR_ELT(mat_sexp, 0);  /* row indices (gene), 0-based */
    SEXP mat_p         = VECTOR_ELT(mat_sexp, 1);  /* col pointers (cell) */
    SEXP mat_x         = VECTOR_ELT(mat_sexp, 2);  /* values */
    SEXP mat_dim       = VECTOR_ELT(mat_sexp, 3);  /* c(n_genes, n_cells) */
    SEXP gene_names    = VECTOR_ELT(mat_sexp, 4);  /* rownames */
    SEXP cell_barcodes = VECTOR_ELT(mat_sexp, 5);  /* colnames */

    int n_genes = INTEGER(mat_dim)[0];
    int n_cells = INTEGER(mat_dim)[1];
    R_xlen_t nnz = XLENGTH(mat_x);

    /* ── Create file ─────────────────────────────────────────────────────────── */
    hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        Rf_warning("C_write_h5ad: cannot create file %s", path);
        return ScalarLogical(FALSE);
    }

    /* Root attributes */
    set_str_attr_ascii(file, "encoding-type", "anndata");
    set_str_attr_ascii(file, "encoding-version", "0.1.0");

    /* ── /X (CSR sparse group) ───────────────────────────────────────────────── */
    {
        hid_t x_grp = H5Gcreate2(file, "X", H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);

        /* data (non-zero values) */
        {
            hid_t dset = write_1d_dataset(x_grp, "data", H5T_NATIVE_DOUBLE,
                                           (hsize_t)nnz, REAL(mat_x), gzip);
            if (dset >= 0) H5Dclose(dset);
        }
        /* indices (gene indices — same as dgCMatrix @i) */
        {
            hid_t dset = write_1d_dataset(x_grp, "indices", H5T_NATIVE_INT,
                                           (hsize_t)XLENGTH(mat_i),
                                           INTEGER(mat_i), gzip);
            if (dset >= 0) H5Dclose(dset);
        }
        /* indptr (cell pointers — same as dgCMatrix @p) */
        {
            hid_t dset = write_1d_dataset(x_grp, "indptr", H5T_NATIVE_INT,
                                           (hsize_t)XLENGTH(mat_p),
                                           INTEGER(mat_p), gzip);
            if (dset >= 0) H5Dclose(dset);
        }

        /* Attributes: encoding-type, encoding-version, shape */
        set_str_attr_ascii(x_grp, "encoding-type", "csr_matrix");
        set_str_attr_ascii(x_grp, "encoding-version", "0.1.0");
        {
            int shape[2] = { n_cells, n_genes };
            hsize_t two = 2;
            hid_t space = H5Screate_simple(1, &two, NULL);
            hid_t attr = H5Acreate2(x_grp, "shape", H5T_NATIVE_INT,
                                     space, H5P_DEFAULT, H5P_DEFAULT);
            if (attr >= 0) {
                H5Awrite(attr, H5T_NATIVE_INT, shape);
                H5Aclose(attr);
            }
            H5Sclose(space);
        }

        H5Gclose(x_grp);
    }

    /* ── /obs (dataframe group) ──────────────────────────────────────────────── */
    {
        hid_t obs_grp = H5Gcreate2(file, "obs", H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);

        /* _index = cell barcodes */
        write_string_dataset(obs_grp, "_index", cell_barcodes, gzip);
        set_string_encoding_attrs(obs_grp, "_index");

        int n_meta_cols = (int)XLENGTH(meta_sexp);
        SEXP meta_names = getAttrib(meta_sexp, R_NamesSymbol);

        /* Gather column names for column-order attribute */
        const char **col_names = NULL;
        if (n_meta_cols > 0)
            col_names = (const char **)malloc(
                (size_t)n_meta_cols * sizeof(char *));

        for (int ci = 0; ci < n_meta_cols; ci++)
            col_names[ci] = CHAR(STRING_ELT(meta_names, ci));

        /* obs group attributes */
        set_str_attr_ascii(obs_grp, "encoding-type", "dataframe");
        set_str_attr_ascii(obs_grp, "encoding-version", "0.2.0");
        set_str_attr_ascii(obs_grp, "_index", "_index");
        if (n_meta_cols > 0) {
            set_str_array_attr_ascii(obs_grp, "column-order",
                                      col_names, (hsize_t)n_meta_cols);
        } else {
            /* Empty column-order: 0-length string array */
            hsize_t zero = 0;
            hid_t strtype = vlen_str_type();
            hid_t space = H5Screate_simple(1, &zero, NULL);
            if (H5Aexists(obs_grp, "column-order") > 0)
                H5Adelete(obs_grp, "column-order");
            hid_t attr = H5Acreate2(obs_grp, "column-order", strtype, space,
                                     H5P_DEFAULT, H5P_DEFAULT);
            if (attr >= 0) H5Aclose(attr);
            H5Sclose(space);
            H5Tclose(strtype);
        }

        /* Write each metadata column */
        for (int ci = 0; ci < n_meta_cols; ci++) {
            const char *cname = col_names[ci];
            SEXP col = VECTOR_ELT(meta_sexp, ci);

            switch (TYPEOF(col)) {
            case REALSXP: {
                hid_t dset = write_1d_dataset(obs_grp, cname,
                                               H5T_NATIVE_DOUBLE,
                                               (hsize_t)XLENGTH(col),
                                               REAL(col), gzip);
                if (dset >= 0) {
                    set_str_attr_ascii(dset, "encoding-type", "array");
                    set_str_attr_ascii(dset, "encoding-version", "0.2.0");
                    H5Dclose(dset);
                }
                break;
            }
            case INTSXP: {
                SEXP levels = getAttrib(col, R_LevelsSymbol);
                if (levels != R_NilValue) {
                    /* Factor -> h5ad categorical (0-based codes) */
                    hid_t cat_grp = H5Gcreate2(obs_grp, cname, H5P_DEFAULT,
                                                 H5P_DEFAULT, H5P_DEFAULT);

                    /* categories (string dataset) */
                    write_string_dataset(cat_grp, "categories", levels, gzip);
                    set_string_encoding_attrs(cat_grp, "categories");

                    /* codes: subtract 1 from R's 1-based to get 0-based */
                    R_xlen_t len = XLENGTH(col);
                    int use_int8 = (XLENGTH(levels) <= 127);

                    if (use_int8) {
                        int8_t *codes = (int8_t *)malloc((size_t)len);
                        if (codes) {
                            int *src = INTEGER(col);
                            for (R_xlen_t k = 0; k < len; k++)
                                codes[k] = (src[k] == NA_INTEGER)
                                    ? (int8_t)-1 : (int8_t)(src[k] - 1);
                            hid_t dset = write_1d_dataset(cat_grp, "codes",
                                                           H5T_NATIVE_INT8,
                                                           (hsize_t)len,
                                                           codes, gzip);
                            if (dset >= 0) H5Dclose(dset);
                            free(codes);
                        }
                    } else {
                        int *codes = (int *)malloc((size_t)len * sizeof(int));
                        if (codes) {
                            int *src = INTEGER(col);
                            for (R_xlen_t k = 0; k < len; k++)
                                codes[k] = (src[k] == NA_INTEGER)
                                    ? -1 : (src[k] - 1);
                            hid_t dset = write_1d_dataset(cat_grp, "codes",
                                                           H5T_NATIVE_INT,
                                                           (hsize_t)len,
                                                           codes, gzip);
                            if (dset >= 0) H5Dclose(dset);
                            free(codes);
                        }
                    }

                    /* Categorical attributes */
                    set_str_attr_ascii(cat_grp, "encoding-type", "categorical");
                    set_str_attr_ascii(cat_grp, "encoding-version", "0.2.0");
                    set_bool_attr(cat_grp, "ordered", 0);

                    H5Gclose(cat_grp);
                } else {
                    /* Plain integer */
                    hid_t dset = write_1d_dataset(obs_grp, cname,
                                                   H5T_NATIVE_INT,
                                                   (hsize_t)XLENGTH(col),
                                                   INTEGER(col), gzip);
                    if (dset >= 0) {
                        set_str_attr_ascii(dset, "encoding-type", "array");
                        set_str_attr_ascii(dset, "encoding-version", "0.2.0");
                        H5Dclose(dset);
                    }
                }
                break;
            }
            case LGLSXP: {
                /* Logical -> boolean (stored as int8) */
                R_xlen_t len = XLENGTH(col);
                int8_t *buf = (int8_t *)malloc((size_t)len);
                if (buf) {
                    int *src = LOGICAL(col);
                    for (R_xlen_t li = 0; li < len; li++)
                        buf[li] = (src[li] == NA_LOGICAL)
                            ? (int8_t)0 : (int8_t)src[li];
                    hid_t dset = write_1d_dataset(obs_grp, cname,
                                                   H5T_NATIVE_INT8,
                                                   (hsize_t)len, buf, gzip);
                    if (dset >= 0) {
                        set_str_attr_ascii(dset, "encoding-type", "array");
                        set_str_attr_ascii(dset, "encoding-version", "0.2.0");
                        H5Dclose(dset);
                    }
                    free(buf);
                }
                break;
            }
            case STRSXP:
                write_string_dataset(obs_grp, cname, col, gzip);
                set_string_encoding_attrs(obs_grp, cname);
                break;
            default:
                break;
            }
        }

        free(col_names);
        H5Gclose(obs_grp);
    }

    /* ── /var (dataframe group) ──────────────────────────────────────────────── */
    {
        hid_t var_grp = H5Gcreate2(file, "var", H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);

        write_string_dataset(var_grp, "_index", gene_names, gzip);
        set_string_encoding_attrs(var_grp, "_index");

        set_str_attr_ascii(var_grp, "encoding-type", "dataframe");
        set_str_attr_ascii(var_grp, "encoding-version", "0.2.0");
        set_str_attr_ascii(var_grp, "_index", "_index");

        /* Empty column-order */
        {
            hsize_t zero = 0;
            hid_t strtype = vlen_str_type();
            hid_t space = H5Screate_simple(1, &zero, NULL);
            hid_t attr = H5Acreate2(var_grp, "column-order", strtype, space,
                                     H5P_DEFAULT, H5P_DEFAULT);
            if (attr >= 0) H5Aclose(attr);
            H5Sclose(space);
            H5Tclose(strtype);
        }

        H5Gclose(var_grp);
    }

    /* ── /obsm (embeddings) ──────────────────────────────────────────────────── */
    {
        hid_t obsm_grp = H5Gcreate2(file, "obsm", H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);
        set_str_attr_ascii(obsm_grp, "encoding-type", "dict");
        set_str_attr_ascii(obsm_grp, "encoding-version", "0.1.0");
        int n_reductions = (int)XLENGTH(reductions_sexp);
        SEXP red_names = getAttrib(reductions_sexp, R_NamesSymbol);

        for (int ri = 0; ri < n_reductions; ri++) {
            const char *red_name = CHAR(STRING_ELT(red_names, ri));
            SEXP embedding = VECTOR_ELT(reductions_sexp, ri);
            SEXP embed_dim = getAttrib(embedding, R_DimSymbol);
            int embed_nrow = INTEGER(embed_dim)[0]; /* n_cells */
            int embed_ncol = INTEGER(embed_dim)[1]; /* n_dims */

            /* Build h5ad key: X_{name} */
            char obsm_key[128];
            snprintf(obsm_key, sizeof(obsm_key), "X_%s", red_name);

            /* R stores column-major [n_cells, n_dims].  anndata expects
             * C row-major [n_cells, n_dims].  Transpose the buffer so
             * that HDF5 dims are [n_cells, n_dims] and values are in
             * row-major order. */
            hsize_t dims2[2] = { (hsize_t)embed_nrow, (hsize_t)embed_ncol };
            hsize_t chunks2[2];
            chunks2[0] = dims2[0] < (1 << 16) ? dims2[0] : (1 << 16);
            chunks2[1] = dims2[1];

            /* Column-major → row-major transpose */
            R_xlen_t total = (R_xlen_t)embed_nrow * embed_ncol;
            double *tbuf = (double *)R_alloc(total, sizeof(double));
            const double *src = REAL(embedding);
            for (int c = 0; c < embed_ncol; c++) {
                for (int r = 0; r < embed_nrow; r++) {
                    tbuf[r * embed_ncol + c] = src[c * embed_nrow + r];
                }
            }

            hid_t space = H5Screate_simple(2, dims2, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(dcpl, 2, chunks2);
            if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);

            hid_t dset = H5Dcreate2(obsm_grp, obsm_key, H5T_NATIVE_DOUBLE,
                                      space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
            if (dset >= 0) {
                H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                         H5P_DEFAULT, tbuf);
                set_str_attr_ascii(dset, "encoding-type", "array");
                set_str_attr_ascii(dset, "encoding-version", "0.2.0");
                H5Dclose(dset);
            }
            H5Pclose(dcpl);
            H5Sclose(space);
        }

        H5Gclose(obsm_grp);
    }

    /* ── /obsp (graphs) ──────────────────────────────────────────────────────── */
    {
        hid_t obsp_grp = H5Gcreate2(file, "obsp", H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);
        set_str_attr_ascii(obsp_grp, "encoding-type", "dict");
        set_str_attr_ascii(obsp_grp, "encoding-version", "0.1.0");
        int n_graphs = (int)XLENGTH(graphs_sexp);
        SEXP graph_names = getAttrib(graphs_sexp, R_NamesSymbol);

        for (int gi = 0; gi < n_graphs; gi++) {
            const char *graph_name = CHAR(STRING_ELT(graph_names, gi));
            SEXP graph = VECTOR_ELT(graphs_sexp, gi);

            SEXP g_i   = VECTOR_ELT(graph, 0);
            SEXP g_p   = VECTOR_ELT(graph, 1);
            SEXP g_x   = VECTOR_ELT(graph, 2);
            SEXP g_dim = VECTOR_ELT(graph, 3);

            hid_t this_graph = H5Gcreate2(obsp_grp, graph_name, H5P_DEFAULT,
                                            H5P_DEFAULT, H5P_DEFAULT);

            /* data */
            {
                hid_t dset = write_1d_dataset(this_graph, "data",
                                               H5T_NATIVE_DOUBLE,
                                               (hsize_t)XLENGTH(g_x),
                                               REAL(g_x), gzip);
                if (dset >= 0) H5Dclose(dset);
            }
            /* indices */
            {
                hid_t dset = write_1d_dataset(this_graph, "indices",
                                               H5T_NATIVE_INT,
                                               (hsize_t)XLENGTH(g_i),
                                               INTEGER(g_i), gzip);
                if (dset >= 0) H5Dclose(dset);
            }
            /* indptr */
            {
                hid_t dset = write_1d_dataset(this_graph, "indptr",
                                               H5T_NATIVE_INT,
                                               (hsize_t)XLENGTH(g_p),
                                               INTEGER(g_p), gzip);
                if (dset >= 0) H5Dclose(dset);
            }

            /* CSR encoding + shape (square: n_cells x n_cells) */
            set_str_attr_ascii(this_graph, "encoding-type", "csr_matrix");
            set_str_attr_ascii(this_graph, "encoding-version", "0.1.0");
            {
                int shape[2] = { INTEGER(g_dim)[0], INTEGER(g_dim)[1] };
                hsize_t two = 2;
                hid_t space = H5Screate_simple(1, &two, NULL);
                hid_t attr = H5Acreate2(this_graph, "shape", H5T_NATIVE_INT,
                                         space, H5P_DEFAULT, H5P_DEFAULT);
                if (attr >= 0) {
                    H5Awrite(attr, H5T_NATIVE_INT, shape);
                    H5Aclose(attr);
                }
                H5Sclose(space);
            }

            H5Gclose(this_graph);
        }

        H5Gclose(obsp_grp);
    }

    /* ── Empty groups: varm, varp, uns, layers, raw ──────────────────────────── */
    {
        const char *empty_groups[] = { "varm", "varp", "uns", "layers" };
        for (int ei = 0; ei < 4; ei++) {
            hid_t grp = H5Gcreate2(file, empty_groups[ei], H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
            if (grp >= 0) {
                set_str_attr_ascii(grp, "encoding-type", "dict");
                set_str_attr_ascii(grp, "encoding-version", "0.1.0");
                H5Gclose(grp);
            }
        }
    }

    H5Fclose(file);
    return ScalarLogical(TRUE);
}

/* Registration is in init.c */
