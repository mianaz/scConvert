/*
 * sc_sparse.c — Streaming sparse matrix read/write for HDF5
 *
 * Handles CSR/CSC matrices in AnnData format: data, indices, indptr arrays
 * stored as separate datasets within an HDF5 group.
 *
 * Key optimization: CSC of A == CSR of A^T, so transposing a dgCMatrix
 * (genes×cells CSC) to AnnData's cells×genes CSR is a zero-copy pointer
 * reinterpretation — no data movement needed, just relabel the arrays.
 */

#include "sc_convert.h"

/* ── Read CSR info from an HDF5 group ───────────────────────────────────────── */

/*
 * Read shape/dims from a 2-element attribute into shape[2].
 * Returns 1 if found, 0 if not.
 */
static int read_shape_attr(hid_t grp, const char *attr_name, int64_t shape[2]) {
    if (H5Aexists(grp, attr_name) <= 0)
        return 0;

    hid_t attr = H5Aopen(grp, attr_name, H5P_DEFAULT);
    hid_t space = H5Aget_space(attr);
    hsize_t dims;
    H5Sget_simple_extent_dims(space, &dims, NULL);

    int ok = 0;
    if (dims == 2) {
        hid_t atype = H5Aget_type(attr);
        if (H5Tget_size(atype) <= 4) {
            int32_t s32[2];
            H5Aread(attr, H5T_NATIVE_INT32, s32);
            shape[0] = s32[0];
            shape[1] = s32[1];
        } else {
            H5Aread(attr, H5T_NATIVE_INT64, shape);
        }
        H5Tclose(atype);
        ok = 1;
    }
    H5Sclose(space);
    H5Aclose(attr);
    return ok;
}

int sc_read_csr_info(hid_t grp, sc_csr_info_t *info) {
    int64_t shape[2] = {0, 0};

    /* Try h5ad "shape" first (n_cells, n_genes for CSR) */
    if (read_shape_attr(grp, "shape", shape)) {
        info->n_rows = (hsize_t)shape[0];
        info->n_cols = (hsize_t)shape[1];
    }
    /* Try h5seurat "dims" (n_genes, n_cells) — note: reversed ordering! */
    else if (read_shape_attr(grp, "dims", shape)) {
        /* h5seurat dims = (n_genes, n_cells), stored as CSC
         * We normalize to (n_pointer_dim, n_index_dim) from indptr perspective */
        info->n_rows = (hsize_t)shape[0];
        info->n_cols = (hsize_t)shape[1];
    }

    /* Get nnz from data dataset */
    if (sc_has_dataset(grp, "data")) {
        hid_t dset = H5Dopen2(grp, "data", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t dims;
        H5Sget_simple_extent_dims(space, &dims, NULL);
        info->nnz = dims;
        H5Sclose(space);
        H5Dclose(dset);
    }

    /* Derive dimensions from indptr if shape/dims was missing */
    if (sc_has_dataset(grp, "indptr")) {
        hid_t dset = H5Dopen2(grp, "indptr", H5P_DEFAULT);
        hid_t space = H5Dget_space(dset);
        hsize_t dims;
        H5Sget_simple_extent_dims(space, &dims, NULL);
        /* indptr has (pointer_dim + 1) entries */
        hsize_t ptr_dim = dims - 1;
        if (info->n_rows == 0 && info->n_cols == 0) {
            /* No shape attr at all — indptr tells us the pointer dimension */
            info->n_rows = ptr_dim;
        }
        H5Sclose(space);
        H5Dclose(dset);
    }

    return SC_OK;
}

/* ── Check if source dataset uses gzip (any level) ───────────────────────────── */

static int src_uses_gzip(hid_t src_dset, hsize_t *src_chunk_sz) {
    hid_t dcpl = H5Dget_create_plist(src_dset);
    if (dcpl < 0) return 0;

    /* Check chunked layout */
    H5D_layout_t layout = H5Pget_layout(dcpl);
    if (layout != H5D_CHUNKED) {
        H5Pclose(dcpl);
        return 0;
    }

    /* Get source chunk size */
    hsize_t chunk;
    if (H5Pget_chunk(dcpl, 1, &chunk) < 0) {
        H5Pclose(dcpl);
        return 0;
    }
    *src_chunk_sz = chunk;

    /* Check filters — need at least one gzip filter (any level) */
    int nfilters = H5Pget_nfilters(dcpl);
    int has_gzip = 0;
    for (int i = 0; i < nfilters; i++) {
        unsigned int flags;
        size_t cd_nelmts = 1;
        unsigned int cd_values[1];
        H5Z_filter_t filt = H5Pget_filter2(dcpl, (unsigned)i, &flags,
                                             &cd_nelmts, cd_values, 0, NULL, NULL);
        if (filt == H5Z_FILTER_DEFLATE) {
            has_gzip = 1;
        }
    }

    H5Pclose(dcpl);
    return has_gzip;
}

/* ── Direct chunk copy: skip decompress/recompress when settings match ──────── */

static int stream_1d_dataset_direct(hid_t src_dset, hid_t dst_dset,
                                      hsize_t total, hsize_t chunk_sz) {
    hsize_t n_chunks = (total + chunk_sz - 1) / chunk_sz;

    for (hsize_t ci = 0; ci < n_chunks; ci++) {
        hsize_t offset = ci * chunk_sz;
        uint32_t filter_mask = 0;
        hsize_t chunk_nbytes = 0;

        /* Get chunk metadata (compressed size) */
        if (H5Dget_chunk_info(src_dset, H5S_ALL, ci, NULL, &filter_mask,
                               NULL, &chunk_nbytes) < 0)
            return SC_ERR_HDF;

        void *buf = malloc(chunk_nbytes);
        if (!buf) return SC_ERR;

        /* Read compressed chunk directly (no decompression).
         * HDF5 2.0+ renamed to H5Dread_chunk2 and added a buf_size in/out
         * param; the 1.x API takes exactly 5 arguments without it. */
        size_t buf_size = (size_t)chunk_nbytes;
#ifndef H5_VERSION_GE
#define H5_VERSION_GE(a,b,c) 0
#endif
#if H5_VERSION_GE(2,0,0)
        if (H5Dread_chunk(src_dset, H5P_DEFAULT, &offset,
                           &filter_mask, buf, &buf_size) < 0) {
#else
        if (H5Dread_chunk(src_dset, H5P_DEFAULT, &offset,
                           &filter_mask, buf) < 0) {
#endif
            free(buf);
            return SC_ERR_HDF;
        }

        /* Write compressed chunk directly (no recompression) */
        if (H5Dwrite_chunk(dst_dset, H5P_DEFAULT, filter_mask,
                            &offset, buf_size, buf) < 0) {
            free(buf);
            return SC_ERR_HDF;
        }

        free(buf);
    }
    return SC_OK;
}

/* ── Stream a single 1D dataset chunk-by-chunk ──────────────────────────────── */

static int stream_1d_dataset(hid_t src_dset, hid_t dst_grp,
                              const char *name, hid_t mem_type,
                              hsize_t total, int gzip_level)
{
    hsize_t chunk_sz = SC_CHUNK_SIZE;
    if (chunk_sz > total) chunk_sz = total;

    /* Check if we can use direct chunk copy.
     * Accept any source gzip level — we copy compressed bytes as-is,
     * which avoids the expensive decompress/recompress cycle.
     * Match destination chunk size to source so chunk boundaries align. */
    hsize_t src_chunk_sz = 0;
    int can_direct = (gzip_level > 0 &&
                      src_uses_gzip(src_dset, &src_chunk_sz));

    /* Create destination dataset with chunking + gzip.
     * For direct copy: use source chunk size so H5Dwrite_chunk works.
     * For standard path: use SC_CHUNK_SIZE.
     * Cap at total — HDF5 disallows chunk_size > dataset_size.
     * If capping changes chunk size, disable direct copy (chunk boundaries won't align). */
    hsize_t dst_chunk = can_direct ? src_chunk_sz : chunk_sz;
    if (dst_chunk > total) {
        dst_chunk = total;
        can_direct = 0;  /* chunk sizes mismatched, fall back to standard path */
    }
    hid_t dst_space = H5Screate_simple(1, &total, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 1, &dst_chunk);
    if (gzip_level > 0)
        H5Pset_deflate(dcpl, (unsigned)gzip_level);

    hid_t src_type = H5Dget_type(src_dset);
    hid_t dst_dset = H5Dcreate2(dst_grp, name, src_type, dst_space,
                                 H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Tclose(src_type);
    H5Pclose(dcpl);

    if (dst_dset < 0) {
        H5Sclose(dst_space);
        return SC_ERR_HDF;
    }

    /* Fast path: direct chunk copy skips decompress/recompress */
    if (can_direct) {
        int rc = stream_1d_dataset_direct(src_dset, dst_dset, total, src_chunk_sz);
        if (rc == SC_OK) {
            H5Dclose(dst_dset);
            H5Sclose(dst_space);
            return SC_OK;
        }
        /* Fall through to standard path on failure */
    }

    /* Standard path: decompress, read into memory, recompress on write */
    size_t elem_size = H5Tget_size(mem_type);
    void *buf = malloc(chunk_sz * elem_size);
    if (!buf) {
        H5Dclose(dst_dset);
        H5Sclose(dst_space);
        return SC_ERR;
    }

    hid_t src_space = H5Dget_space(src_dset);
    hid_t mem_space;
    hsize_t offset = 0;

    while (offset < total) {
        hsize_t count = chunk_sz;
        if (offset + count > total)
            count = total - offset;

        /* Select hyperslab in source */
        H5Sselect_hyperslab(src_space, H5S_SELECT_SET, &offset, NULL, &count, NULL);

        /* Memory space for this chunk */
        mem_space = H5Screate_simple(1, &count, NULL);

        /* Read from source */
        if (H5Dread(src_dset, mem_type, mem_space, src_space, H5P_DEFAULT, buf) < 0) {
            H5Sclose(mem_space);
            goto fail;
        }

        /* Select hyperslab in destination */
        hid_t dst_fspace = H5Dget_space(dst_dset);
        H5Sselect_hyperslab(dst_fspace, H5S_SELECT_SET, &offset, NULL, &count, NULL);

        /* Write to destination */
        if (H5Dwrite(dst_dset, mem_type, mem_space, dst_fspace, H5P_DEFAULT, buf) < 0) {
            H5Sclose(dst_fspace);
            H5Sclose(mem_space);
            goto fail;
        }

        H5Sclose(dst_fspace);
        H5Sclose(mem_space);
        offset += count;
    }

    free(buf);
    H5Sclose(src_space);
    H5Dclose(dst_dset);
    H5Sclose(dst_space);
    return SC_OK;

fail:
    free(buf);
    H5Sclose(src_space);
    H5Dclose(dst_dset);
    H5Sclose(dst_space);
    return SC_ERR_HDF;
}

/* ── Copy CSR group without transposing (direct copy) ───────────────────────── */

int sc_stream_csr_copy(hid_t src_grp, hid_t dst_grp, int gzip_level) {
    sc_csr_info_t info = {0};
    int rc = sc_read_csr_info(src_grp, &info);
    if (rc != SC_OK) return rc;

    /* NOTE: This function only copies data/indices/indptr datasets.
     * Callers must set the appropriate encoding attributes (shape vs dims,
     * encoding-type, etc.) for the target format. */

    /* Stream data (float64) */
    hid_t src_data = H5Dopen2(src_grp, "data", H5P_DEFAULT);
    rc = stream_1d_dataset(src_data, dst_grp, "data", H5T_NATIVE_DOUBLE,
                           info.nnz, gzip_level);
    H5Dclose(src_data);
    if (rc != SC_OK) return rc;

    /* Stream indices — detect native type to handle int32 or int64 */
    hid_t src_idx = H5Dopen2(src_grp, "indices", H5P_DEFAULT);
    {
        hid_t idx_type = H5Dget_type(src_idx);
        hid_t idx_mem = H5Tget_native_type(idx_type, H5T_DIR_ASCEND);
        rc = stream_1d_dataset(src_idx, dst_grp, "indices", idx_mem,
                               info.nnz, gzip_level);
        H5Tclose(idx_mem); H5Tclose(idx_type);
    }
    H5Dclose(src_idx);
    if (rc != SC_OK) return rc;

    /* Stream indptr — detect native type to handle int32 or int64 */
    hid_t src_ptr = H5Dopen2(src_grp, "indptr", H5P_DEFAULT);
    hid_t ptr_space = H5Dget_space(src_ptr);
    hsize_t ptr_len;
    H5Sget_simple_extent_dims(ptr_space, &ptr_len, NULL);
    H5Sclose(ptr_space);

    {
        hid_t ptr_type = H5Dget_type(src_ptr);
        hid_t ptr_mem = H5Tget_native_type(ptr_type, H5T_DIR_ASCEND);
        rc = stream_1d_dataset(src_ptr, dst_grp, "indptr", ptr_mem,
                               ptr_len, gzip_level);
        H5Tclose(ptr_mem); H5Tclose(ptr_type);
    }
    H5Dclose(src_ptr);
    return rc;
}

/* ── Copy CSR group with transpose (CSR ↔ CSC reinterpretation) ─────────────
 *
 * AnnData stores cells×genes CSR.  h5Seurat stores genes×cells CSC.
 * CSR of A == CSC of A^T — same arrays, different interpretation.
 *
 * When going h5ad→h5seurat:
 *   - Read CSR (cells×genes): data, indices (col idx), indptr (row ptr)
 *   - Write as CSC (genes×cells): data, indices (row idx = same), indptr (col ptr = same)
 *   - Just need to swap n_rows/n_cols in the dims attribute
 *
 * When going h5seurat→h5ad:
 *   - Read CSC (genes×cells) from h5Seurat
 *   - Write as CSR (cells×genes) — same arrays, swap dims
 *
 * In both cases, the raw data/indices/indptr arrays are identical.
 * We just stream-copy them and fix the shape/dims metadata.
 */
int sc_stream_csr_transpose(hid_t src_grp, hid_t dst_grp, int gzip_level) {
    sc_csr_info_t info = {0};
    int rc = sc_read_csr_info(src_grp, &info);
    if (rc != SC_OK) return rc;

    /* NOTE: This function only copies data/indices/indptr datasets.
     * The caller must set shape/dims and encoding attributes for the target. */
    (void)info; /* n_rows/n_cols used by caller via sc_read_csr_info */

    /* Stream data (values are identical in CSR/CSC) */
    hid_t src_data = H5Dopen2(src_grp, "data", H5P_DEFAULT);
    rc = stream_1d_dataset(src_data, dst_grp, "data", H5T_NATIVE_DOUBLE,
                           info.nnz, gzip_level);
    H5Dclose(src_data);
    if (rc != SC_OK) return rc;

    /* Stream indices — detect native type to handle int32 or int64 */
    hid_t src_idx = H5Dopen2(src_grp, "indices", H5P_DEFAULT);
    {
        hid_t idx_type = H5Dget_type(src_idx);
        hid_t idx_mem = H5Tget_native_type(idx_type, H5T_DIR_ASCEND);
        rc = stream_1d_dataset(src_idx, dst_grp, "indices", idx_mem,
                               info.nnz, gzip_level);
        H5Tclose(idx_mem); H5Tclose(idx_type);
    }
    H5Dclose(src_idx);
    if (rc != SC_OK) return rc;

    /* Stream indptr — detect native type to handle int32 or int64 */
    hid_t src_ptr = H5Dopen2(src_grp, "indptr", H5P_DEFAULT);
    hid_t ptr_space = H5Dget_space(src_ptr);
    hsize_t ptr_len;
    H5Sget_simple_extent_dims(ptr_space, &ptr_len, NULL);
    H5Sclose(ptr_space);

    {
        hid_t ptr_type = H5Dget_type(src_ptr);
        hid_t ptr_mem = H5Tget_native_type(ptr_type, H5T_DIR_ASCEND);
        rc = stream_1d_dataset(src_ptr, dst_grp, "indptr", ptr_mem,
                               ptr_len, gzip_level);
        H5Tclose(ptr_mem); H5Tclose(ptr_type);
    }
    H5Dclose(src_ptr);
    return rc;
}

/* ── Copy a 2D dataset with transpose ([M,N] → [N,M]) ──────────────────────── */

int sc_copy_2d_transposed(hid_t src_dset, hid_t dst_grp,
                            const char *name, int gzip_level)
{
    hid_t src_space = H5Dget_space(src_dset);
    int ndims = H5Sget_simple_extent_ndims(src_space);
    if (ndims != 2) {
        H5Sclose(src_space);
        return SC_ERR_ARG;
    }

    hsize_t src_dims[2];
    H5Sget_simple_extent_dims(src_space, src_dims, NULL);

    hid_t src_type = H5Dget_type(src_dset);
    hid_t mem_type = H5Tget_native_type(src_type, H5T_DIR_ASCEND);
    size_t elem_size = H5Tget_size(mem_type);

    hsize_t M = src_dims[0], N = src_dims[1];
    hsize_t total = M * N;

    /* Overflow guard: total * elem_size must fit in size_t */
    if (elem_size > 0 && total > (hsize_t)(SIZE_MAX / elem_size)) {
        fprintf(stderr, "Error: matrix too large for transpose (%llu x %llu)\n",
                (unsigned long long)M, (unsigned long long)N);
        H5Tclose(mem_type); H5Tclose(src_type); H5Sclose(src_space);
        return SC_ERR;
    }

    /* Read source */
    void *src_buf = malloc(total * elem_size);
    if (!src_buf) {
        H5Tclose(mem_type); H5Tclose(src_type); H5Sclose(src_space);
        return SC_ERR;
    }
    H5Dread(src_dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_buf);

    /* Transpose: src[i*N + j] → dst[j*M + i] */
    void *dst_buf = malloc(total * elem_size);
    if (!dst_buf) {
        free(src_buf);
        H5Tclose(mem_type); H5Tclose(src_type); H5Sclose(src_space);
        return SC_ERR;
    }

    unsigned char *s = (unsigned char *)src_buf;
    unsigned char *d = (unsigned char *)dst_buf;
    for (hsize_t i = 0; i < M; i++) {
        for (hsize_t j = 0; j < N; j++) {
            memcpy(d + (j * M + i) * elem_size,
                   s + (i * N + j) * elem_size,
                   elem_size);
        }
    }

    /* Create transposed dataset [N, M] */
    hsize_t dst_dims[2] = {N, M};
    hid_t dst_space = H5Screate_simple(2, dst_dims, NULL);
    hsize_t chunks[2] = {N, M};
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunks);
    if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);

    hid_t dst_dset = H5Dcreate2(dst_grp, name, src_type, dst_space,
                                 H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dwrite(dst_dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst_buf);

    free(src_buf);
    free(dst_buf);
    H5Dclose(dst_dset);
    H5Pclose(dcpl);
    H5Sclose(dst_space);
    H5Tclose(mem_type);
    H5Tclose(src_type);
    H5Sclose(src_space);
    return SC_OK;
}

/* ── Copy a generic 1D dataset (used by dataframe and group copiers) ────────── */

int sc_copy_dataset_chunked(hid_t src_dset, hid_t dst_grp,
                             const char *name, int gzip_level)
{
    hid_t src_space = H5Dget_space(src_dset);
    int ndims = H5Sget_simple_extent_ndims(src_space);

    if (ndims == 0) {
        /* Scalar dataset — direct copy */
        hid_t src_type = H5Dget_type(src_dset);
        size_t sz = H5Tget_size(src_type);
        hid_t dst_space_s = H5Screate(H5S_SCALAR);
        hid_t dst_dset = H5Dcreate2(dst_grp, name, src_type, dst_space_s,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        void *buf = malloc(sz);
        if (!buf) {
            H5Dclose(dst_dset);
            H5Sclose(dst_space_s);
            H5Tclose(src_type);
            H5Sclose(src_space);
            return SC_ERR;
        }
        H5Dread(src_dset, src_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dwrite(dst_dset, src_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        free(buf);
        H5Dclose(dst_dset);
        H5Sclose(dst_space_s);
        H5Tclose(src_type);
        H5Sclose(src_space);
        return SC_OK;
    }

    if (ndims == 1) {
        hsize_t dims;
        H5Sget_simple_extent_dims(src_space, &dims, NULL);
        H5Sclose(src_space);

        hid_t src_type = H5Dget_type(src_dset);

        /* Check if it's a string type */
        if (H5Tget_class(src_type) == H5T_STRING) {
            /* String dataset — read all at once (strings are typically small) */
            hid_t memtype = sc_create_vlen_str_type();
            char **strs = (char **)calloc(dims, sizeof(char *));
            if (!strs) {
                H5Tclose(memtype);
                H5Tclose(src_type);
                return SC_ERR;
            }
            H5Dread(src_dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);

            /* Write */
            hsize_t chunk = dims;
            hid_t dst_space_d = H5Screate_simple(1, &dims, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(dcpl, 1, &chunk);
            if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);

            hid_t dst_dset = H5Dcreate2(dst_grp, name, memtype, dst_space_d,
                                         H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Dwrite(dst_dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);

            /* Cleanup */
            hid_t reclaim_space = H5Screate_simple(1, &dims, NULL);
            H5Treclaim(memtype, reclaim_space, H5P_DEFAULT, strs);
            H5Sclose(reclaim_space);
            free(strs);

            H5Dclose(dst_dset);
            H5Pclose(dcpl);
            H5Sclose(dst_space_d);
            H5Tclose(memtype);
            H5Tclose(src_type);
            return SC_OK;
        }

        /* Numeric dataset — stream */
        hid_t mem_type = H5Tget_native_type(src_type, H5T_DIR_ASCEND);
        int rc = stream_1d_dataset(src_dset, dst_grp, name, mem_type,
                                   dims, gzip_level);
        H5Tclose(mem_type);
        H5Tclose(src_type);
        return rc;
    }

    /* 2D dataset — read/write whole (embeddings are typically small) */
    hsize_t dims2[2];
    H5Sget_simple_extent_dims(src_space, dims2, NULL);

    hid_t src_type = H5Dget_type(src_dset);
    hid_t mem_type = H5Tget_native_type(src_type, H5T_DIR_ASCEND);
    size_t elem_size = H5Tget_size(mem_type);

    hsize_t total = dims2[0] * dims2[1];

    /* Overflow guard: total * elem_size must fit in size_t */
    if (elem_size > 0 && total > (hsize_t)(SIZE_MAX / elem_size)) {
        fprintf(stderr, "Error: 2D dataset too large for copy (%llu x %llu)\n",
                (unsigned long long)dims2[0], (unsigned long long)dims2[1]);
        H5Tclose(mem_type); H5Tclose(src_type); H5Sclose(src_space);
        return SC_ERR;
    }

    void *buf = malloc(total * elem_size);
    if (!buf) {
        H5Tclose(mem_type);
        H5Tclose(src_type);
        H5Sclose(src_space);
        return SC_ERR;
    }

    H5Dread(src_dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    hsize_t chunks2[2] = {dims2[0], dims2[1]};
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunks2);
    if (gzip_level > 0) H5Pset_deflate(dcpl, (unsigned)gzip_level);

    hid_t dst_dset = H5Dcreate2(dst_grp, name, src_type, src_space,
                                 H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dwrite(dst_dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    free(buf);
    H5Dclose(dst_dset);
    H5Pclose(dcpl);
    H5Tclose(mem_type);
    H5Tclose(src_type);
    H5Sclose(src_space);
    return SC_OK;
}
