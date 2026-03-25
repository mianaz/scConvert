/*
 * sc_zarr.c — Zarr v2 ↔ h5seurat conversion for scConvert CLI
 *
 * Reads/writes zarr AnnData stores (directory-based format).
 * Composite conversions via temp h5seurat for zarr ↔ {h5ad, h5mu, loom}.
 */

#include "sc_convert.h"
#include <zlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

/* ══════════════════════════════════════════════════════════════════════════════
 *  Directory utilities
 * ══════════════════════════════════════════════════════════════════════════════ */

int sc_mkdir_p(const char *path) {
    char buf[2048];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char *p = buf + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            mkdir(buf, 0755);
            *p = '/';
        }
    }
    return mkdir(buf, 0755) == 0 || errno == EEXIST ? SC_OK : SC_ERR_IO;
}

int sc_rmdir_recursive(const char *path) {
    DIR *d = opendir(path);
    if (!d) { unlink(path); return SC_OK; }
    struct dirent *ent;
    while ((ent = readdir(d))) {
        if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0)
            continue;
        char child[2048];
        snprintf(child, sizeof(child), "%s/%s", path, ent->d_name);
        struct stat st;
        if (stat(child, &st) == 0 && S_ISDIR(st.st_mode))
            sc_rmdir_recursive(child);
        else
            unlink(child);
    }
    closedir(d);
    rmdir(path);
    return SC_OK;
}

/* ── List child directories (excluding . files) ─────────────────────────────── */

static int zarr_list_children(const char *dir, char ***out_names, int *out_n) {
    DIR *d = opendir(dir);
    if (!d) { *out_names = NULL; *out_n = 0; return SC_ERR_IO; }

    int cap = 64, n = 0;
    char **names = (char **)calloc((size_t)cap, sizeof(char *));
    if (!names) { closedir(d); *out_names = NULL; *out_n = 0; return SC_ERR; }

    struct dirent *ent;
    while ((ent = readdir(d))) {
        if (ent->d_name[0] == '.') continue;
        /* Must be a directory or a file — check for .zarray or .zgroup
         * to distinguish zarr children from chunk files */
        char test[2048];
        struct stat st;
        snprintf(test, sizeof(test), "%s/%s", dir, ent->d_name);
        if (stat(test, &st) != 0) continue;
        if (!S_ISDIR(st.st_mode)) continue;
        if (n >= cap) {
            cap *= 2;
            char **tmp = (char **)realloc(names, (size_t)cap * sizeof(char *));
            if (!tmp) {
                for (int i = 0; i < n; i++) free(names[i]);
                free(names);
                closedir(d);
                *out_names = NULL;
                *out_n = 0;
                return SC_ERR;
            }
            names = tmp;
        }
        names[n++] = strdup(ent->d_name);
    }
    closedir(d);
    *out_names = names;
    *out_n = n;
    return SC_OK;
}

static void zarr_free_children(char **names, int n) {
    if (!names) return;
    for (int i = 0; i < n; i++) free(names[i]);
    free(names);
}

/* ── Check zarr node types ──────────────────────────────────────────────────── */

static int zarr_exists(const char *store, const char *rel) {
    char path[2048];
    snprintf(path, sizeof(path), "%s/%s", store, rel);
    struct stat st;
    return stat(path, &st) == 0;
}

static int zarr_is_group(const char *store, const char *rel) {
    char path[2048];
    snprintf(path, sizeof(path), "%s/%s/.zgroup", store, rel);
    return access(path, F_OK) == 0;
}

static int zarr_is_array(const char *store, const char *rel) {
    char path[2048];
    snprintf(path, sizeof(path), "%s/%s/.zarray", store, rel);
    return access(path, F_OK) == 0;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Zarr chunk I/O with zlib
 * ══════════════════════════════════════════════════════════════════════════════ */

/* Read and decompress a zarr chunk file. Returns malloc'd buffer. */
static void *zarr_read_chunk(const char *chunk_path, const sc_zarr_meta_t *meta,
                              size_t *out_size) {
    FILE *f = fopen(chunk_path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long clen = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (clen <= 0) { fclose(f); return NULL; }

    void *cbuf = malloc((size_t)clen);
    if (!cbuf) { fclose(f); return NULL; }
    fread(cbuf, 1, (size_t)clen, f);
    fclose(f);

    int use_zlib = (strcmp(meta->compressor_id, "zlib") == 0 ||
                    strcmp(meta->compressor_id, "gzip") == 0);
    if (!use_zlib) {
        /* uncompressed */
        *out_size = (size_t)clen;
        return cbuf;
    }

    /* Decompress using zlib */
    /* Estimate: start with 4x, grow if needed */
    size_t dbuf_sz = (size_t)clen * 4;
    if (dbuf_sz < 4096) dbuf_sz = 4096;
    void *dbuf = malloc(dbuf_sz);
    if (!dbuf) { free(cbuf); return NULL; }

    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    /* Use raw deflate for zlib codec, auto-detect for gzip */
    int wbits = (strcmp(meta->compressor_id, "gzip") == 0) ? 15 + 16 : 15;
    if (inflateInit2(&zs, wbits) != Z_OK) {
        free(cbuf); free(dbuf); return NULL;
    }
    zs.next_in = (Bytef *)cbuf;
    zs.avail_in = (uInt)clen;
    zs.next_out = (Bytef *)dbuf;
    zs.avail_out = (uInt)dbuf_sz;

    int ret;
    while ((ret = inflate(&zs, Z_NO_FLUSH)) != Z_STREAM_END) {
        if (ret == Z_BUF_ERROR || ret == Z_OK) {
            size_t done = dbuf_sz - zs.avail_out;
            dbuf_sz *= 2;
            void *tmp = realloc(dbuf, dbuf_sz);
            if (!tmp) {
                inflateEnd(&zs);
                free(cbuf); free(dbuf);
                return NULL;
            }
            dbuf = tmp;
            zs.next_out = (Bytef *)dbuf + done;
            zs.avail_out = (uInt)(dbuf_sz - done);
        } else {
            inflateEnd(&zs);
            free(cbuf); free(dbuf);
            return NULL;
        }
    }
    *out_size = dbuf_sz - zs.avail_out;
    inflateEnd(&zs);
    free(cbuf);
    return dbuf;
}

/* Compress and write a zarr chunk file */
static int zarr_write_chunk(const char *chunk_path, const void *data,
                             size_t data_size, int gzip_level) {
    if (gzip_level <= 0) {
        /* Uncompressed */
        FILE *f = fopen(chunk_path, "wb");
        if (!f) return SC_ERR_IO;
        fwrite(data, 1, data_size, f);
        fclose(f);
        return SC_OK;
    }

    /* Compress with zlib (deflate, not gzip wrapper) */
    uLongf clen = compressBound((uLong)data_size);
    void *cbuf = malloc(clen);
    if (!cbuf) return SC_ERR;
    if (compress2((Bytef *)cbuf, &clen, (const Bytef *)data,
                   (uLong)data_size, gzip_level) != Z_OK) {
        free(cbuf);
        return SC_ERR;
    }
    FILE *f = fopen(chunk_path, "wb");
    if (!f) { free(cbuf); return SC_ERR_IO; }
    fwrite(cbuf, 1, clen, f);
    fclose(f);
    free(cbuf);
    return SC_OK;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Typed array I/O (single-chunk)
 * ══════════════════════════════════════════════════════════════════════════════ */

/* Read .zarray metadata from a zarr array directory */
static int zarr_read_meta(const char *arr_dir, sc_zarr_meta_t *meta) {
    char path[2048];
    snprintf(path, sizeof(path), "%s/.zarray", arr_dir);
    char *json = sc_json_read_file(path);
    if (!json) return SC_ERR_IO;
    int rc = sc_json_parse_zarray(json, meta);
    free(json);
    return rc;
}

/* Element size from dtype */
static int dtype_size(const char *dtype) {
    if (strcmp(dtype, "<f8") == 0 || strcmp(dtype, "<f4") == 0 ||
        strcmp(dtype, "<i8") == 0)
        return (dtype[2] - '0');
    if (strcmp(dtype, "<i4") == 0 || strcmp(dtype, "<i2") == 0)
        return (dtype[2] - '0');
    if (strcmp(dtype, "|i1") == 0 || strcmp(dtype, "|b1") == 0 ||
        strcmp(dtype, "|u1") == 0)
        return 1;
    if (strcmp(dtype, "|O") == 0) return 0; /* vlen-utf8, variable */
    return 8; /* default */
}

/* Read 1D numeric array. Returns malloc'd buffer. Sets *out_n. */
static void *zarr_read_numeric(const char *arr_dir, int64_t *out_n) {
    sc_zarr_meta_t meta;
    if (zarr_read_meta(arr_dir, &meta) != 0) return NULL;

    int64_t total = meta.shape[0];
    for (int i = 1; i < meta.ndim; i++) total *= meta.shape[i];
    *out_n = total;
    if (total == 0) return calloc(1, 1);

    int elem_sz = dtype_size(meta.dtype);
    if (elem_sz == 0) return NULL; /* vlen-utf8 -- use zarr_read_strings */

    /* Single-chunk fast path */
    char chunk_path[2048];
    if (meta.ndim == 1)
        snprintf(chunk_path, sizeof(chunk_path), "%s/0", arr_dir);
    else
        snprintf(chunk_path, sizeof(chunk_path), "%s/0.0", arr_dir);

    /* Check for multi-chunk */
    int n_chunks = 1;
    for (int i = 0; i < meta.ndim; i++) {
        int64_t nc = (meta.shape[i] + meta.chunks[i] - 1) / meta.chunks[i];
        n_chunks *= (int)nc;
    }

    if (n_chunks == 1) {
        size_t out_sz;
        void *data = zarr_read_chunk(chunk_path, &meta, &out_sz);
        return data;
    }

    /* Multi-chunk: allocate full buffer, read chunk by chunk */
    size_t full_sz = (size_t)total * (size_t)elem_sz;
    void *full = calloc(1, full_sz);
    if (!full) return NULL;

    if (meta.ndim == 1) {
        int64_t chunk_elems = meta.chunks[0];
        int64_t nc = (meta.shape[0] + chunk_elems - 1) / chunk_elems;
        for (int64_t ci = 0; ci < nc; ci++) {
            snprintf(chunk_path, sizeof(chunk_path), "%s/%lld", arr_dir, (long long)ci);
            size_t out_sz;
            void *cdata = zarr_read_chunk(chunk_path, &meta, &out_sz);
            if (!cdata) continue;
            int64_t offset = ci * chunk_elems;
            int64_t elems = (offset + chunk_elems > meta.shape[0]) ?
                            meta.shape[0] - offset : chunk_elems;
            memcpy((char *)full + offset * elem_sz, cdata, (size_t)elems * elem_sz);
            free(cdata);
        }
    } else if (meta.ndim == 2) {
        /* 2D multi-chunk: iterate over (ci, cj) in row-major order */
        int64_t nchunks_0 = (meta.shape[0] + meta.chunks[0] - 1) / meta.chunks[0];
        int64_t nchunks_1 = (meta.shape[1] + meta.chunks[1] - 1) / meta.chunks[1];
        for (int64_t ci = 0; ci < nchunks_0; ci++) {
            for (int64_t cj = 0; cj < nchunks_1; cj++) {
                snprintf(chunk_path, sizeof(chunk_path), "%s/%lld.%lld",
                         arr_dir, (long long)ci, (long long)cj);
                size_t out_sz;
                void *cdata = zarr_read_chunk(chunk_path, &meta, &out_sz);
                if (!cdata) continue;

                /* Row range covered by this chunk */
                int64_t row_start = ci * meta.chunks[0];
                int64_t row_end   = row_start + meta.chunks[0];
                if (row_end > meta.shape[0]) row_end = meta.shape[0];

                /* Col range covered by this chunk */
                int64_t col_start = cj * meta.chunks[1];
                int64_t col_end   = col_start + meta.chunks[1];
                if (col_end > meta.shape[1]) col_end = meta.shape[1];

                int64_t chunk_cols = meta.chunks[1]; /* full chunk width (stride) */
                int64_t out_cols   = meta.shape[1];  /* full output width */

                for (int64_t r = row_start; r < row_end; r++) {
                    int64_t src_row   = r - row_start;
                    int64_t col_count = col_end - col_start;
                    memcpy((char *)full + (r * out_cols + col_start) * elem_sz,
                           (char *)cdata + src_row * chunk_cols * elem_sz,
                           (size_t)col_count * (size_t)elem_sz);
                }
                free(cdata);
            }
        }
    } else {
        /* ndim > 2: not expected in AnnData; return what we have (zeros) */
        SC_MSG("Warning: zarr_read_numeric: ndim=%d not fully supported\n",
                meta.ndim);
    }

    return full;
}

/* ── Read vlen-utf8 string array ────────────────────────────────────────────── */

static char **zarr_read_strings(const char *arr_dir, int64_t *out_n) {
    sc_zarr_meta_t meta;
    if (zarr_read_meta(arr_dir, &meta) != 0) return NULL;

    *out_n = meta.shape[0];
    if (meta.shape[0] == 0) return NULL;

    int64_t n           = meta.shape[0];
    int64_t chunk_size  = meta.chunks[0];
    int64_t n_chunks    = (n + chunk_size - 1) / chunk_size;

    char **strs = (char **)calloc((size_t)n, sizeof(char *));
    if (!strs) return NULL;

    char chunk_path[2048];
    int64_t global_idx = 0; /* next string slot to fill */

    for (int64_t ci = 0; ci < n_chunks; ci++) {
        snprintf(chunk_path, sizeof(chunk_path), "%s/%lld", arr_dir, (long long)ci);

        size_t raw_sz;
        void *raw = zarr_read_chunk(chunk_path, &meta, &raw_sz);
        if (!raw) continue; /* chunk missing — leave NULLs for these entries */

        /* Number of strings in this chunk */
        int64_t chunk_n = chunk_size;
        if (global_idx + chunk_n > n) chunk_n = n - global_idx;

        const uint8_t *p   = (const uint8_t *)raw;
        const uint8_t *end = p + raw_sz;

        /* Check for numcodecs header: first uint32 == chunk_n means header present */
        if (raw_sz >= 4) {
            uint32_t first;
            memcpy(&first, p, 4);
            if ((int64_t)first == chunk_n || (int64_t)first == chunk_size) {
                p += 4; /* skip item count header */
            }
        }

        for (int64_t i = 0; i < chunk_n && p + 4 <= end; i++) {
            uint32_t slen;
            memcpy(&slen, p, 4);
            p += 4;
            if (p + slen > end) break;
            strs[(size_t)(global_idx + i)] = (char *)malloc(slen + 1);
            if (!strs[(size_t)(global_idx + i)]) {
                free(raw);
                /* Return partial result — caller handles NULL entries */
                break;
            }
            memcpy(strs[(size_t)(global_idx + i)], p, slen);
            strs[(size_t)(global_idx + i)][slen] = '\0';
            p += slen;
        }

        free(raw);
        global_idx += chunk_n;
    }

    return strs;
}

/* ── Free string array ──────────────────────────────────────────────────────── */

static void free_strings(char **strs, int64_t n) {
    if (!strs) return;
    for (int64_t i = 0; i < n; i++) free(strs[(size_t)i]);
    free(strs);
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Zarr write helpers
 * ══════════════════════════════════════════════════════════════════════════════ */

static int zarr_create_group(const char *dir) {
    sc_mkdir_p(dir);
    return sc_json_write_zgroup(dir);
}

/* Write a 1D numeric array as a zarr array */
static int zarr_write_numeric_1d(const char *dir, const void *data,
                                  int64_t n, const char *dtype, int gzip) {
    sc_mkdir_p(dir);
    sc_zarr_meta_t meta = {0};
    meta.zarr_format = 2;
    meta.ndim = 1;
    meta.shape[0] = n;
    meta.chunks[0] = n;
    strncpy(meta.dtype, dtype, sizeof(meta.dtype) - 1);
    if (gzip > 0) {
        strcpy(meta.compressor_id, "zlib");
        meta.compressor_level = gzip;
    }
    meta.order = 'C';
    sc_json_write_zarray(dir, &meta);

    int elem_sz = dtype_size(dtype);
    char chunk_path[2048];
    snprintf(chunk_path, sizeof(chunk_path), "%s/0", dir);
    return zarr_write_chunk(chunk_path, data, (size_t)n * (size_t)elem_sz, gzip);
}

/* Write a 2D numeric array as a zarr array (C order) */
static int zarr_write_numeric_2d(const char *dir, const void *data,
                                  int64_t rows, int64_t cols,
                                  const char *dtype, int gzip) {
    sc_mkdir_p(dir);
    sc_zarr_meta_t meta = {0};
    meta.zarr_format = 2;
    meta.ndim = 2;
    meta.shape[0] = rows;
    meta.shape[1] = cols;
    meta.chunks[0] = rows;
    meta.chunks[1] = cols;
    strncpy(meta.dtype, dtype, sizeof(meta.dtype) - 1);
    if (gzip > 0) {
        strcpy(meta.compressor_id, "zlib");
        meta.compressor_level = gzip;
    }
    meta.order = 'C';
    sc_json_write_zarray(dir, &meta);

    int elem_sz = dtype_size(dtype);
    char chunk_path[2048];
    snprintf(chunk_path, sizeof(chunk_path), "%s/0.0", dir);
    return zarr_write_chunk(chunk_path, data, (size_t)(rows * cols) * elem_sz, gzip);
}

/* Write a vlen-utf8 string array */
static int zarr_write_strings(const char *dir, const char **strings,
                               int64_t n, int gzip) {
    sc_mkdir_p(dir);
    sc_zarr_meta_t meta = {0};
    meta.zarr_format = 2;
    meta.ndim = 1;
    meta.shape[0] = n;
    meta.chunks[0] = n;
    strcpy(meta.dtype, "|O");
    meta.has_vlen_utf8 = 1;
    if (gzip > 0) {
        strcpy(meta.compressor_id, "zlib");
        meta.compressor_level = gzip;
    }
    meta.order = 'C';
    sc_json_write_zarray(dir, &meta);

    /* Encode vlen-utf8: [uint32 count] [uint32 len][bytes] ... */
    size_t total = 4; /* header */
    for (int64_t i = 0; i < n; i++)
        total += 4 + (strings[(size_t)i] ? strlen(strings[(size_t)i]) : 0);

    uint8_t *buf = (uint8_t *)malloc(total);
    if (!buf) return SC_ERR;
    uint8_t *p = buf;
    uint32_t count = (uint32_t)n;
    memcpy(p, &count, 4); p += 4;
    for (int64_t i = 0; i < n; i++) {
        const char *s = strings[(size_t)i] ? strings[(size_t)i] : "";
        uint32_t slen = (uint32_t)strlen(s);
        memcpy(p, &slen, 4); p += 4;
        memcpy(p, s, slen); p += slen;
    }

    char chunk_path[2048];
    snprintf(chunk_path, sizeof(chunk_path), "%s/0", dir);
    int rc = zarr_write_chunk(chunk_path, buf, total, gzip);
    free(buf);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  zarr → h5seurat
 * ══════════════════════════════════════════════════════════════════════════════ */

/* Read a zarr sparse matrix group and write as h5seurat CSC group.
 * Zarr CSR (cells × genes) → h5seurat CSC (genes × cells): same arrays, swap dims.
 * Zarr CSC (genes × cells) → h5seurat CSC (genes × cells): direct copy.
 */
static int zarr_sparse_to_h5seurat(const char *zarr_dir, hid_t h5_grp,
                                    int64_t *out_nrows, int64_t *out_ncols,
                                    int gzip, const char *enc_type) {
    /* Read shape from .zattrs */
    char attrs_path[2048];
    snprintf(attrs_path, sizeof(attrs_path), "%s/.zattrs", zarr_dir);
    char *attrs_json = sc_json_read_file(attrs_path);
    int64_t shape[2] = {0, 0};
    if (attrs_json) {
        sc_json_get_int_array(attrs_json, "shape", shape, 2);
        free(attrs_json);
    }

    int is_csr = (enc_type && strcmp(enc_type, "csr_matrix") == 0);

    /* Read data, indices, indptr from zarr sub-arrays */
    char sub[2048];
    int64_t n_data, n_idx, n_ptr;

    snprintf(sub, sizeof(sub), "%s/data", zarr_dir);
    void *data_buf = zarr_read_numeric(sub, &n_data);
    snprintf(sub, sizeof(sub), "%s/indices", zarr_dir);
    void *idx_buf = zarr_read_numeric(sub, &n_idx);
    snprintf(sub, sizeof(sub), "%s/indptr", zarr_dir);
    void *ptr_buf = zarr_read_numeric(sub, &n_ptr);

    if (!data_buf || !idx_buf || !ptr_buf) {
        free(data_buf); free(idx_buf); free(ptr_buf);
        return SC_ERR;
    }

    /* Detect dtypes for correct interpretation */
    sc_zarr_meta_t data_meta, idx_meta, ptr_meta;
    snprintf(sub, sizeof(sub), "%s/data", zarr_dir);
    zarr_read_meta(sub, &data_meta);
    snprintf(sub, sizeof(sub), "%s/indices", zarr_dir);
    zarr_read_meta(sub, &idx_meta);
    snprintf(sub, sizeof(sub), "%s/indptr", zarr_dir);
    zarr_read_meta(sub, &ptr_meta);

    /* Write to HDF5 */
    int rc = SC_OK;
    hsize_t nnz = (hsize_t)n_data;

    /* data */
    {
        hsize_t dim = nnz;
        hid_t sp = H5Screate_simple(1, &dim, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = dim < 65536 ? dim : 65536;
        if (chunk > 0) H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip > 0 && chunk > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
        hid_t ds = H5Dcreate2(h5_grp, "data", H5T_NATIVE_DOUBLE, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        if (strcmp(data_meta.dtype, "<f4") == 0) {
            /* float32 -> float64 conversion */
            double *d64 = (double *)malloc(nnz * sizeof(double));
            if (d64) {
                float *f32 = (float *)data_buf;
                for (hsize_t i = 0; i < nnz; i++) d64[i] = (double)f32[i];
                H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d64);
                free(d64);
            }
        } else {
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data_buf);
        }
        H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
    }

    /* indices */
    {
        hsize_t dim = (hsize_t)n_idx;
        hid_t sp = H5Screate_simple(1, &dim, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = dim < 65536 ? dim : 65536;
        if (chunk > 0) H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip > 0 && chunk > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
        hid_t ds = H5Dcreate2(h5_grp, "indices", H5T_NATIVE_INT32, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        if (dtype_size(idx_meta.dtype) == 8) {
            int32_t *i32 = (int32_t *)malloc(n_idx * sizeof(int32_t));
            if (i32) {
                int64_t *i64 = (int64_t *)idx_buf;
                for (int64_t k = 0; k < n_idx; k++) i32[k] = (int32_t)i64[k];
                H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, i32);
                free(i32);
            }
        } else {
            H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, idx_buf);
        }
        H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
    }

    /* indptr */
    {
        hsize_t dim = (hsize_t)n_ptr;
        hid_t sp = H5Screate_simple(1, &dim, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk = dim < 65536 ? dim : 65536;
        if (chunk > 0) H5Pset_chunk(dcpl, 1, &chunk);
        if (gzip > 0 && chunk > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
        hid_t ds = H5Dcreate2(h5_grp, "indptr", H5T_NATIVE_INT64, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        if (dtype_size(ptr_meta.dtype) == 4) {
            int64_t *i64 = (int64_t *)malloc(n_ptr * sizeof(int64_t));
            if (i64) {
                int32_t *i32 = (int32_t *)ptr_buf;
                for (int64_t k = 0; k < n_ptr; k++) i64[k] = (int64_t)i32[k];
                H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, i64);
                free(i64);
            }
        } else {
            H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptr_buf);
        }
        H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
    }

    /* dims attr: h5seurat stores [n_genes, n_cells] */
    if (is_csr) {
        /* CSR(cells × genes) → CSC(genes × cells): swap dims */
        int32_t dims[2] = {(int32_t)shape[1], (int32_t)shape[0]};
        hsize_t two = 2;
        hid_t sp = H5Screate_simple(1, &two, NULL);
        hid_t attr = H5Acreate2(h5_grp, "dims", H5T_NATIVE_INT32, sp,
                                  H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT32, dims);
        H5Aclose(attr); H5Sclose(sp);
        if (out_nrows) *out_nrows = shape[1]; /* genes */
        if (out_ncols) *out_ncols = shape[0]; /* cells */
    } else {
        /* CSC(genes × cells): same orientation */
        int32_t dims[2] = {(int32_t)shape[0], (int32_t)shape[1]};
        hsize_t two = 2;
        hid_t sp = H5Screate_simple(1, &two, NULL);
        hid_t attr = H5Acreate2(h5_grp, "dims", H5T_NATIVE_INT32, sp,
                                  H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT32, dims);
        H5Aclose(attr); H5Sclose(sp);
        if (out_nrows) *out_nrows = shape[0];
        if (out_ncols) *out_ncols = shape[1];
    }

    free(data_buf); free(idx_buf); free(ptr_buf);
    return rc;
}

/* Read zarr categorical group → h5seurat factor group (0-based→1-based) */
static int zarr_categorical_to_h5seurat(const char *cat_dir, hid_t h5_grp,
                                          int gzip) {
    char sub[2048];

    /* Read categories (strings) */
    snprintf(sub, sizeof(sub), "%s/categories", cat_dir);
    int64_t n_cats;
    char **cats = zarr_read_strings(sub, &n_cats);
    if (!cats) return SC_ERR;

    /* Read codes (int8 array, 0-based) */
    snprintf(sub, sizeof(sub), "%s/codes", cat_dir);
    int64_t n_codes;
    void *codes_raw = zarr_read_numeric(sub, &n_codes);
    if (!codes_raw) { free_strings(cats, n_cats); return SC_ERR; }

    /* Detect codes dtype */
    sc_zarr_meta_t codes_meta;
    zarr_read_meta(sub, &codes_meta);
    int codes_elem = dtype_size(codes_meta.dtype);

    /* Write levels (string array) */
    {
        hid_t str_type = sc_create_vlen_str_type();
        hsize_t dim = (hsize_t)n_cats;
        hid_t sp = H5Screate_simple(1, &dim, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        if (dim > 0) {
            hsize_t chunk = dim < 65536 ? dim : 65536;
            H5Pset_chunk(dcpl, 1, &chunk);
            if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
        }
        hid_t ds = H5Dcreate2(h5_grp, "levels", str_type, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, cats);
        H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp); H5Tclose(str_type);
    }

    /* Write values (int32 array, 1-based) */
    {
        int32_t *vals = (int32_t *)malloc((size_t)n_codes * sizeof(int32_t));
        if (!vals) { free_strings(cats, n_cats); free(codes_raw); return SC_ERR; }
        for (int64_t i = 0; i < n_codes; i++) {
            int32_t code;
            if (codes_elem == 1)
                code = (int32_t)((int8_t *)codes_raw)[i];
            else if (codes_elem == 2)
                code = (int32_t)((int16_t *)codes_raw)[i];
            else
                code = ((int32_t *)codes_raw)[i];
            vals[i] = (code < 0) ? 0 : code + 1; /* 0-based → 1-based */
        }
        hsize_t dim = (hsize_t)n_codes;
        hid_t sp = H5Screate_simple(1, &dim, NULL);
        hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
        if (dim > 0) {
            hsize_t chunk = dim < 65536 ? dim : 65536;
            H5Pset_chunk(dcpl, 1, &chunk);
            if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
        }
        hid_t ds = H5Dcreate2(h5_grp, "values", H5T_NATIVE_INT32, sp,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
        H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
        free(vals);
    }

    free_strings(cats, n_cats);
    free(codes_raw);
    return SC_OK;
}

/* Write zarr obs/var DataFrame → h5seurat meta.data / var */
static int zarr_df_to_h5seurat_md(const char *zarr_df_dir, hid_t h5_md_grp,
                                    int gzip) {
    /* Read .zattrs for column-order */
    char attrs_path[2048];
    snprintf(attrs_path, sizeof(attrs_path), "%s/.zattrs", zarr_df_dir);
    char *attrs_json = sc_json_read_file(attrs_path);
    if (!attrs_json) return SC_ERR_IO;

    char **col_order = NULL;
    int n_cols = 0;
    sc_json_get_str_array(attrs_json, "column-order", &col_order, &n_cols);
    free(attrs_json);

    for (int ci = 0; ci < n_cols; ci++) {
        const char *col_name = col_order[ci];
        char col_dir[2048];
        snprintf(col_dir, sizeof(col_dir), "%s/%s", zarr_df_dir, col_name);

        /* Determine if categorical (has .zgroup and codes/ subdir) */
        char test[2048];
        snprintf(test, sizeof(test), "%s/.zgroup", col_dir);
        if (access(test, F_OK) == 0) {
            /* Categorical group */
            hid_t factor_grp = H5Gcreate2(h5_md_grp, col_name, H5P_DEFAULT,
                                            H5P_DEFAULT, H5P_DEFAULT);
            zarr_categorical_to_h5seurat(col_dir, factor_grp, gzip);
            H5Gclose(factor_grp);
            continue;
        }

        /* Regular array */
        snprintf(test, sizeof(test), "%s/.zarray", col_dir);
        if (access(test, F_OK) != 0) continue;

        sc_zarr_meta_t col_meta;
        zarr_read_meta(col_dir, &col_meta);

        if (col_meta.has_vlen_utf8 || strcmp(col_meta.dtype, "|O") == 0) {
            /* String array */
            int64_t n_str;
            char **strs = zarr_read_strings(col_dir, &n_str);
            if (strs) {
                hid_t str_type = sc_create_vlen_str_type();
                hsize_t dim = (hsize_t)n_str;
                hid_t sp = H5Screate_simple(1, &dim, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                if (dim > 0) {
                    hsize_t chunk = dim < 65536 ? dim : 65536;
                    H5Pset_chunk(dcpl, 1, &chunk);
                    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
                }
                hid_t ds = H5Dcreate2(h5_md_grp, col_name, str_type, sp,
                                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
                H5Dwrite(ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);
                H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp); H5Tclose(str_type);
                free_strings(strs, n_str);
            }
        } else if (strcmp(col_meta.dtype, "|b1") == 0) {
            /* Boolean -> int32 */
            int64_t nn;
            void *raw = zarr_read_numeric(col_dir, &nn);
            if (raw) {
                int32_t *ivals = (int32_t *)malloc((size_t)nn * sizeof(int32_t));
                if (!ivals) { free(raw); sc_json_free_str_array(col_order, n_cols); return SC_ERR; }
                int8_t *b = (int8_t *)raw;
                for (int64_t k = 0; k < nn; k++) ivals[k] = (int32_t)b[k];
                hsize_t dim = (hsize_t)nn;
                hid_t sp = H5Screate_simple(1, &dim, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                if (dim > 0) {
                    hsize_t chunk = dim < 65536 ? dim : 65536;
                    H5Pset_chunk(dcpl, 1, &chunk);
                    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
                }
                hid_t ds = H5Dcreate2(h5_md_grp, col_name, H5T_NATIVE_INT32, sp,
                                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
                H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, ivals);
                H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
                free(ivals); free(raw);
            }
        } else {
            /* Numeric (float64 or int) */
            int64_t nn;
            void *raw = zarr_read_numeric(col_dir, &nn);
            if (raw) {
                hsize_t dim = (hsize_t)nn;
                hid_t sp = H5Screate_simple(1, &dim, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                if (dim > 0) {
                    hsize_t chunk = dim < 65536 ? dim : 65536;
                    H5Pset_chunk(dcpl, 1, &chunk);
                    if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
                }
                hid_t h5type = (col_meta.dtype[1] == 'i') ?
                                H5T_NATIVE_INT32 : H5T_NATIVE_DOUBLE;
                hid_t ds = H5Dcreate2(h5_md_grp, col_name, h5type, sp,
                                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
                H5Dwrite(ds, h5type, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw);
                H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
                free(raw);
            }
        }
    }

    /* Set colnames attribute */
    if (n_cols > 0) {
        sc_set_str_array_attr(h5_md_grp, "colnames",
                               (const char **)col_order, (hsize_t)n_cols);
    }

    sc_json_free_str_array(col_order, n_cols);
    return SC_OK;
}

/* ── Main zarr → h5seurat ──────────────────────────────────────────────────── */

int sc_zarr_to_h5seurat(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (zarr → h5seurat)\n",
                opts->input_path, opts->output_path);

    const char *store = opts->input_path;
    const char *assay = opts->assay_name;
    int gzip = opts->gzip_level;

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC,
                            H5P_DEFAULT, H5P_DEFAULT);
    if (dst < 0) return SC_ERR_HDF;

    int64_t n_cells = 0, n_genes = 0;

    /* ── 1. X → assays/{assay}/layers/data ─────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [1/6] Transferring X...\n");
    {
        char assay_path[512];
        snprintf(assay_path, sizeof(assay_path), "assays/%s", assay);
        hid_t assays_grp = H5Gcreate2(dst, "assays", H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
        hid_t a_grp = H5Gcreate2(assays_grp, assay, H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);
        sc_set_str_attr(a_grp, "key", "rna_");
        sc_set_str_attr(a_grp, "s4class", "SeuratObject::Assay5");
        hid_t layers = H5Gcreate2(a_grp, "layers", H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);

        /* Read encoding type of X */
        char x_enc[64] = "csr_matrix";
        char x_attrs_path[2048];
        snprintf(x_attrs_path, sizeof(x_attrs_path), "%s/X/.zattrs", store);
        char *x_attrs = sc_json_read_file(x_attrs_path);
        if (x_attrs) {
            sc_json_get_string(x_attrs, "encoding-type", x_enc, sizeof(x_enc));
            free(x_attrs);
        }

        char x_dir[2048];
        snprintf(x_dir, sizeof(x_dir), "%s/X", store);

        if (zarr_is_group(store, "X") &&
            (strcmp(x_enc, "csr_matrix") == 0 || strcmp(x_enc, "csc_matrix") == 0)) {
            hid_t data_grp = H5Gcreate2(layers, "data", H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
            zarr_sparse_to_h5seurat(x_dir, data_grp, &n_genes, &n_cells,
                                     gzip, x_enc);
            H5Gclose(data_grp);
        } else if (zarr_is_array(store, "X")) {
            /* Dense X */
            int64_t nn;
            void *data = zarr_read_numeric(x_dir, &nn);
            if (data) {
                sc_zarr_meta_t xm;
                zarr_read_meta(x_dir, &xm);
                n_cells = xm.shape[0];
                n_genes = xm.shape[1];
                /* TODO: store dense or convert to sparse */
                free(data);
            }
        }

        /* raw/X → layers/counts */
        char raw_x[2048];
        snprintf(raw_x, sizeof(raw_x), "%s/raw/X", store);
        if (zarr_exists(store, "raw/X")) {
            if (opts->verbose) SC_MSG("  [1b] Transferring raw/X...\n");
            char raw_enc[64] = "csr_matrix";
            char raw_attrs_path[2048];
            snprintf(raw_attrs_path, sizeof(raw_attrs_path), "%s/.zattrs", raw_x);
            char *ra = sc_json_read_file(raw_attrs_path);
            if (ra) {
                sc_json_get_string(ra, "encoding-type", raw_enc, sizeof(raw_enc));
                free(ra);
            }
            hid_t counts_grp = H5Gcreate2(layers, "counts", H5P_DEFAULT,
                                            H5P_DEFAULT, H5P_DEFAULT);
            zarr_sparse_to_h5seurat(raw_x, counts_grp, NULL, NULL, gzip, raw_enc);
            H5Gclose(counts_grp);
        }

        /* features (var/_index) */
        char var_idx[2048];
        snprintf(var_idx, sizeof(var_idx), "%s/var/_index", store);
        if (access(var_idx, F_OK) == 0 || zarr_is_array(store, "var/_index")) {
            int64_t nf;
            char **feats = zarr_read_strings(var_idx, &nf);
            if (feats) {
                hid_t str_type = sc_create_vlen_str_type();
                hsize_t dim = (hsize_t)nf;
                hid_t sp = H5Screate_simple(1, &dim, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                if (dim > 0) {
                    hsize_t chunk = dim < 65536 ? dim : 65536;
                    H5Pset_chunk(dcpl, 1, &chunk);
                }
                hid_t ds = H5Dcreate2(a_grp, "features", str_type, sp,
                                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
                H5Dwrite(ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, feats);
                H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp); H5Tclose(str_type);
                if (n_genes == 0) n_genes = nf;
                free_strings(feats, nf);
            }
        }

        H5Gclose(layers);
        H5Gclose(a_grp);
        H5Gclose(assays_grp);
    }

    /* ── 2. obs → cell.names + meta.data ───────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [2/6] Transferring obs...\n");
    {
        /* Read cell names */
        char idx_dir[2048];
        snprintf(idx_dir, sizeof(idx_dir), "%s/obs/_index", store);
        int64_t nc;
        char **cells = zarr_read_strings(idx_dir, &nc);
        if (cells) {
            if (n_cells == 0) n_cells = nc;
            hid_t str_type = sc_create_vlen_str_type();
            hsize_t dim = (hsize_t)nc;
            hid_t sp = H5Screate_simple(1, &dim, NULL);
            hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
            if (dim > 0) {
                hsize_t chunk = dim < 65536 ? dim : 65536;
                H5Pset_chunk(dcpl, 1, &chunk);
            }
            hid_t ds = H5Dcreate2(dst, "cell.names", str_type, sp,
                                   H5P_DEFAULT, dcpl, H5P_DEFAULT);
            H5Dwrite(ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, cells);
            H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp); H5Tclose(str_type);
            free_strings(cells, nc);
        }

        /* meta.data columns */
        char obs_dir[2048];
        snprintf(obs_dir, sizeof(obs_dir), "%s/obs", store);
        hid_t md_grp = H5Gcreate2(dst, "meta.data", H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);
        zarr_df_to_h5seurat_md(obs_dir, md_grp, gzip);
        H5Gclose(md_grp);
    }

    /* ── 3. var metadata ───────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [3/6] Transferring var...\n");
    /* Variable features etc. are stored in var columns — skip for now,
     * they'll be in the assay already from features above */

    /* ── 4. obsm → reductions ──────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [4/6] Transferring obsm...\n");
    {
        hid_t red_grp = H5Gcreate2(dst, "reductions", H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
        char obsm_dir[2048];
        snprintf(obsm_dir, sizeof(obsm_dir), "%s/obsm", store);
        if (zarr_exists(store, "obsm")) {
            char **children;
            int n_children;
            zarr_list_children(obsm_dir, &children, &n_children);
            for (int ci = 0; ci < n_children; ci++) {
                char emb_dir[2048];
                snprintf(emb_dir, sizeof(emb_dir), "%s/%s", obsm_dir, children[ci]);
                if (!zarr_is_array(obsm_dir, children[ci])) continue;

                sc_zarr_meta_t emb_meta;
                zarr_read_meta(emb_dir, &emb_meta);
                if (emb_meta.ndim != 2) continue;

                int64_t nn;
                void *flat = zarr_read_numeric(emb_dir, &nn);
                if (!flat) continue;

                /* Strip X_ prefix for reduction name */
                const char *rname = children[ci];
                if (rname[0] == 'X' && rname[1] == '_') rname += 2;

                /* Transpose: zarr C [n_cells, n_comp] → h5seurat [n_comp, n_cells] */
                int64_t nr = emb_meta.shape[0], nc_e = emb_meta.shape[1];
                double *transposed = (double *)malloc((size_t)(nr * nc_e) * sizeof(double));
                if (!transposed) { free(flat); continue; }
                double *src = (double *)flat;
                for (int64_t r = 0; r < nr; r++)
                    for (int64_t c = 0; c < nc_e; c++)
                        transposed[c * nr + r] = src[r * nc_e + c];

                hid_t r_grp = H5Gcreate2(red_grp, rname, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
                sc_set_str_attr(r_grp, "active.assay", assay);
                char key_buf[64];
                snprintf(key_buf, sizeof(key_buf), "%s_", rname);
                sc_set_str_attr(r_grp, "key", key_buf);
                {
                    int32_t global = 0;
                    hid_t sp = H5Screate(H5S_SCALAR);
                    hid_t attr = H5Acreate2(r_grp, "global", H5T_NATIVE_INT32,
                                             sp, H5P_DEFAULT, H5P_DEFAULT);
                    H5Awrite(attr, H5T_NATIVE_INT32, &global);
                    H5Aclose(attr); H5Sclose(sp);
                }
                { hid_t g = H5Gcreate2(r_grp, "misc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                  if (g >= 0) H5Gclose(g); }

                /* Write cell.embeddings [n_comp, n_cells] */
                hsize_t dims[2] = {(hsize_t)nc_e, (hsize_t)nr};
                hid_t sp = H5Screate_simple(2, dims, NULL);
                hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl, 2, dims);
                if (gzip > 0) H5Pset_deflate(dcpl, (unsigned)gzip);
                hid_t ds = H5Dcreate2(r_grp, "cell.embeddings",
                                       H5T_NATIVE_DOUBLE, sp,
                                       H5P_DEFAULT, dcpl, H5P_DEFAULT);
                H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, transposed);
                H5Dclose(ds); H5Pclose(dcpl); H5Sclose(sp);
                H5Gclose(r_grp);

                free(transposed);
                free(flat);
            }
            zarr_free_children(children, n_children);
        }
        H5Gclose(red_grp);
    }

    /* ── 5. obsp → graphs ──────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [5/6] Transferring obsp...\n");
    {
        hid_t graphs_grp = H5Gcreate2(dst, "graphs", H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT);
        char obsp_dir[2048];
        snprintf(obsp_dir, sizeof(obsp_dir), "%s/obsp", store);
        if (zarr_exists(store, "obsp")) {
            char **children;
            int n_children;
            zarr_list_children(obsp_dir, &children, &n_children);
            for (int ci = 0; ci < n_children; ci++) {
                char g_dir[2048];
                snprintf(g_dir, sizeof(g_dir), "%s/%s", obsp_dir, children[ci]);
                if (!zarr_is_group(obsp_dir, children[ci])) continue;

                /* Read encoding type */
                char g_attrs[2048];
                snprintf(g_attrs, sizeof(g_attrs), "%s/.zattrs", g_dir);
                char *gaj = sc_json_read_file(g_attrs);
                char enc[64] = "csr_matrix";
                if (gaj) {
                    sc_json_get_string(gaj, "encoding-type", enc, sizeof(enc));
                    free(gaj);
                }

                /* Add assay prefix for h5seurat convention */
                char gname[512];
                snprintf(gname, sizeof(gname), "%s_%s", assay, children[ci]);

                hid_t g_grp = H5Gcreate2(graphs_grp, gname, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
                zarr_sparse_to_h5seurat(g_dir, g_grp, NULL, NULL, gzip, enc);
                sc_set_str_attr(g_grp, "assay.used", assay);
                H5Gclose(g_grp);
            }
            zarr_free_children(children, n_children);
        }
        H5Gclose(graphs_grp);
    }

    /* ── 6. scaffold ───────────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [6/6] Writing scaffold...\n");
    {
        sc_set_str_attr(dst, "active.assay", assay);
        sc_set_str_attr(dst, "project", "SeuratProject");
        sc_set_str_attr(dst, "version", "5.3.0");

        /* Required empty groups */
        const char *empty_grps[] = {"tools", "commands", "images", "neighbors"};
        for (int i = 0; i < 4; i++) {
            if (!sc_has_group(dst, empty_grps[i])) {
                hid_t g = H5Gcreate2(dst, empty_grps[i], H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
                if (g >= 0) H5Gclose(g);
            }
        }
        if (!sc_has_group(dst, "misc")) {
            hid_t g = H5Gcreate2(dst, "misc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (g >= 0) H5Gclose(g);
        }

        /* active.ident */
        if (n_cells > 0) {
            hid_t ident_grp = H5Gcreate2(dst, "active.ident", H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
            const char *lvl = "SeuratProject";
            hid_t str_type = sc_create_vlen_str_type();
            hsize_t one = 1;
            hid_t sp2 = H5Screate_simple(1, &one, NULL);
            hid_t ds = H5Dcreate2(ident_grp, "levels", str_type, sp2,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &lvl);
            H5Dclose(ds); H5Sclose(sp2);

            hsize_t nc_dim = (hsize_t)n_cells;
            int32_t *ones = (int32_t *)calloc((size_t)n_cells, sizeof(int32_t));
            if (!ones) { H5Tclose(str_type); H5Gclose(ident_grp); H5Fclose(dst); return SC_ERR; }
            for (int64_t i = 0; i < n_cells; i++) ones[i] = 1;
            sp2 = H5Screate_simple(1, &nc_dim, NULL);
            ds = H5Dcreate2(ident_grp, "values", H5T_NATIVE_INT32, sp2,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, ones);
            H5Dclose(ds); H5Sclose(sp2);
            H5Tclose(str_type);
            H5Gclose(ident_grp);
            free(ones);
        }
    }

    H5Fclose(dst);
    if (opts->verbose) SC_MSG("[scConvert] Done.\n");
    return SC_OK;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  h5seurat → zarr
 * ══════════════════════════════════════════════════════════════════════════════ */

/* Write h5seurat CSC sparse group → zarr CSR (swap shape = zero-copy trick) */
static int h5seurat_sparse_to_zarr(hid_t h5_grp, const char *zarr_dir,
                                     int gzip, const char *enc_type) {
    zarr_create_group(zarr_dir);

    /* Read dims attr */
    int32_t dims[2] = {0, 0};
    if (H5Aexists(h5_grp, "dims")) {
        hid_t attr = H5Aopen(h5_grp, "dims", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_INT32, dims);
        H5Aclose(attr);
    }

    int is_csr = (strcmp(enc_type, "csr_matrix") == 0);

    /* data */
    if (sc_has_dataset(h5_grp, "data")) {
        hid_t ds = H5Dopen2(h5_grp, "data", H5P_DEFAULT);
        hid_t sp = H5Dget_space(ds);
        hsize_t nnz;
        H5Sget_simple_extent_dims(sp, &nnz, NULL);
        H5Sclose(sp);

        double *buf = (double *)malloc(nnz * sizeof(double));
        if (!buf) { H5Dclose(ds); return SC_ERR; }
        H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dclose(ds);

        char dir[2048];
        snprintf(dir, sizeof(dir), "%s/data", zarr_dir);
        zarr_write_numeric_1d(dir, buf, (int64_t)nnz, "<f8", gzip);
        free(buf);
    }

    /* indices */
    if (sc_has_dataset(h5_grp, "indices")) {
        hid_t ds = H5Dopen2(h5_grp, "indices", H5P_DEFAULT);
        hid_t sp = H5Dget_space(ds);
        hsize_t n;
        H5Sget_simple_extent_dims(sp, &n, NULL);
        H5Sclose(sp);

        int32_t *buf = (int32_t *)malloc(n * sizeof(int32_t));
        if (!buf) { H5Dclose(ds); return SC_ERR; }
        H5Dread(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dclose(ds);

        char dir[2048];
        snprintf(dir, sizeof(dir), "%s/indices", zarr_dir);
        zarr_write_numeric_1d(dir, buf, (int64_t)n, "<i4", gzip);
        free(buf);
    }

    /* indptr */
    if (sc_has_dataset(h5_grp, "indptr")) {
        hid_t ds = H5Dopen2(h5_grp, "indptr", H5P_DEFAULT);
        hid_t sp = H5Dget_space(ds);
        hsize_t n;
        H5Sget_simple_extent_dims(sp, &n, NULL);
        H5Sclose(sp);

        int64_t *buf = (int64_t *)malloc(n * sizeof(int64_t));
        if (!buf) { H5Dclose(ds); return SC_ERR; }
        H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        H5Dclose(ds);

        char dir[2048];
        snprintf(dir, sizeof(dir), "%s/indptr", zarr_dir);
        zarr_write_numeric_1d(dir, buf, (int64_t)n, "<i8", gzip);
        free(buf);
    }

    /* .zattrs with encoding-type and shape */
    char attrs_buf[512];
    if (is_csr) {
        /* CSC [n_genes, n_cells] → CSR shape [n_cells, n_genes] */
        snprintf(attrs_buf, sizeof(attrs_buf),
                 "{\"encoding-type\": \"%s\", \"encoding-version\": \"0.1.0\", "
                 "\"shape\": [%d, %d]}",
                 enc_type, dims[1], dims[0]);
    } else {
        snprintf(attrs_buf, sizeof(attrs_buf),
                 "{\"encoding-type\": \"%s\", \"encoding-version\": \"0.1.0\", "
                 "\"shape\": [%d, %d]}",
                 enc_type, dims[0], dims[1]);
    }
    sc_json_write_zattrs_str(zarr_dir, attrs_buf);
    return SC_OK;
}

/* Write h5seurat factor → zarr categorical */
static int h5seurat_factor_to_zarr(hid_t h5_grp, const char *cat_dir, int gzip) {
    zarr_create_group(cat_dir);

    /* Read levels */
    hid_t lev_ds = H5Dopen2(h5_grp, "levels", H5P_DEFAULT);
    hid_t lev_sp = H5Dget_space(lev_ds);
    hsize_t n_lev;
    H5Sget_simple_extent_dims(lev_sp, &n_lev, NULL);
    H5Sclose(lev_sp);

    hid_t str_type = sc_create_vlen_str_type();
    char **levels = (char **)calloc(n_lev, sizeof(char *));
    if (!levels) { H5Tclose(str_type); H5Dclose(lev_ds); return SC_ERR; }
    H5Dread(lev_ds, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, levels);
    H5Dclose(lev_ds);

    /* Write categories */
    char cat_path[2048];
    snprintf(cat_path, sizeof(cat_path), "%s/categories", cat_dir);
    zarr_write_strings(cat_path, (const char **)levels, (int64_t)n_lev, gzip);

    /* Read values (1-based int32) */
    hid_t val_ds = H5Dopen2(h5_grp, "values", H5P_DEFAULT);
    hid_t val_sp = H5Dget_space(val_ds);
    hsize_t n_val;
    H5Sget_simple_extent_dims(val_sp, &n_val, NULL);
    H5Sclose(val_sp);

    int32_t *vals = (int32_t *)malloc(n_val * sizeof(int32_t));
    if (!vals) {
        for (hsize_t i = 0; i < n_lev; i++) free(levels[i]);
        free(levels);
        H5Tclose(str_type);
        H5Dclose(val_ds);
        return SC_ERR;
    }
    H5Dread(val_ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
    H5Dclose(val_ds);

    /* Convert 1-based → 0-based int8 */
    int8_t *codes = (int8_t *)malloc(n_val * sizeof(int8_t));
    if (!codes) {
        for (hsize_t i = 0; i < n_lev; i++) free(levels[i]);
        free(levels);
        H5Tclose(str_type);
        free(vals);
        return SC_ERR;
    }
    for (hsize_t i = 0; i < n_val; i++)
        codes[i] = (vals[i] <= 0) ? -1 : (int8_t)(vals[i] - 1);

    /* Write codes */
    char codes_path[2048];
    snprintf(codes_path, sizeof(codes_path), "%s/codes", cat_dir);
    zarr_write_numeric_1d(codes_path, codes, (int64_t)n_val, "|i1", gzip);

    /* .zattrs */
    sc_json_write_zattrs_str(cat_dir,
        "{\"encoding-type\": \"categorical\", \"encoding-version\": \"0.2.0\", "
        "\"ordered\": false}");

    /* Cleanup */
    for (hsize_t i = 0; i < n_lev; i++) free(levels[i]);
    free(levels);
    H5Tclose(str_type);
    free(vals);
    free(codes);
    return SC_OK;
}

/* Write h5seurat meta.data → zarr obs DataFrame */
static int h5seurat_md_to_zarr_df(hid_t h5_grp, const char *zarr_dir,
                                    const char **cell_names, int64_t n_cells,
                                    int gzip) {
    zarr_create_group(zarr_dir);

    /* Write _index */
    char idx_dir[2048];
    snprintf(idx_dir, sizeof(idx_dir), "%s/_index", zarr_dir);
    zarr_write_strings(idx_dir, cell_names, n_cells, gzip);

    /* Get column names from colnames attr */
    char **col_names = NULL;
    hsize_t n_cols = 0;
    if (H5Aexists(h5_grp, "colnames")) {
        hid_t attr = H5Aopen(h5_grp, "colnames", H5P_DEFAULT);
        hid_t sp = H5Aget_space(attr);
        H5Sget_simple_extent_dims(sp, &n_cols, NULL);
        H5Sclose(sp);
        hid_t stype = sc_create_vlen_str_type();
        col_names = (char **)calloc(n_cols, sizeof(char *));
        if (!col_names) { H5Aclose(attr); H5Tclose(stype); return SC_ERR; }
        H5Aread(attr, stype, col_names);
        H5Aclose(attr);
        H5Tclose(stype);
    }

    /* Write each column */
    for (hsize_t ci = 0; ci < n_cols; ci++) {
        const char *cname = col_names[ci];
        if (!cname) continue;

        if (sc_has_group(h5_grp, cname)) {
            /* Factor group */
            hid_t fg = H5Gopen2(h5_grp, cname, H5P_DEFAULT);
            char col_dir[2048];
            snprintf(col_dir, sizeof(col_dir), "%s/%s", zarr_dir, cname);
            h5seurat_factor_to_zarr(fg, col_dir, gzip);
            H5Gclose(fg);
        } else if (sc_has_dataset(h5_grp, cname)) {
            hid_t ds = H5Dopen2(h5_grp, cname, H5P_DEFAULT);
            hid_t dt = H5Dget_type(ds);
            hid_t sp = H5Dget_space(ds);
            hsize_t nn;
            H5Sget_simple_extent_dims(sp, &nn, NULL);
            H5Sclose(sp);

            char col_dir[2048];
            snprintf(col_dir, sizeof(col_dir), "%s/%s", zarr_dir, cname);

            if (H5Tget_class(dt) == H5T_STRING) {
                hid_t stype = sc_create_vlen_str_type();
                char **strs = (char **)calloc(nn, sizeof(char *));
                if (strs) {
                    H5Dread(ds, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, strs);
                    zarr_write_strings(col_dir, (const char **)strs, (int64_t)nn, gzip);
                    for (hsize_t k = 0; k < nn; k++) free(strs[k]);
                    free(strs);
                }
                H5Tclose(stype);
            } else if (H5Tget_class(dt) == H5T_FLOAT) {
                double *buf = (double *)malloc(nn * sizeof(double));
                if (buf) {
                    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
                    zarr_write_numeric_1d(col_dir, buf, (int64_t)nn, "<f8", gzip);
                    free(buf);
                }
            } else {
                int32_t *buf = (int32_t *)malloc(nn * sizeof(int32_t));
                if (buf) {
                    H5Dread(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
                    zarr_write_numeric_1d(col_dir, buf, (int64_t)nn, "<i4", gzip);
                    free(buf);
                }
            }
            H5Tclose(dt);
            H5Dclose(ds);
        }
    }

    /* Write .zattrs */
    /* Build column-order JSON array */
    size_t json_sz = 256;
    for (hsize_t i = 0; i < n_cols; i++)
        json_sz += strlen(col_names[i]) + 6;
    char *json = (char *)malloc(json_sz);
    if (!json) {
        for (hsize_t i = 0; i < n_cols; i++) free(col_names[i]);
        free(col_names);
        return SC_ERR;
    }
    char *jp = json;
    size_t rem = json_sz;
    int w;
    w = snprintf(jp, rem, "{\"encoding-type\": \"dataframe\", "
                           "\"encoding-version\": \"0.2.0\", "
                           "\"_index\": \"_index\", \"column-order\": [");
    jp += w; rem -= (size_t)w;
    for (hsize_t i = 0; i < n_cols; i++) {
        if (i > 0) { w = snprintf(jp, rem, ", "); jp += w; rem -= (size_t)w; }
        w = snprintf(jp, rem, "\"%s\"", col_names[i]); jp += w; rem -= (size_t)w;
    }
    snprintf(jp, rem, "]}");
    sc_json_write_zattrs_str(zarr_dir, json);
    free(json);

    /* Cleanup */
    for (hsize_t i = 0; i < n_cols; i++) free(col_names[i]);
    free(col_names);

    return SC_OK;
}

/* ── Main h5seurat → zarr ──────────────────────────────────────────────────── */

int sc_h5seurat_to_zarr(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5seurat → zarr)\n",
                opts->input_path, opts->output_path);

    const char *assay = opts->assay_name;
    int gzip = opts->gzip_level;
    const char *out = opts->output_path;

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) return SC_ERR_HDF;

    /* Detect active assay */
    char active_assay[256] = "RNA";
    sc_get_str_attr(src, "active.assay", active_assay, sizeof(active_assay));
    if (strcmp(assay, "RNA") == 0 && strcmp(active_assay, "RNA") != 0)
        assay = active_assay;

    /* Create root zarr group */
    sc_mkdir_p(out);
    sc_json_write_zgroup(out);
    sc_json_write_zattrs_str(out,
        "{\"encoding-type\": \"anndata\", \"encoding-version\": \"0.1.0\"}");

    int64_t n_cells = 0, n_genes = 0;

    /* ── 1. X (layers/data → X as CSR) ─────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [1/6] Transferring X...\n");
    {
        char data_path[256];
        snprintf(data_path, sizeof(data_path), "assays/%s/layers/data", assay);
        if (!sc_has_group(src, data_path))
            snprintf(data_path, sizeof(data_path), "assays/%s/data", assay);

        if (sc_has_group(src, data_path)) {
            hid_t grp = H5Gopen2(src, data_path, H5P_DEFAULT);
            sc_csr_info_t info = {0};
            sc_read_csr_info(grp, &info);
            n_genes = (int64_t)info.n_rows;
            n_cells = (int64_t)info.n_cols;

            char x_dir[2048];
            snprintf(x_dir, sizeof(x_dir), "%s/X", out);
            h5seurat_sparse_to_zarr(grp, x_dir, gzip, "csr_matrix");
            H5Gclose(grp);
        }

        /* raw/X (layers/counts) */
        char counts_path[256];
        snprintf(counts_path, sizeof(counts_path), "assays/%s/layers/counts", assay);
        if (!sc_has_group(src, counts_path))
            snprintf(counts_path, sizeof(counts_path), "assays/%s/counts", assay);

        if (sc_has_group(src, counts_path)) {
            if (opts->verbose) SC_MSG("  [1b] Transferring raw/X...\n");
            hid_t grp = H5Gopen2(src, counts_path, H5P_DEFAULT);

            char raw_dir[2048], raw_x_dir[2048];
            snprintf(raw_dir, sizeof(raw_dir), "%s/raw", out);
            zarr_create_group(raw_dir);
            sc_json_write_zattrs_str(raw_dir,
                "{\"encoding-type\": \"raw\", \"encoding-version\": \"0.1.0\"}");

            snprintf(raw_x_dir, sizeof(raw_x_dir), "%s/raw/X", out);
            h5seurat_sparse_to_zarr(grp, raw_x_dir, gzip, "csr_matrix");

            /* raw/var with _index */
            char feat_path[256];
            snprintf(feat_path, sizeof(feat_path), "assays/%s/features", assay);
            if (sc_has_dataset(src, feat_path)) {
                hid_t fds = H5Dopen2(src, feat_path, H5P_DEFAULT);
                hid_t fsp = H5Dget_space(fds);
                hsize_t nf;
                H5Sget_simple_extent_dims(fsp, &nf, NULL);
                H5Sclose(fsp);
                hid_t stype = sc_create_vlen_str_type();
                char **feats = (char **)calloc(nf, sizeof(char *));
                H5Dread(fds, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, feats);
                H5Dclose(fds);

                char rv_dir[2048];
                snprintf(rv_dir, sizeof(rv_dir), "%s/raw/var", out);
                zarr_create_group(rv_dir);
                char rv_idx[2048];
                snprintf(rv_idx, sizeof(rv_idx), "%s/raw/var/_index", out);
                zarr_write_strings(rv_idx, (const char **)feats, (int64_t)nf, gzip);
                sc_json_write_zattrs_str(rv_dir,
                    "{\"encoding-type\": \"dataframe\", \"encoding-version\": \"0.2.0\", "
                    "\"_index\": \"_index\", \"column-order\": []}");

                for (hsize_t k = 0; k < nf; k++) free(feats[k]);
                free(feats);
                H5Tclose(stype);
            }

            H5Gclose(grp);
        }
    }

    /* ── 2. obs ────────────────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [2/6] Transferring obs...\n");
    {
        /* Read cell names */
        char **cells = NULL;
        if (sc_has_dataset(src, "cell.names")) {
            hid_t ds = H5Dopen2(src, "cell.names", H5P_DEFAULT);
            hid_t sp = H5Dget_space(ds);
            hsize_t nc;
            H5Sget_simple_extent_dims(sp, &nc, NULL);
            H5Sclose(sp);
            hid_t stype = sc_create_vlen_str_type();
            cells = (char **)calloc(nc, sizeof(char *));
            H5Dread(ds, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cells);
            H5Dclose(ds);
            H5Tclose(stype);
            if (n_cells == 0) n_cells = (int64_t)nc;
        }

        char obs_dir[2048];
        snprintf(obs_dir, sizeof(obs_dir), "%s/obs", out);
        if (sc_has_group(src, "meta.data")) {
            hid_t md = H5Gopen2(src, "meta.data", H5P_DEFAULT);
            h5seurat_md_to_zarr_df(md, obs_dir,
                                    (const char **)cells, n_cells, gzip);
            H5Gclose(md);
        } else {
            zarr_create_group(obs_dir);
            char idx_dir[2048];
            snprintf(idx_dir, sizeof(idx_dir), "%s/_index", obs_dir);
            if (cells)
                zarr_write_strings(idx_dir, (const char **)cells, n_cells, gzip);
            sc_json_write_zattrs_str(obs_dir,
                "{\"encoding-type\": \"dataframe\", \"encoding-version\": \"0.2.0\", "
                "\"_index\": \"_index\", \"column-order\": []}");
        }

        if (cells) {
            for (int64_t i = 0; i < n_cells; i++) free(cells[i]);
            free(cells);
        }
    }

    /* ── 3. var ────────────────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [3/6] Transferring var...\n");
    {
        char feat_path[256];
        snprintf(feat_path, sizeof(feat_path), "assays/%s/features", assay);
        char var_dir[2048];
        snprintf(var_dir, sizeof(var_dir), "%s/var", out);
        zarr_create_group(var_dir);

        if (sc_has_dataset(src, feat_path)) {
            hid_t ds = H5Dopen2(src, feat_path, H5P_DEFAULT);
            hid_t sp = H5Dget_space(ds);
            hsize_t nf;
            H5Sget_simple_extent_dims(sp, &nf, NULL);
            H5Sclose(sp);
            hid_t stype = sc_create_vlen_str_type();
            char **feats = (char **)calloc(nf, sizeof(char *));
            H5Dread(ds, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, feats);
            H5Dclose(ds);

            char idx_dir[2048];
            snprintf(idx_dir, sizeof(idx_dir), "%s/var/_index", out);
            zarr_write_strings(idx_dir, (const char **)feats, (int64_t)nf, gzip);
            if (n_genes == 0) n_genes = (int64_t)nf;

            for (hsize_t k = 0; k < nf; k++) free(feats[k]);
            free(feats);
            H5Tclose(stype);
        }
        sc_json_write_zattrs_str(var_dir,
            "{\"encoding-type\": \"dataframe\", \"encoding-version\": \"0.2.0\", "
            "\"_index\": \"_index\", \"column-order\": []}");
    }

    /* ── 4. obsm (reductions) ──────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [4/6] Transferring obsm...\n");
    {
        char obsm_dir[2048];
        snprintf(obsm_dir, sizeof(obsm_dir), "%s/obsm", out);
        zarr_create_group(obsm_dir);

        if (sc_has_group(src, "reductions")) {
            hid_t red = H5Gopen2(src, "reductions", H5P_DEFAULT);
            hsize_t n;
            H5Gget_num_objs(red, &n);
            for (hsize_t i = 0; i < n; i++) {
                char rname[256];
                H5Gget_objname_by_idx(red, i, rname, sizeof(rname));
                char rpath[512];
                snprintf(rpath, sizeof(rpath), "reductions/%s/cell.embeddings", rname);
                if (!sc_has_dataset(src, rpath)) continue;

                hid_t ds = H5Dopen2(src, rpath, H5P_DEFAULT);
                hid_t sp = H5Dget_space(ds);
                hsize_t dims[2];
                H5Sget_simple_extent_dims(sp, dims, NULL);
                H5Sclose(sp);

                /* h5seurat: [n_comp, n_cells] → zarr: [n_cells, n_comp] */
                /* h5seurat stores as HDF5 [n_comp, n_cells] */
                int64_t n_comp = (int64_t)dims[0];
                int64_t n_obs = (int64_t)dims[1];
                double *flat = (double *)malloc((size_t)(n_comp * n_obs) * sizeof(double));
                if (!flat) { H5Dclose(ds); continue; }
                H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat);
                H5Dclose(ds);

                /* Transpose [n_comp, n_obs] → [n_obs, n_comp] C-order */
                double *transposed = (double *)malloc((size_t)(n_comp * n_obs) * sizeof(double));
                if (!transposed) { free(flat); continue; }
                for (int64_t r = 0; r < n_comp; r++)
                    for (int64_t c = 0; c < n_obs; c++)
                        transposed[c * n_comp + r] = flat[r * n_obs + c];

                char emb_name[256];
                snprintf(emb_name, sizeof(emb_name), "X_%s", rname);
                char emb_dir[2048];
                snprintf(emb_dir, sizeof(emb_dir), "%s/%s", obsm_dir, emb_name);
                zarr_write_numeric_2d(emb_dir, transposed, n_obs, n_comp, "<f8", gzip);

                /* .zattrs for obsm arrays */
                sc_json_write_zattrs_str(emb_dir,
                    "{\"encoding-type\": \"array\", \"encoding-version\": \"0.2.0\"}");

                free(flat);
                free(transposed);
            }
            H5Gclose(red);
        }
    }

    /* ── 5. obsp (graphs) ──────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [5/6] Transferring obsp...\n");
    {
        char obsp_dir[2048];
        snprintf(obsp_dir, sizeof(obsp_dir), "%s/obsp", out);
        zarr_create_group(obsp_dir);

        if (sc_has_group(src, "graphs")) {
            hid_t g = H5Gopen2(src, "graphs", H5P_DEFAULT);
            hsize_t n;
            H5Gget_num_objs(g, &n);
            for (hsize_t i = 0; i < n; i++) {
                char gname[256];
                H5Gget_objname_by_idx(g, i, gname, sizeof(gname));

                char gpath[512];
                snprintf(gpath, sizeof(gpath), "graphs/%s", gname);
                if (!sc_has_group(src, gpath)) continue;

                hid_t gg = H5Gopen2(src, gpath, H5P_DEFAULT);
                if (!sc_has_dataset(gg, "data")) { H5Gclose(gg); continue; }

                /* Strip assay prefix (RNA_nn → nn) */
                char *short_name = gname;
                size_t alen = strlen(assay);
                if (strncmp(gname, assay, alen) == 0 && gname[alen] == '_')
                    short_name = gname + alen + 1;

                char g_dir[2048];
                snprintf(g_dir, sizeof(g_dir), "%s/%s", obsp_dir, short_name);
                h5seurat_sparse_to_zarr(gg, g_dir, gzip, "csc_matrix");
                H5Gclose(gg);
            }
            H5Gclose(g);
        }
    }

    /* ── 6. Empty groups ───────────────────────────────────────────────────── */
    if (opts->verbose) SC_MSG("  [6/6] Writing empty groups...\n");
    {
        const char *empty[] = {"layers", "uns", "varm", "varp"};
        for (int i = 0; i < 4; i++) {
            char dir[2048];
            snprintf(dir, sizeof(dir), "%s/%s", out, empty[i]);
            zarr_create_group(dir);
        }
    }

    H5Fclose(src);
    if (opts->verbose) SC_MSG("[scConvert] Done.\n");
    return SC_OK;
}

/* ══════════════════════════════════════════════════════════════════════════════
 *  Composite conversions via temp h5seurat
 * ══════════════════════════════════════════════════════════════════════════════ */

static int sc_convert_via_temp_h5seurat_zarr(
    const sc_opts_t *opts,
    int (*to_h5seurat)(const sc_opts_t *),
    int (*from_h5seurat)(const sc_opts_t *)
) {
    char tmp[1024];
    snprintf(tmp, sizeof(tmp), "%s.tmp.h5seurat", opts->output_path);
    sc_opts_t step1 = *opts;
    step1.output_path = tmp;
    int rc = to_h5seurat(&step1);
    if (rc != SC_OK) { unlink(tmp); return rc; }
    sc_opts_t step2 = *opts;
    step2.input_path = tmp;
    rc = from_h5seurat(&step2);
    unlink(tmp);
    return rc;
}

int sc_zarr_to_h5ad(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (zarr → h5ad via h5seurat)\n",
                opts->input_path, opts->output_path);
    return sc_convert_via_temp_h5seurat_zarr(opts,
                                              sc_zarr_to_h5seurat,
                                              sc_h5seurat_to_h5ad);
}

int sc_h5ad_to_zarr(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5ad → zarr via h5seurat)\n",
                opts->input_path, opts->output_path);
    return sc_convert_via_temp_h5seurat_zarr(opts,
                                              sc_h5ad_to_h5seurat,
                                              sc_h5seurat_to_zarr);
}

int sc_zarr_to_h5mu(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (zarr → h5mu via h5seurat)\n",
                opts->input_path, opts->output_path);
    return sc_convert_via_temp_h5seurat_zarr(opts,
                                              sc_zarr_to_h5seurat,
                                              sc_h5seurat_to_h5mu);
}

int sc_h5mu_to_zarr(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5mu → zarr via h5seurat)\n",
                opts->input_path, opts->output_path);
    return sc_convert_via_temp_h5seurat_zarr(opts,
                                              sc_h5mu_to_h5seurat,
                                              sc_h5seurat_to_zarr);
}

int sc_zarr_to_loom(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (zarr → loom via h5seurat)\n",
                opts->input_path, opts->output_path);
    return sc_convert_via_temp_h5seurat_zarr(opts,
                                              sc_zarr_to_h5seurat,
                                              sc_h5seurat_to_loom);
}

int sc_loom_to_zarr(const sc_opts_t *opts) {
    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (loom → zarr via h5seurat)\n",
                opts->input_path, opts->output_path);
    return sc_convert_via_temp_h5seurat_zarr(opts,
                                              sc_loom_to_h5seurat,
                                              sc_h5seurat_to_zarr);
}
