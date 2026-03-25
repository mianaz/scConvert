/*
 * sc_h5ad.c — Top-level h5ad ↔ h5seurat conversion orchestration
 *
 * Each conversion function opens source (read-only) and destination (create),
 * then orchestrates streaming transfer of all components.
 */

#include "sc_convert.h"
#include <ctype.h>

/* ── h5ad → h5seurat ────────────────────────────────────────────────────────── */

int sc_h5ad_to_h5seurat(const sc_opts_t *opts) {
    const char *assay = opts->assay_name ? opts->assay_name : "RNA";
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5ad → h5seurat)\n",
                opts->input_path, opts->output_path);

    /* Open source h5ad (read-only) */
    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        SC_MSG("Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    /* Create destination h5seurat */
    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (dst < 0) {
        SC_MSG("Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Set required h5seurat root attributes */
    sc_set_str_attr(dst, "version", "5.3.0");
    sc_set_str_attr(dst, "project", "scConvert");
    sc_set_str_attr(dst, "active.assay", assay);

    /* 1. Transfer X (expression matrix) */
    if (opts->verbose) SC_MSG("  [1/6] Transferring X...\n");
    int has_raw = sc_has_group(src, "raw");
    if (sc_has_group(src, "X")) {
        /* X is a sparse group (CSR) */
        hid_t assays = sc_create_or_open_group(dst, "assays");
        hid_t assay_grp = sc_create_or_open_group(assays, assay);

        /* Set required assay attributes */
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
        /* V5 s4class attribute — required for LoadH5Seurat to use the
         * V5 code path that looks for layers/ subdirectory */
        sc_set_str_attr(assay_grp, "s4class", "SeuratObject::Assay5");

        hid_t src_x = H5Gopen2(src, "X", H5P_DEFAULT);
        sc_csr_info_t info = {0};
        sc_read_csr_info(src_x, &info);
        int64_t dims[2] = {(int64_t)info.n_cols, (int64_t)info.n_rows};

        /* Seurat v5 expects layers under assays/{assay}/layers/ */
        hid_t layers_grp = sc_create_or_open_group(assay_grp, "layers");

        if (has_raw) {
            /* raw/X exists: X -> data (normalized), raw/X -> counts (handled below) */
            hid_t dst_data = sc_create_or_open_group(layers_grp, "data");
            rc = sc_stream_csr_copy(src_x, dst_data, gzip);
            sc_set_int_array_attr(dst_data, "dims", dims, 2);
            H5Gclose(dst_data);
        } else {
            /* No raw/X: X is the only matrix — write to both counts and data
             * so Seurat has proper layer access for both GetAssayData(layer="counts")
             * and GetAssayData(layer="data") */
            hid_t dst_counts = sc_create_or_open_group(layers_grp, "counts");
            rc = sc_stream_csr_copy(src_x, dst_counts, gzip);
            sc_set_int_array_attr(dst_counts, "dims", dims, 2);
            H5Gclose(dst_counts);
            if (rc == SC_OK) {
                /* Re-open source and copy to data as well */
                hid_t dst_data = sc_create_or_open_group(layers_grp, "data");
                rc = sc_stream_csr_copy(src_x, dst_data, gzip);
                sc_set_int_array_attr(dst_data, "dims", dims, 2);
                H5Gclose(dst_data);
            }
        }

        H5Gclose(layers_grp);

        H5Gclose(src_x);
        H5Gclose(assay_grp);
        H5Gclose(assays);

        if (rc != SC_OK) goto cleanup;
    } else if (sc_has_dataset(src, "X")) {
        /* X is a dense dataset — store in layers/data (or layers/counts
         * if no raw/X present, since it may be the only expression layer) */
        if (opts->verbose)
            SC_MSG("  Note: dense X matrix detected\n");
        hid_t assays = sc_create_or_open_group(dst, "assays");
        hid_t assay_grp = sc_create_or_open_group(assays, assay);

        /* Set required assay attributes (same as sparse path) */
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
        hid_t src_x = H5Dopen2(src, "X", H5P_DEFAULT);
        const char *layer_name = has_raw ? "data" : "counts";
        rc = sc_copy_dataset_chunked(src_x, layers_grp, layer_name, gzip);
        H5Dclose(src_x);
        H5Gclose(layers_grp);
        H5Gclose(assay_grp);
        H5Gclose(assays);
        if (rc != SC_OK) goto cleanup;
    }

    /* 2. Transfer raw/X (counts) if present */
    if (sc_has_group(src, "raw")) {
        if (opts->verbose) SC_MSG("  [1b] Transferring raw/X...\n");
        hid_t raw = H5Gopen2(src, "raw", H5P_DEFAULT);
        if (sc_has_group(raw, "X")) {
            /* Sparse raw/X */
            hid_t assays = sc_create_or_open_group(dst, "assays");
            hid_t assay_grp = sc_create_or_open_group(assays, assay);
            hid_t layers_grp = sc_create_or_open_group(assay_grp, "layers");

            hid_t src_raw_x = H5Gopen2(raw, "X", H5P_DEFAULT);
            hid_t dst_counts = sc_create_or_open_group(layers_grp, "counts");
            rc = sc_stream_csr_copy(src_raw_x, dst_counts, gzip);

            sc_csr_info_t info = {0};
            sc_read_csr_info(src_raw_x, &info);
            int64_t dims[2] = {(int64_t)info.n_cols, (int64_t)info.n_rows};
            sc_set_int_array_attr(dst_counts, "dims", dims, 2);

            H5Gclose(dst_counts);
            H5Gclose(src_raw_x);
            H5Gclose(layers_grp);
            H5Gclose(assay_grp);
            H5Gclose(assays);
        } else if (sc_has_dataset(raw, "X")) {
            /* Dense raw/X */
            hid_t assays = sc_create_or_open_group(dst, "assays");
            hid_t assay_grp = sc_create_or_open_group(assays, assay);
            hid_t layers_grp = sc_create_or_open_group(assay_grp, "layers");
            hid_t src_raw_x = H5Dopen2(raw, "X", H5P_DEFAULT);
            rc = sc_copy_dataset_chunked(src_raw_x, layers_grp, "counts", gzip);
            H5Dclose(src_raw_x);
            H5Gclose(layers_grp);
            H5Gclose(assay_grp);
            H5Gclose(assays);
        }
        H5Gclose(raw);
        if (rc != SC_OK) goto cleanup;
    }

    /* 2b. Copy cell names from obs/_index to top-level cell.names */
    if (sc_has_group(src, "obs")) {
        hid_t src_obs = H5Gopen2(src, "obs", H5P_DEFAULT);
        if (sc_has_dataset(src_obs, "_index")) {
            hid_t src_idx = H5Dopen2(src_obs, "_index", H5P_DEFAULT);
            sc_copy_dataset_chunked(src_idx, dst, "cell.names", gzip);
            H5Dclose(src_idx);
        }
        H5Gclose(src_obs);
    }

    /* 3. Transfer obs → meta.data */
    if (opts->verbose) SC_MSG("  [2/6] Transferring obs...\n");
    rc = sc_stream_obs_h5ad_to_h5seurat(src, dst, assay);
    if (rc != SC_OK) goto cleanup;

    /* 4. Transfer var → assays/<assay>/features */
    if (opts->verbose) SC_MSG("  [3/6] Transferring var...\n");
    rc = sc_stream_var_h5ad_to_h5seurat(src, dst, assay);
    if (rc != SC_OK) goto cleanup;

    /* 5. Transfer obsm → reductions */
    if (opts->verbose) SC_MSG("  [4/6] Transferring obsm...\n");
    rc = sc_stream_obsm(src, dst, SC_H5AD_TO_H5SEURAT);
    if (rc != SC_OK) goto cleanup;

    /* 6. Transfer obsp → graphs */
    if (opts->verbose) SC_MSG("  [5/6] Transferring obsp...\n");
    rc = sc_stream_obsp(src, dst, SC_H5AD_TO_H5SEURAT, assay);
    if (rc != SC_OK) goto cleanup;

    /* 7. Transfer uns → misc */
    if (opts->verbose) SC_MSG("  [6/6] Transferring uns...\n");
    rc = sc_stream_uns(src, dst, SC_H5AD_TO_H5SEURAT);
    if (rc != SC_OK) goto cleanup;

    /* 8. Transfer layers */
    rc = sc_stream_layers(src, dst, SC_H5AD_TO_H5SEURAT, gzip);
    if (rc != SC_OK) goto cleanup;

    /* 9. Create required empty groups for h5seurat */
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

    /* 10. Populate active.ident with default single-level identity */
    {
        hid_t ai = sc_create_or_open_group(dst, "active.ident");

        /* levels: one string "scConvert" */
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

        /* values: n_cells integers, all = 1 (1-based index) */
        if (!sc_has_dataset(ai, "values") && sc_has_dataset(dst, "cell.names")) {
            hid_t cn = H5Dopen2(dst, "cell.names", H5P_DEFAULT);
            hid_t cn_sp = H5Dget_space(cn);
            hsize_t n_cells;
            H5Sget_simple_extent_dims(cn_sp, &n_cells, NULL);
            H5Sclose(cn_sp);
            H5Dclose(cn);

            int *vals = (int *)malloc(n_cells * sizeof(int));
            if (!vals) { H5Gclose(ai); goto cleanup; }
            for (hsize_t c = 0; c < n_cells; c++)
                vals[c] = 1;

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

/* ── h5seurat → h5ad ────────────────────────────────────────────────────────── */

int sc_h5seurat_to_h5ad(const sc_opts_t *opts) {
    const char *assay = opts->assay_name ? opts->assay_name : "RNA";
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    if (opts->verbose)
        SC_MSG("[scConvert] Converting %s → %s (h5seurat → h5ad)\n",
                opts->input_path, opts->output_path);

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        SC_MSG("Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (dst < 0) {
        SC_MSG("Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Set root anndata encoding */
    sc_set_str_attr(dst, "encoding-type", "anndata");
    sc_set_str_attr(dst, "encoding-version", "0.1.0");

    /* 1. Transfer X (expression matrix) */
    if (opts->verbose) SC_MSG("  [1/6] Transferring X...\n");
    {
        /* Seurat v5 stores layers under assays/{assay}/layers/{name},
         * Seurat v4 stores directly under assays/{assay}/{name}.
         * Priority: layers/data > data > layers/counts > counts
         * When only counts exists (no data layer), counts becomes X. */
        char data_path[512];
        int x_from_counts = 0;
        snprintf(data_path, sizeof(data_path), "assays/%s/layers/data", assay);
        if (!sc_has_group(src, data_path) && !sc_has_dataset(src, data_path)) {
            snprintf(data_path, sizeof(data_path), "assays/%s/data", assay);
            if (!sc_has_group(src, data_path) && !sc_has_dataset(src, data_path)) {
                /* No data layer at all — use counts as X */
                snprintf(data_path, sizeof(data_path), "assays/%s/layers/counts", assay);
                if (!sc_has_group(src, data_path) && !sc_has_dataset(src, data_path))
                    snprintf(data_path, sizeof(data_path), "assays/%s/counts", assay);
                x_from_counts = 1;
            }
        }

        if (sc_has_group(src, data_path)) {
            hid_t src_data = H5Gopen2(src, data_path, H5P_DEFAULT);
            hid_t dst_x = H5Gcreate2(dst, "X", H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);

            /* h5seurat stores sparse as CSC (genes×cells)
             * h5ad needs CSR (cells×genes)
             * CSC of A == CSR of A^T — just stream-copy and swap dims */
            rc = sc_stream_csr_transpose(src_data, dst_x, gzip);

            /* Read h5seurat dims = (n_genes, n_cells) and write h5ad
             * shape = (n_cells, n_genes) — the swapped version */
            sc_csr_info_t info = {0};
            sc_read_csr_info(src_data, &info);
            /* info from h5seurat: n_rows=n_genes, n_cols=n_cells */
            int64_t shape[2] = {(int64_t)info.n_cols, (int64_t)info.n_rows};
            sc_set_str_attr(dst_x, "encoding-type", "csr_matrix");
            sc_set_str_attr(dst_x, "encoding-version", "0.1.0");
            sc_set_int_array_attr(dst_x, "shape", shape, 2);

            H5Gclose(dst_x);
            H5Gclose(src_data);
            if (rc != SC_OK) goto cleanup;
        } else if (sc_has_dataset(src, data_path)) {
            /* Dense data layer — copy as dense X dataset */
            if (opts->verbose)
                SC_MSG("  Note: dense data layer detected\n");
            hid_t src_data = H5Dopen2(src, data_path, H5P_DEFAULT);
            rc = sc_copy_dataset_chunked(src_data, dst, "X", gzip);
            H5Dclose(src_data);
            if (rc != SC_OK) goto cleanup;
        }

        /* Transfer counts to raw/X if present (skip if counts was already used as X) */
        char counts_path[512];
        int counts_is_sparse = 0;
        if (x_from_counts) {
            counts_path[0] = '\0';  /* already used as X */
        } else {
            snprintf(counts_path, sizeof(counts_path), "assays/%s/layers/counts", assay);
            if (sc_has_group(src, counts_path)) {
                counts_is_sparse = 1;
            } else if (sc_has_dataset(src, counts_path)) {
                counts_is_sparse = 0;
            } else {
                snprintf(counts_path, sizeof(counts_path), "assays/%s/counts", assay);
                if (sc_has_group(src, counts_path))
                    counts_is_sparse = 1;
                else if (sc_has_dataset(src, counts_path))
                    counts_is_sparse = 0;
                else
                    counts_path[0] = '\0';  /* no counts found */
            }
        }
        if (counts_path[0] != '\0' && counts_is_sparse) {
            if (opts->verbose)
                SC_MSG("  [1b] Transferring counts → raw/X...\n");
            hid_t raw = sc_create_or_open_group(dst, "raw");
            hid_t src_counts = H5Gopen2(src, counts_path, H5P_DEFAULT);
            hid_t dst_raw_x = H5Gcreate2(raw, "X", H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);

            rc = sc_stream_csr_transpose(src_counts, dst_raw_x, gzip);

            /* Set h5ad CSR encoding and shape */
            sc_csr_info_t cnt_info = {0};
            sc_read_csr_info(src_counts, &cnt_info);
            int64_t cnt_shape[2] = {(int64_t)cnt_info.n_cols,
                                    (int64_t)cnt_info.n_rows};
            sc_set_str_attr(dst_raw_x, "encoding-type", "csr_matrix");
            sc_set_str_attr(dst_raw_x, "encoding-version", "0.1.0");
            sc_set_int_array_attr(dst_raw_x, "shape", cnt_shape, 2);

            H5Gclose(dst_raw_x);
            H5Gclose(src_counts);

            /* Create raw/var with _index */
            hid_t raw_var = sc_create_or_open_group(raw, "var");
            sc_set_str_attr(raw_var, "encoding-type", "dataframe");
            sc_set_str_attr(raw_var, "encoding-version", "0.2.0");
            sc_set_str_attr(raw_var, "_index", "_index");

            /* Copy feature names to raw/var/_index */
            char feat_path[512];
            snprintf(feat_path, sizeof(feat_path), "assays/%s/features", assay);
            if (sc_has_dataset(src, feat_path)) {
                hid_t src_feat = H5Dopen2(src, feat_path, H5P_DEFAULT);
                sc_copy_dataset_chunked(src_feat, raw_var, "_index",
                                        SC_GZIP_LEVEL);
                hid_t rv_idx = H5Dopen2(raw_var, "_index", H5P_DEFAULT);
                sc_set_str_attr(rv_idx, "encoding-type", "string-array");
                sc_set_str_attr(rv_idx, "encoding-version", "0.2.0");
                H5Dclose(rv_idx);
                H5Dclose(src_feat);
            }

            /* Add empty column-order attribute (required by anndata) */
            {
                const char *empty = NULL;
                hsize_t zero = 0;
                hid_t space = H5Screate_simple(1, &zero, NULL);
                hid_t tid = sc_create_vlen_str_type();
                hid_t attr = H5Acreate2(raw_var, "column-order", tid, space,
                                         H5P_DEFAULT, H5P_DEFAULT);
                H5Awrite(attr, tid, &empty);
                H5Aclose(attr);
                H5Tclose(tid);
                H5Sclose(space);
            }

            H5Gclose(raw_var);
            H5Gclose(raw);
            if (rc != SC_OK) goto cleanup;
        } else if (counts_path[0] != '\0' && !counts_is_sparse) {
            /* Dense counts — copy as dense dataset to raw/X */
            if (opts->verbose)
                SC_MSG("  [1b] Transferring dense counts → raw/X...\n");
            hid_t raw = sc_create_or_open_group(dst, "raw");
            hid_t src_counts = H5Dopen2(src, counts_path, H5P_DEFAULT);
            rc = sc_copy_dataset_chunked(src_counts, raw, "X", gzip);
            H5Dclose(src_counts);
            H5Gclose(raw);
            if (rc != SC_OK) goto cleanup;
        }
    }

    /* 2. Transfer obs */
    if (opts->verbose) SC_MSG("  [2/6] Transferring obs...\n");
    rc = sc_stream_obs_h5seurat_to_h5ad(src, dst, assay);
    if (rc != SC_OK) goto cleanup;

    /* 3. Transfer var */
    if (opts->verbose) SC_MSG("  [3/6] Transferring var...\n");
    rc = sc_stream_var_h5seurat_to_h5ad(src, dst, assay);
    if (rc != SC_OK) goto cleanup;

    /* 4. Transfer obsm (reductions → obsm) */
    if (opts->verbose) SC_MSG("  [4/6] Transferring obsm...\n");
    rc = sc_stream_obsm(src, dst, SC_H5SEURAT_TO_H5AD);
    if (rc != SC_OK) goto cleanup;

    /* 5. Transfer obsp (graphs → obsp) */
    if (opts->verbose) SC_MSG("  [5/6] Transferring obsp...\n");
    rc = sc_stream_obsp(src, dst, SC_H5SEURAT_TO_H5AD, NULL);
    if (rc != SC_OK) goto cleanup;

    /* 6. Transfer uns */
    if (opts->verbose) SC_MSG("  [6/6] Transferring uns...\n");
    rc = sc_stream_uns(src, dst, SC_H5SEURAT_TO_H5AD);
    if (rc != SC_OK) goto cleanup;

    /* 7. Transfer layers */
    rc = sc_stream_layers(src, dst, SC_H5SEURAT_TO_H5AD, gzip);
    if (rc != SC_OK) goto cleanup;

    /* 8. Ensure empty groups exist */
    rc = sc_ensure_empty_groups(dst, SC_H5SEURAT_TO_H5AD);

    if (opts->verbose && rc == SC_OK)
        SC_MSG("[scConvert] Done.\n");

cleanup:
    H5Fclose(dst);
    H5Fclose(src);
    return rc;
}
