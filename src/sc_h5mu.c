/*
 * sc_h5mu.c — h5mu ↔ h5seurat and h5mu ↔ h5ad conversion
 *
 * Orchestrates streaming conversion of MuData (.h5mu) files, which store
 * multi-assay single-cell data as separate AnnData modalities under /mod/.
 *
 * Reuses existing per-group streaming primitives from sc_sparse.c,
 * sc_dataframe.c, and sc_groups.c — just changes the source/dest group handles.
 */

#include "sc_convert.h"
#include <ctype.h>

/* ── Internal helpers ──────────────────────────────────────────────────────── */

/*
 * Enumerate modality names from /mod/ group.
 * First tries the mod-order attribute; falls back to enumerating children.
 * Returns number of modalities found. names[] array is populated with
 * allocated strings that must be freed by the caller.
 */
static int sc_enumerate_modalities(hid_t h5mu, char **names, int max_n) {
    hid_t mod_grp = H5Gopen2(h5mu, "mod", H5P_DEFAULT);
    if (mod_grp < 0) return 0;

    int n_mod = 0;

    /* Try mod-order attribute first */
    if (H5Aexists(mod_grp, "mod-order") > 0) {
        char **order = NULL;
        hsize_t n_order = 0;
        if (sc_get_str_array_attr(mod_grp, "mod-order", &order, &n_order) == SC_OK
            && order && n_order > 0) {
            for (hsize_t i = 0; i < n_order && n_mod < max_n; i++) {
                names[n_mod] = strdup(order[i]);
                n_mod++;
            }
            sc_free_str_array(order, n_order);
            H5Gclose(mod_grp);
            return n_mod;
        }
        if (order) sc_free_str_array(order, n_order);
    }

    /* Fallback: enumerate children of /mod/ */
    hsize_t n_members;
    H5Gget_num_objs(mod_grp, &n_members);
    for (hsize_t i = 0; i < n_members && n_mod < max_n; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(mod_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(mod_grp, i);
        if (obj_type == H5G_GROUP) {
            names[n_mod] = strdup(name);
            n_mod++;
        }
    }

    H5Gclose(mod_grp);
    return n_mod;
}

/*
 * Stream a single modality from h5mu into an h5seurat assay group.
 * Handles X, layers/counts, var, and assay key attribute.
 */
static int sc_stream_modality_to_assay(
    hid_t mod_grp,          /* /mod/{modality}/ */
    hid_t h5seurat,         /* destination h5seurat file */
    const char *assay_name, /* e.g. "RNA" */
    int gzip,
    int verbose
) {
    int rc = SC_OK;
    char path_buf[512];

    hid_t assays = sc_create_or_open_group(h5seurat, "assays");
    hid_t assay_grp = sc_create_or_open_group(assays, assay_name);

    /* Set assay key attribute (lowercase + "_") */
    {
        char key_buf[SC_MAX_NAME_LEN];
        size_t len = strlen(assay_name);
        for (size_t k = 0; k < len && k < sizeof(key_buf) - 2; k++)
            key_buf[k] = (char)tolower((unsigned char)assay_name[k]);
        key_buf[len] = '_';
        key_buf[len + 1] = '\0';
        sc_set_str_attr(assay_grp, "key", key_buf);
    }
    /* V5 s4class attribute — required for LoadH5Seurat to use the
     * V5 code path that looks for layers/ subdirectory */
    sc_set_str_attr(assay_grp, "s4class", "SeuratObject::Assay5");

    /* X → assays/{assay}/data */
    if (sc_has_group(mod_grp, "X")) {
        if (verbose)
            fprintf(stderr, "    X → assays/%s/data\n", assay_name);
        hid_t src_x = H5Gopen2(mod_grp, "X", H5P_DEFAULT);
        hid_t dst_data = sc_create_or_open_group(assay_grp, "data");

        rc = sc_stream_csr_copy(src_x, dst_data, gzip);

        /* Set dims: h5seurat uses (n_genes, n_cells) = reversed shape */
        sc_csr_info_t info = {0};
        sc_read_csr_info(src_x, &info);
        int64_t dims[2] = {(int64_t)info.n_cols, (int64_t)info.n_rows};
        sc_set_int_array_attr(dst_data, "dims", dims, 2);

        H5Gclose(dst_data);
        H5Gclose(src_x);
        if (rc != SC_OK) goto done;
    } else if (sc_has_dataset(mod_grp, "X")) {
        if (verbose)
            fprintf(stderr, "    X (dense) → assays/%s/layers/data\n", assay_name);
        hid_t layers_grp = sc_create_or_open_group(assay_grp, "layers");
        hid_t src_x = H5Dopen2(mod_grp, "X", H5P_DEFAULT);
        rc = sc_copy_dataset_chunked(src_x, layers_grp, "data", gzip);
        H5Dclose(src_x);
        H5Gclose(layers_grp);
        if (rc != SC_OK) goto done;
    }

    /* layers/counts → assays/{assay}/counts */
    if (sc_has_group(mod_grp, "layers")) {
        hid_t layers = H5Gopen2(mod_grp, "layers", H5P_DEFAULT);
        if (sc_has_group(layers, "counts")) {
            if (verbose)
                fprintf(stderr, "    layers/counts → assays/%s/counts\n", assay_name);
            hid_t src_counts = H5Gopen2(layers, "counts", H5P_DEFAULT);
            hid_t dst_counts = sc_create_or_open_group(assay_grp, "counts");

            rc = sc_stream_csr_copy(src_counts, dst_counts, gzip);

            sc_csr_info_t cinfo = {0};
            sc_read_csr_info(src_counts, &cinfo);
            int64_t cdims[2] = {(int64_t)cinfo.n_cols, (int64_t)cinfo.n_rows};
            sc_set_int_array_attr(dst_counts, "dims", cdims, 2);

            H5Gclose(dst_counts);
            H5Gclose(src_counts);
            if (rc != SC_OK) { H5Gclose(layers); goto done; }
        }
        H5Gclose(layers);
    }

    /* var/_index → assays/{assay}/features */
    if (sc_has_group(mod_grp, "var")) {
        hid_t var_grp = H5Gopen2(mod_grp, "var", H5P_DEFAULT);
        if (sc_has_dataset(var_grp, "_index")) {
            if (verbose)
                fprintf(stderr, "    var/_index → assays/%s/features\n", assay_name);
            hid_t src_idx = H5Dopen2(var_grp, "_index", H5P_DEFAULT);
            sc_copy_dataset_chunked(src_idx, assay_grp, "features", gzip);
            H5Dclose(src_idx);
        }
        H5Gclose(var_grp);
    }

    /* obsm → per-modality reductions */
    if (sc_has_group(mod_grp, "obsm")) {
        hid_t src_obsm = H5Gopen2(mod_grp, "obsm", H5P_DEFAULT);
        hid_t dst_reducs = sc_create_or_open_group(h5seurat, "reductions");

        hsize_t n;
        H5Gget_num_objs(src_obsm, &n);
        for (hsize_t i = 0; i < n; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_obsm, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(src_obsm, i);

            if (obj_type == H5G_DATASET) {
                /* Strip X_ prefix for reduction name */
                const char *red_name = name;
                if (strncmp(name, "X_", 2) == 0)
                    red_name = name + 2;

                if (verbose)
                    fprintf(stderr, "    obsm/%s → reductions/%s\n", name, red_name);

                hid_t red_grp = sc_create_or_open_group(dst_reducs, red_name);
                hid_t src_dset = H5Dopen2(src_obsm, name, H5P_DEFAULT);
                sc_copy_dataset_chunked(src_dset, red_grp, "cell.embeddings", gzip);
                H5Dclose(src_dset);

                sc_set_str_attr(red_grp, "active.assay", assay_name);
                {
                    char key[SC_MAX_NAME_LEN];
                    if (strcmp(red_name, "pca") == 0)
                        strcpy(key, "PC_");
                    else if (strcmp(red_name, "tsne") == 0)
                        strcpy(key, "tSNE_");
                    else {
                        size_t k;
                        for (k = 0; red_name[k] && k < sizeof(key) - 2; k++)
                            key[k] = (char)toupper((unsigned char)red_name[k]);
                        key[k] = '_';
                        key[k + 1] = '\0';
                    }
                    sc_set_str_attr(red_grp, "key", key);
                }
                {
                    int global_val = 0;
                    hid_t s = H5Screate(H5S_SCALAR);
                    hid_t a = H5Acreate2(red_grp, "global", H5T_NATIVE_INT,
                                          s, H5P_DEFAULT, H5P_DEFAULT);
                    H5Awrite(a, H5T_NATIVE_INT, &global_val);
                    H5Aclose(a);
                    H5Sclose(s);
                }
                {
                    hid_t misc = H5Gcreate2(red_grp, "misc", H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
                    H5Gclose(misc);
                }
                H5Gclose(red_grp);
            }
        }
        H5Gclose(dst_reducs);
        H5Gclose(src_obsm);
    }

    /* obsp → graphs with assay-prefixed names */
    if (sc_has_group(mod_grp, "obsp")) {
        hid_t src_obsp = H5Gopen2(mod_grp, "obsp", H5P_DEFAULT);
        hid_t dst_graphs = sc_create_or_open_group(h5seurat, "graphs");

        hsize_t n;
        H5Gget_num_objs(src_obsp, &n);
        for (hsize_t i = 0; i < n; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_obsp, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(src_obsp, i);

            if (obj_type == H5G_GROUP) {
                /* Prefix with assay name: e.g. "connectivities" → "RNA_snn" */
                snprintf(path_buf, sizeof(path_buf), "%s_%s", assay_name, name);
                if (verbose)
                    fprintf(stderr, "    obsp/%s → graphs/%s\n", name, path_buf);

                hid_t src_child = H5Gopen2(src_obsp, name, H5P_DEFAULT);
                hid_t dst_child = sc_create_or_open_group(dst_graphs, path_buf);

                sc_copy_group_recursive(src_child, dst_child, gzip);

                /* Set assay.used attribute */
                sc_set_str_attr(dst_child, "assay.used", assay_name);

                H5Gclose(dst_child);
                H5Gclose(src_child);
            }
        }
        H5Gclose(dst_graphs);
        H5Gclose(src_obsp);
    }

done:
    H5Gclose(assay_grp);
    H5Gclose(assays);
    return rc;
}

/*
 * Stream a single h5seurat assay into an h5mu modality group.
 * Creates a complete AnnData structure under /mod/{modality}/.
 */
static int sc_stream_assay_to_modality(
    hid_t h5seurat,          /* source h5seurat file */
    hid_t mod_grp,           /* /mod/{modality}/ */
    const char *assay_name,  /* e.g. "RNA" */
    const char *modality,    /* e.g. "rna" */
    int gzip,
    int verbose
) {
    int rc = SC_OK;
    char path_buf[512];

    /* Set AnnData encoding on modality group */
    sc_set_str_attr(mod_grp, "encoding-type", "anndata");
    sc_set_str_attr(mod_grp, "encoding-version", "0.1.0");

    /* assays/{assay}/data → X (CSR, transposing dims) */
    snprintf(path_buf, sizeof(path_buf), "assays/%s/data", assay_name);
    if (sc_has_group(h5seurat, path_buf)) {
        if (verbose)
            fprintf(stderr, "    assays/%s/data → mod/%s/X\n", assay_name, modality);
        hid_t src_data = H5Gopen2(h5seurat, path_buf, H5P_DEFAULT);
        hid_t dst_x = H5Gcreate2(mod_grp, "X", H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);

        rc = sc_stream_csr_transpose(src_data, dst_x, gzip);

        /* h5seurat dims = (n_genes, n_cells), h5ad shape = (n_cells, n_genes) */
        sc_csr_info_t info = {0};
        sc_read_csr_info(src_data, &info);
        int64_t shape[2] = {(int64_t)info.n_cols, (int64_t)info.n_rows};
        sc_set_str_attr(dst_x, "encoding-type", "csr_matrix");
        sc_set_str_attr(dst_x, "encoding-version", "0.1.0");
        sc_set_int_array_attr(dst_x, "shape", shape, 2);

        H5Gclose(dst_x);
        H5Gclose(src_data);
        if (rc != SC_OK) return rc;
    }

    /* assays/{assay}/counts → layers/counts */
    snprintf(path_buf, sizeof(path_buf), "assays/%s/counts", assay_name);
    if (sc_has_group(h5seurat, path_buf)) {
        if (verbose)
            fprintf(stderr, "    assays/%s/counts → mod/%s/layers/counts\n",
                    assay_name, modality);
        hid_t layers = sc_create_or_open_group(mod_grp, "layers");
        sc_set_str_attr(layers, "encoding-type", "dict");
        sc_set_str_attr(layers, "encoding-version", "0.1.0");

        hid_t src_counts = H5Gopen2(h5seurat, path_buf, H5P_DEFAULT);
        hid_t dst_counts = H5Gcreate2(layers, "counts", H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);

        rc = sc_stream_csr_transpose(src_counts, dst_counts, gzip);

        sc_csr_info_t cinfo = {0};
        sc_read_csr_info(src_counts, &cinfo);
        int64_t cshape[2] = {(int64_t)cinfo.n_cols, (int64_t)cinfo.n_rows};
        sc_set_str_attr(dst_counts, "encoding-type", "csr_matrix");
        sc_set_str_attr(dst_counts, "encoding-version", "0.1.0");
        sc_set_int_array_attr(dst_counts, "shape", cshape, 2);

        H5Gclose(dst_counts);
        H5Gclose(src_counts);
        H5Gclose(layers);
        if (rc != SC_OK) return rc;
    }

    /* assays/{assay}/features → var/_index */
    {
        hid_t var = sc_create_or_open_group(mod_grp, "var");
        sc_set_str_attr(var, "encoding-type", "dataframe");
        sc_set_str_attr(var, "encoding-version", "0.2.0");
        sc_set_str_attr(var, "_index", "_index");

        snprintf(path_buf, sizeof(path_buf), "assays/%s/features", assay_name);
        if (sc_has_dataset(h5seurat, path_buf)) {
            if (verbose)
                fprintf(stderr, "    assays/%s/features → mod/%s/var/_index\n",
                        assay_name, modality);
            hid_t src_feat = H5Dopen2(h5seurat, path_buf, H5P_DEFAULT);
            sc_copy_dataset_chunked(src_feat, var, "_index", gzip);
            H5Dclose(src_feat);
        }
        H5Gclose(var);
    }

    /* Per-modality obs: just _index (cell names) for now */
    {
        hid_t obs = sc_create_or_open_group(mod_grp, "obs");
        sc_set_str_attr(obs, "encoding-type", "dataframe");
        sc_set_str_attr(obs, "encoding-version", "0.2.0");
        sc_set_str_attr(obs, "_index", "_index");

        if (sc_has_dataset(h5seurat, "cell.names")) {
            hid_t src_cn = H5Dopen2(h5seurat, "cell.names", H5P_DEFAULT);
            sc_copy_dataset_chunked(src_cn, obs, "_index", gzip);
            H5Dclose(src_cn);
        }
        H5Gclose(obs);
    }

    /* Filter reductions → obsm: only those with active.assay matching */
    if (sc_has_group(h5seurat, "reductions")) {
        hid_t src_reducs = H5Gopen2(h5seurat, "reductions", H5P_DEFAULT);
        hid_t dst_obsm = sc_create_or_open_group(mod_grp, "obsm");
        sc_set_str_attr(dst_obsm, "encoding-type", "dict");
        sc_set_str_attr(dst_obsm, "encoding-version", "0.1.0");

        hsize_t n;
        H5Gget_num_objs(src_reducs, &n);
        for (hsize_t i = 0; i < n; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_reducs, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(src_reducs, i);

            if (obj_type != H5G_GROUP) continue;

            hid_t red_grp = H5Gopen2(src_reducs, name, H5P_DEFAULT);

            /* Check active.assay attribute */
            char aa_buf[SC_MAX_NAME_LEN] = {0};
            sc_get_str_attr(red_grp, "active.assay", aa_buf, sizeof(aa_buf));

            /* Include if assay matches or attribute is missing/empty */
            if (aa_buf[0] == '\0' || strcmp(aa_buf, assay_name) == 0) {
                if (sc_has_dataset(red_grp, "cell.embeddings")) {
                    char dst_name[SC_MAX_NAME_LEN];
                    snprintf(dst_name, sizeof(dst_name), "X_%s", name);

                    if (verbose)
                        fprintf(stderr, "    reductions/%s → mod/%s/obsm/%s\n",
                                name, modality, dst_name);

                    hid_t src_emb = H5Dopen2(red_grp, "cell.embeddings",
                                              H5P_DEFAULT);
                    sc_copy_dataset_chunked(src_emb, dst_obsm, dst_name, gzip);

                    hid_t dst_dset = H5Dopen2(dst_obsm, dst_name, H5P_DEFAULT);
                    sc_set_str_attr(dst_dset, "encoding-type", "array");
                    sc_set_str_attr(dst_dset, "encoding-version", "0.2.0");
                    H5Dclose(dst_dset);
                    H5Dclose(src_emb);
                }
            }
            H5Gclose(red_grp);
        }
        H5Gclose(dst_obsm);
        H5Gclose(src_reducs);
    }

    /* Filter graphs → obsp: only those with assay.used matching */
    if (sc_has_group(h5seurat, "graphs")) {
        hid_t src_graphs = H5Gopen2(h5seurat, "graphs", H5P_DEFAULT);
        hid_t dst_obsp = sc_create_or_open_group(mod_grp, "obsp");
        sc_set_str_attr(dst_obsp, "encoding-type", "dict");
        sc_set_str_attr(dst_obsp, "encoding-version", "0.1.0");

        hsize_t n;
        H5Gget_num_objs(src_graphs, &n);
        for (hsize_t i = 0; i < n; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_graphs, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(src_graphs, i);

            if (obj_type != H5G_GROUP) continue;

            hid_t graph_grp = H5Gopen2(src_graphs, name, H5P_DEFAULT);

            /* Check assay.used attribute */
            char au_buf[SC_MAX_NAME_LEN] = {0};
            sc_get_str_attr(graph_grp, "assay.used", au_buf, sizeof(au_buf));

            if (au_buf[0] == '\0' || strcmp(au_buf, assay_name) == 0) {
                /* Strip assay prefix if present (e.g. "RNA_snn" → "snn") */
                const char *dst_name = name;
                size_t prefix_len = strlen(assay_name);
                if (strncmp(name, assay_name, prefix_len) == 0
                    && name[prefix_len] == '_') {
                    dst_name = name + prefix_len + 1;
                }

                if (verbose)
                    fprintf(stderr, "    graphs/%s → mod/%s/obsp/%s\n",
                            name, modality, dst_name);

                hid_t dst_child = sc_create_or_open_group(dst_obsp, dst_name);
                sc_copy_group_recursive(graph_grp, dst_child, gzip);

                /* Set CSR encoding */
                sc_set_str_attr(dst_child, "encoding-type", "csr_matrix");
                sc_set_str_attr(dst_child, "encoding-version", "0.1.0");

                H5Gclose(dst_child);
            }
            H5Gclose(graph_grp);
        }
        H5Gclose(dst_obsp);
        H5Gclose(src_graphs);
    }

    /* Ensure empty dict groups for completeness */
    {
        const char *empty_groups[] = {"obsm", "obsp", "varm", "varp", "layers", "uns"};
        for (int i = 0; i < 6; i++) {
            if (!sc_has_group(mod_grp, empty_groups[i])) {
                hid_t g = sc_create_or_open_group(mod_grp, empty_groups[i]);
                sc_set_str_attr(g, "encoding-type", "dict");
                sc_set_str_attr(g, "encoding-version", "0.1.0");
                H5Gclose(g);
            }
        }
    }

    return rc;
}

/* ── Public API ────────────────────────────────────────────────────────────── */

/*
 * sc_h5mu_to_h5seurat — Convert .h5mu to .h5seurat
 *
 * Reads each modality from /mod/{name}/ and writes corresponding assay
 * groups into the h5seurat file. Global /obs/ → meta.data.
 */
int sc_h5mu_to_h5seurat(const sc_opts_t *opts) {
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    if (opts->verbose)
        fprintf(stderr, "[scConvert] Converting %s → %s (h5mu → h5seurat)\n",
                opts->input_path, opts->output_path);

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        fprintf(stderr, "Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    /* Verify MuData encoding */
    {
        char enc[SC_MAX_NAME_LEN] = {0};
        sc_get_str_attr(src, "encoding-type", enc, sizeof(enc));
        if (enc[0] != '\0' && strcmp(enc, "MuData") != 0 &&
            strcmp(enc, "mudata") != 0) {
            fprintf(stderr, "Warning: expected MuData encoding, got '%s'\n", enc);
        }
    }

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dst < 0) {
        fprintf(stderr, "Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Enumerate modalities */
    char *mod_names[SC_MAX_MODALITIES];
    int n_mod = sc_enumerate_modalities(src, mod_names, SC_MAX_MODALITIES);
    if (n_mod == 0) {
        fprintf(stderr, "Error: no modalities found in /mod/\n");
        rc = SC_ERR_ARG;
        goto cleanup;
    }

    if (opts->verbose) {
        fprintf(stderr, "  Found %d modalities:", n_mod);
        for (int i = 0; i < n_mod; i++)
            fprintf(stderr, " %s", mod_names[i]);
        fprintf(stderr, "\n");
    }

    /* Determine active assay: prefer "rna", else first modality */
    const char *active_assay = NULL;
    for (int i = 0; i < n_mod; i++) {
        const char *assay = sc_modality_to_assay(mod_names[i]);
        if (strcmp(assay, "RNA") == 0) {
            active_assay = assay;
            break;
        }
    }
    if (!active_assay)
        active_assay = sc_modality_to_assay(mod_names[0]);

    /* Set h5seurat root attributes */
    sc_set_str_attr(dst, "version", "5.3.0");
    sc_set_str_attr(dst, "project", "scConvert");
    sc_set_str_attr(dst, "active.assay", active_assay);

    /* Process each modality */
    for (int i = 0; i < n_mod; i++) {
        const char *assay = sc_modality_to_assay(mod_names[i]);

        if (opts->verbose)
            fprintf(stderr, "  [%d/%d] Modality '%s' → assay '%s'\n",
                    i + 1, n_mod, mod_names[i], assay);

        char mod_path[512];
        snprintf(mod_path, sizeof(mod_path), "mod/%s", mod_names[i]);
        hid_t mod_grp = H5Gopen2(src, mod_path, H5P_DEFAULT);
        if (mod_grp < 0) {
            fprintf(stderr, "Warning: cannot open modality '%s'\n", mod_names[i]);
            continue;
        }

        rc = sc_stream_modality_to_assay(mod_grp, dst, assay, gzip,
                                          opts->verbose);
        H5Gclose(mod_grp);
        if (rc != SC_OK) goto cleanup;
    }

    /* Global obs → meta.data */
    if (opts->verbose) fprintf(stderr, "  Transferring global obs → meta.data\n");
    if (sc_has_group(src, "obs")) {
        hid_t src_obs = H5Gopen2(src, "obs", H5P_DEFAULT);
        hid_t dst_meta = sc_create_or_open_group(dst, "meta.data");
        sc_stream_df_group(src_obs, dst_meta, gzip);

        /* Copy cell names */
        if (sc_has_dataset(src_obs, "_index")) {
            hid_t src_idx = H5Dopen2(src_obs, "_index", H5P_DEFAULT);
            sc_copy_dataset_chunked(src_idx, dst, "cell.names", gzip);
            H5Dclose(src_idx);
        }

        H5Gclose(dst_meta);
        H5Gclose(src_obs);
    }

    /* Create active.ident */
    {
        hid_t ai = sc_create_or_open_group(dst, "active.ident");

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

        if (sc_has_dataset(dst, "cell.names")) {
            hid_t cn = H5Dopen2(dst, "cell.names", H5P_DEFAULT);
            hid_t cn_sp = H5Dget_space(cn);
            hsize_t n_cells;
            H5Sget_simple_extent_dims(cn_sp, &n_cells, NULL);
            H5Sclose(cn_sp);
            H5Dclose(cn);

            int *vals = (int *)calloc(n_cells, sizeof(int));
            if (vals) {
                for (hsize_t c = 0; c < n_cells; c++)
                    vals[c] = 1;

                hid_t vs = H5Screate_simple(1, &n_cells, NULL);
                hid_t vdcpl = H5Pcreate(H5P_DATASET_CREATE);
                hsize_t chunk = n_cells;
                H5Pset_chunk(vdcpl, 1, &chunk);
                hid_t vds = H5Dcreate2(ai, "values", H5T_NATIVE_INT, vs,
                                        H5P_DEFAULT, vdcpl, H5P_DEFAULT);
                H5Dwrite(vds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
                H5Dclose(vds);
                H5Pclose(vdcpl);
                H5Sclose(vs);
                free(vals);
            }
        }

        H5Gclose(ai);
    }

    /* Create required empty groups */
    {
        const char *required[] = {"tools", "commands", "images",
                                   "neighbors", "misc", "graphs"};
        for (int i = 0; i < 6; i++) {
            if (!sc_has_group(dst, required[i])) {
                hid_t g = H5Gcreate2(dst, required[i], H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
                H5Gclose(g);
            }
        }
    }

    if (opts->verbose && rc == SC_OK)
        fprintf(stderr, "[scConvert] Done.\n");

cleanup:
    H5Fclose(dst);
    H5Fclose(src);
    for (int i = 0; i < n_mod; i++) free(mod_names[i]);
    return rc;
}

/*
 * sc_h5seurat_to_h5mu — Convert .h5seurat to .h5mu
 *
 * Enumerates assays in h5seurat and writes each as a separate modality.
 */
int sc_h5seurat_to_h5mu(const sc_opts_t *opts) {
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    if (opts->verbose)
        fprintf(stderr, "[scConvert] Converting %s → %s (h5seurat → h5mu)\n",
                opts->input_path, opts->output_path);

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        fprintf(stderr, "Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dst < 0) {
        fprintf(stderr, "Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Set MuData root encoding */
    sc_set_str_attr(dst, "encoding-type", "MuData");
    sc_set_str_attr(dst, "encoding-version", "0.1.0");

    /* Enumerate assays from /assays/ */
    char *assay_names[SC_MAX_MODALITIES];
    int n_assays = 0;

    if (sc_has_group(src, "assays")) {
        hid_t assays = H5Gopen2(src, "assays", H5P_DEFAULT);
        hsize_t n_members;
        H5Gget_num_objs(assays, &n_members);

        for (hsize_t i = 0; i < n_members && n_assays < SC_MAX_MODALITIES; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(assays, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(assays, i);
            if (obj_type == H5G_GROUP) {
                assay_names[n_assays] = strdup(name);
                n_assays++;
            }
        }
        H5Gclose(assays);
    }

    if (n_assays == 0) {
        fprintf(stderr, "Error: no assays found in /assays/\n");
        rc = SC_ERR_ARG;
        goto cleanup;
    }

    if (opts->verbose) {
        fprintf(stderr, "  Found %d assays:", n_assays);
        for (int i = 0; i < n_assays; i++)
            fprintf(stderr, " %s", assay_names[i]);
        fprintf(stderr, "\n");
    }

    /* Create /mod/ group with mod-order attribute */
    {
        hid_t mod_root = H5Gcreate2(dst, "mod", H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
        const char *mod_order[SC_MAX_MODALITIES];
        for (int i = 0; i < n_assays; i++)
            mod_order[i] = sc_assay_to_modality(assay_names[i]);

        sc_set_str_array_attr(mod_root, "mod-order", mod_order,
                               (hsize_t)n_assays);
        H5Gclose(mod_root);
    }

    /* Process each assay */
    for (int i = 0; i < n_assays; i++) {
        const char *modality = sc_assay_to_modality(assay_names[i]);

        if (opts->verbose)
            fprintf(stderr, "  [%d/%d] Assay '%s' → modality '%s'\n",
                    i + 1, n_assays, assay_names[i], modality);

        char mod_path[512];
        snprintf(mod_path, sizeof(mod_path), "mod/%s", modality);
        hid_t mod_grp = H5Gcreate2(dst, mod_path, H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);

        rc = sc_stream_assay_to_modality(src, mod_grp, assay_names[i],
                                          modality, gzip, opts->verbose);
        H5Gclose(mod_grp);
        if (rc != SC_OK) goto cleanup;
    }

    /* Global obs from meta.data */
    if (opts->verbose) fprintf(stderr, "  Transferring meta.data → global obs\n");
    {
        hid_t obs = sc_create_or_open_group(dst, "obs");
        sc_set_str_attr(obs, "encoding-type", "dataframe");
        sc_set_str_attr(obs, "encoding-version", "0.2.0");
        sc_set_str_attr(obs, "_index", "_index");

        if (sc_has_group(src, "meta.data")) {
            hid_t src_meta = H5Gopen2(src, "meta.data", H5P_DEFAULT);
            sc_stream_df_group(src_meta, obs, gzip);
            H5Gclose(src_meta);
        }

        /* Write cell names as _index */
        if (sc_has_dataset(src, "cell.names")) {
            hid_t src_cn = H5Dopen2(src, "cell.names", H5P_DEFAULT);
            sc_copy_dataset_chunked(src_cn, obs, "_index", gzip);
            H5Dclose(src_cn);
        }
        H5Gclose(obs);
    }

    /* Create empty global obsm group */
    {
        hid_t obsm = sc_create_or_open_group(dst, "obsm");
        sc_set_str_attr(obsm, "encoding-type", "dict");
        sc_set_str_attr(obsm, "encoding-version", "0.1.0");
        H5Gclose(obsm);
    }

    if (opts->verbose && rc == SC_OK)
        fprintf(stderr, "[scConvert] Done.\n");

cleanup:
    H5Fclose(dst);
    H5Fclose(src);
    for (int i = 0; i < n_assays; i++) free(assay_names[i]);
    return rc;
}

/*
 * sc_h5mu_to_h5ad — Extract a single modality from h5mu as h5ad.
 *
 * If --assay is specified, extracts that modality. Otherwise extracts "rna"
 * (or the first modality).
 */
int sc_h5mu_to_h5ad(const sc_opts_t *opts) {
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    /* Determine which modality to extract */
    const char *target_mod = NULL;
    if (opts->assay_name && strcmp(opts->assay_name, "RNA") != 0) {
        /* User specified an assay — map to modality name */
        target_mod = sc_assay_to_modality(opts->assay_name);
    }

    if (opts->verbose)
        fprintf(stderr, "[scConvert] Converting %s → %s (h5mu → h5ad, modality: %s)\n",
                opts->input_path, opts->output_path,
                target_mod ? target_mod : "auto");

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        fprintf(stderr, "Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    /* Find the target modality */
    char *mod_names[SC_MAX_MODALITIES];
    int n_mod = sc_enumerate_modalities(src, mod_names, SC_MAX_MODALITIES);
    if (n_mod == 0) {
        fprintf(stderr, "Error: no modalities found\n");
        H5Fclose(src);
        return SC_ERR_ARG;
    }

    /* Resolve target modality */
    int target_idx = -1;
    if (target_mod) {
        for (int i = 0; i < n_mod; i++) {
            if (strcasecmp(mod_names[i], target_mod) == 0) {
                target_idx = i;
                break;
            }
        }
    } else {
        /* Prefer "rna", else first */
        for (int i = 0; i < n_mod; i++) {
            if (strcasecmp(mod_names[i], "rna") == 0) {
                target_idx = i;
                break;
            }
        }
        if (target_idx < 0) target_idx = 0;
    }

    if (target_idx < 0) {
        fprintf(stderr, "Error: modality '%s' not found in h5mu file\n", target_mod);
        for (int i = 0; i < n_mod; i++) free(mod_names[i]);
        H5Fclose(src);
        return SC_ERR_ARG;
    }

    if (opts->verbose)
        fprintf(stderr, "  Extracting modality '%s'\n", mod_names[target_idx]);

    /* Open modality group and copy as complete h5ad */
    char mod_path[512];
    snprintf(mod_path, sizeof(mod_path), "mod/%s", mod_names[target_idx]);
    hid_t mod_grp = H5Gopen2(src, mod_path, H5P_DEFAULT);

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dst < 0) {
        fprintf(stderr, "Error: cannot create output file: %s\n", opts->output_path);
        H5Gclose(mod_grp);
        H5Fclose(src);
        for (int i = 0; i < n_mod; i++) free(mod_names[i]);
        return SC_ERR_IO;
    }

    /* The modality group IS a complete AnnData — recursive copy */
    rc = sc_copy_group_recursive(mod_grp, dst, gzip);

    /* Merge global obs columns if present */
    if (rc == SC_OK && sc_has_group(src, "obs")) {
        hid_t src_obs = H5Gopen2(src, "obs", H5P_DEFAULT);
        hid_t dst_obs;
        if (sc_has_group(dst, "obs"))
            dst_obs = H5Gopen2(dst, "obs", H5P_DEFAULT);
        else
            dst_obs = sc_create_or_open_group(dst, "obs");

        /* Copy columns from global obs that don't exist in modality obs */
        hsize_t n_members;
        H5Gget_num_objs(src_obs, &n_members);
        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_obs, i, name, sizeof(name));
            /* Skip _index (already in modality obs) */
            if (strcmp(name, "_index") == 0) continue;
            /* Skip if already exists in destination obs */
            if (sc_has_dataset(dst_obs, name) || sc_has_group(dst_obs, name))
                continue;

            int obj_type = H5Gget_objtype_by_idx(src_obs, i);
            if (obj_type == H5G_DATASET) {
                hid_t src_dset = H5Dopen2(src_obs, name, H5P_DEFAULT);
                sc_copy_dataset_chunked(src_dset, dst_obs, name, gzip);
                H5Dclose(src_dset);
            } else if (obj_type == H5G_GROUP) {
                hid_t src_child = H5Gopen2(src_obs, name, H5P_DEFAULT);
                hid_t dst_child = H5Gcreate2(dst_obs, name, H5P_DEFAULT,
                                              H5P_DEFAULT, H5P_DEFAULT);
                sc_copy_group_recursive(src_child, dst_child, gzip);
                H5Gclose(dst_child);
                H5Gclose(src_child);
            }
        }
        H5Gclose(dst_obs);
        H5Gclose(src_obs);
    }

    /* Set root anndata encoding if not already present */
    sc_set_str_attr(dst, "encoding-type", "anndata");
    sc_set_str_attr(dst, "encoding-version", "0.1.0");

    /* Ensure empty dict groups */
    rc = sc_ensure_empty_groups(dst, SC_H5SEURAT_TO_H5AD);

    if (opts->verbose && rc == SC_OK)
        fprintf(stderr, "[scConvert] Done.\n");

    H5Gclose(mod_grp);
    H5Fclose(dst);
    H5Fclose(src);
    for (int i = 0; i < n_mod; i++) free(mod_names[i]);
    return rc;
}

/*
 * sc_h5ad_to_h5mu — Wrap a single h5ad as a one-modality h5mu.
 *
 * Creates /mod/{modality}/ and copies the entire h5ad structure into it.
 */
int sc_h5ad_to_h5mu(const sc_opts_t *opts) {
    int gzip = opts->gzip_level > 0 ? opts->gzip_level : SC_GZIP_LEVEL;
    int rc = SC_OK;

    const char *assay = opts->assay_name ? opts->assay_name : "RNA";
    const char *modality = sc_assay_to_modality(assay);

    if (opts->verbose)
        fprintf(stderr, "[scConvert] Converting %s → %s (h5ad → h5mu, modality: %s)\n",
                opts->input_path, opts->output_path, modality);

    hid_t src = H5Fopen(opts->input_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (src < 0) {
        fprintf(stderr, "Error: cannot open input file: %s\n", opts->input_path);
        return SC_ERR_IO;
    }

    hid_t dst = H5Fcreate(opts->output_path, H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dst < 0) {
        fprintf(stderr, "Error: cannot create output file: %s\n", opts->output_path);
        H5Fclose(src);
        return SC_ERR_IO;
    }

    /* Set MuData root encoding */
    sc_set_str_attr(dst, "encoding-type", "MuData");
    sc_set_str_attr(dst, "encoding-version", "0.1.0");

    /* Create /mod/ with mod-order */
    hid_t mod_root = H5Gcreate2(dst, "mod", H5P_DEFAULT,
                                  H5P_DEFAULT, H5P_DEFAULT);
    sc_set_str_array_attr(mod_root, "mod-order", &modality, 1);

    /* Create modality group and copy entire h5ad into it */
    char mod_path[512];
    snprintf(mod_path, sizeof(mod_path), "mod/%s", modality);
    hid_t mod_grp = H5Gcreate2(dst, mod_path, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);

    if (opts->verbose)
        fprintf(stderr, "  Copying h5ad → mod/%s\n", modality);

    rc = sc_copy_group_recursive(src, mod_grp, gzip);
    H5Gclose(mod_grp);
    H5Gclose(mod_root);

    /* Copy global obs from h5ad obs */
    if (rc == SC_OK && sc_has_group(src, "obs")) {
        hid_t src_obs = H5Gopen2(src, "obs", H5P_DEFAULT);
        hid_t dst_obs = sc_create_or_open_group(dst, "obs");
        sc_copy_group_recursive(src_obs, dst_obs, gzip);
        H5Gclose(dst_obs);
        H5Gclose(src_obs);
    } else {
        hid_t obs = sc_create_or_open_group(dst, "obs");
        sc_set_str_attr(obs, "encoding-type", "dataframe");
        sc_set_str_attr(obs, "encoding-version", "0.2.0");
        sc_set_str_attr(obs, "_index", "_index");
        H5Gclose(obs);
    }

    /* Create empty global obsm */
    {
        hid_t obsm = sc_create_or_open_group(dst, "obsm");
        sc_set_str_attr(obsm, "encoding-type", "dict");
        sc_set_str_attr(obsm, "encoding-version", "0.1.0");
        H5Gclose(obsm);
    }

    if (opts->verbose && rc == SC_OK)
        fprintf(stderr, "[scConvert] Done.\n");

    H5Fclose(dst);
    H5Fclose(src);
    return rc;
}
