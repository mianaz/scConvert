/*
 * sc_dataframe.c — Streaming obs/var dataframe transfer between h5ad and h5seurat
 *
 * AnnData dataframes (obs, var) are stored as HDF5 groups:
 *   /obs/
 *     attrs: encoding-type="dataframe", _index="_index", column-order=[...]
 *     _index: vlen UTF-8 string dataset (cell barcodes / gene names)
 *     <col>:  numeric, string, or categorical (group with codes + categories)
 *
 * h5Seurat stores metadata as:
 *   /meta.data/       — cell metadata (obs equivalent)
 *   /assays/<assay>/features  — gene names
 *   /assays/<assay>/data      — expression matrix (CSC sparse)
 */

#include "sc_convert.h"

/* Forward declarations */
static int convert_categoricals_to_factors(hid_t meta_grp);
static int convert_factors_in_group(hid_t obs_grp);

/* ── Internal: copy all attributes from src to dst ──────────────────────────── */

static herr_t copy_attr_cb(hid_t src_loc, const char *name,
                            const H5A_info_t *ainfo, void *op_data)
{
    (void)ainfo;
    hid_t dst_loc = *(hid_t *)op_data;

    hid_t attr = H5Aopen(src_loc, name, H5P_DEFAULT);
    hid_t space = H5Aget_space(attr);
    hid_t type = H5Aget_type(attr);

    if (H5Aexists(dst_loc, name) > 0)
        H5Adelete(dst_loc, name);

    hid_t dst_attr = H5Acreate2(dst_loc, name, type, space,
                                 H5P_DEFAULT, H5P_DEFAULT);

    size_t sz = H5Aget_storage_size(attr);
    void *buf = malloc(sz > 0 ? sz : 1);
    if (!buf) {
        H5Aclose(dst_attr);
        H5Tclose(type);
        H5Sclose(space);
        H5Aclose(attr);
        return SC_ERR_HDF;
    }

    /* For variable-length strings, need special handling */
    if (H5Tget_class(type) == H5T_STRING && H5Tis_variable_str(type) > 0) {
        hid_t memtype = sc_create_vlen_str_type();
        int ndims = H5Sget_simple_extent_ndims(space);
        if (ndims == 0) {
            /* Scalar string */
            char *str = NULL;
            H5Aread(attr, memtype, &str);
            H5Awrite(dst_attr, memtype, &str);
            if (str) H5free_memory(str);
        } else {
            /* Array of strings */
            hsize_t dims;
            H5Sget_simple_extent_dims(space, &dims, NULL);
            char **strs = (char **)calloc(dims, sizeof(char *));
            if (strs) {
                H5Aread(attr, memtype, strs);
                H5Awrite(dst_attr, memtype, strs);
                for (hsize_t i = 0; i < dims; i++)
                    if (strs[i]) H5free_memory(strs[i]);
                free(strs);
            }
        }
        H5Tclose(memtype);
    } else {
        H5Aread(attr, type, buf);
        H5Awrite(dst_attr, type, buf);
    }

    free(buf);
    H5Aclose(dst_attr);
    H5Tclose(type);
    H5Sclose(space);
    H5Aclose(attr);
    return 0;
}

static int copy_all_attrs(hid_t src, hid_t dst) {
    hsize_t idx = 0;
    herr_t err = H5Aiterate2(src, H5_INDEX_NAME, H5_ITER_NATIVE,
                              &idx, copy_attr_cb, &dst);
    return (err < 0) ? SC_ERR_HDF : SC_OK;
}

/* ── Copy a categorical column (group with codes + categories) ──────────────── */

static int copy_categorical(hid_t src_grp, hid_t dst_parent,
                             const char *name, int gzip_level)
{
    hid_t src_cat = H5Gopen2(src_grp, name, H5P_DEFAULT);
    if (src_cat < 0) return SC_ERR_HDF;

    hid_t dst_cat = H5Gcreate2(dst_parent, name, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
    if (dst_cat < 0) { H5Gclose(src_cat); return SC_ERR_HDF; }

    /* Copy encoding attributes */
    copy_all_attrs(src_cat, dst_cat);

    int rc = SC_OK;

    /* Copy codes dataset */
    if (sc_has_dataset(src_cat, "codes")) {
        hid_t src_codes = H5Dopen2(src_cat, "codes", H5P_DEFAULT);
        rc = sc_copy_dataset_chunked(src_codes, dst_cat, "codes", gzip_level);
        /* Copy codes attributes */
        if (rc == SC_OK) {
            hid_t dst_codes = H5Dopen2(dst_cat, "codes", H5P_DEFAULT);
            copy_all_attrs(src_codes, dst_codes);
            H5Dclose(dst_codes);
        }
        H5Dclose(src_codes);
        if (rc != SC_OK) goto done;
    }

    /* Copy categories dataset */
    if (sc_has_dataset(src_cat, "categories")) {
        hid_t src_cats = H5Dopen2(src_cat, "categories", H5P_DEFAULT);
        rc = sc_copy_dataset_chunked(src_cats, dst_cat, "categories", gzip_level);
        if (rc == SC_OK) {
            hid_t dst_cats = H5Dopen2(dst_cat, "categories", H5P_DEFAULT);
            copy_all_attrs(src_cats, dst_cats);
            H5Dclose(dst_cats);
        }
        H5Dclose(src_cats);
    }

done:
    H5Gclose(dst_cat);
    H5Gclose(src_cat);
    return rc;
}

/* ── Stream an entire AnnData dataframe group ───────────────────────────────── */

int sc_stream_df_group(hid_t src_grp, hid_t dst_grp, int gzip_level) {
    /* Copy group-level attributes (encoding-type, _index, column-order, etc.) */
    copy_all_attrs(src_grp, dst_grp);

    /* Iterate members of the source group */
    hsize_t n_members;
    H5Gget_num_objs(src_grp, &n_members);

    for (hsize_t i = 0; i < n_members; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(src_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(src_grp, i);

        if (obj_type == H5G_DATASET) {
            /* Regular column or _index */
            hid_t src_dset = H5Dopen2(src_grp, name, H5P_DEFAULT);
            int rc = sc_copy_dataset_chunked(src_dset, dst_grp, name, gzip_level);
            if (rc == SC_OK) {
                /* Copy dataset-level attributes */
                hid_t dst_dset = H5Dopen2(dst_grp, name, H5P_DEFAULT);
                copy_all_attrs(src_dset, dst_dset);
                H5Dclose(dst_dset);
            }
            H5Dclose(src_dset);
            if (rc != SC_OK) return rc;
        }
        else if (obj_type == H5G_GROUP) {
            /* Categorical column or nested group — check encoding type */
            hid_t child = H5Gopen2(src_grp, name, H5P_DEFAULT);
            char enc[64] = "";
            sc_get_encoding_type(child, enc, sizeof(enc));
            H5Gclose(child);

            if (strcmp(enc, "categorical") == 0) {
                int rc = copy_categorical(src_grp, dst_grp, name, gzip_level);
                if (rc != SC_OK) return rc;
            } else {
                /* Generic nested group — recursive copy */
                hid_t src_child = H5Gopen2(src_grp, name, H5P_DEFAULT);
                hid_t dst_child = H5Gcreate2(dst_grp, name, H5P_DEFAULT,
                                              H5P_DEFAULT, H5P_DEFAULT);
                int rc = sc_copy_group_recursive(src_child, dst_child, gzip_level);
                H5Gclose(dst_child);
                H5Gclose(src_child);
                if (rc != SC_OK) return rc;
            }
        }
    }

    return SC_OK;
}

/* ── obs: h5ad → h5seurat ───────────────────────────────────────────────────── */

int sc_stream_obs_h5ad_to_h5seurat(hid_t src, hid_t dst, const char *assay) {
    (void)assay;
    /* h5ad: /obs (dataframe group)
     * h5seurat: /meta.data (same structure — just copy the group) */
    if (!sc_has_group(src, "obs")) return SC_OK;

    hid_t src_obs = H5Gopen2(src, "obs", H5P_DEFAULT);
    hid_t dst_meta = sc_create_or_open_group(dst, "meta.data");

    int rc = sc_stream_df_group(src_obs, dst_meta, SC_GZIP_LEVEL);

    /* Fix attributes: h5seurat uses "colnames" instead of "column-order",
     * and does not have AnnData encoding-type/version attributes */
    if (rc == SC_OK) {
        /* Rename column-order → colnames via HDF5's native H5Arename,
         * which preserves the attribute's on-disk type/encoding exactly
         * and avoids the read-deep-copy-write path that segfaulted on
         * Linux libhdf5 1.10 serial when re-converting vlen strings
         * (strlen on HDF5-private buffer layouts across lib versions). */
        if (H5Aexists(dst_meta, "column-order") > 0) {
            if (H5Aexists(dst_meta, "colnames") > 0)
                H5Adelete(dst_meta, "colnames");
            H5Arename(dst_meta, "column-order", "colnames");
        }
        /* Remove AnnData-specific attributes */
        if (H5Aexists(dst_meta, "encoding-type") > 0)
            H5Adelete(dst_meta, "encoding-type");
        if (H5Aexists(dst_meta, "encoding-version") > 0)
            H5Adelete(dst_meta, "encoding-version");

        /* Convert h5ad categoricals (categories/codes, 0-based) to
         * h5seurat factors (levels/values, 1-based) */
        convert_categoricals_to_factors(dst_meta);
    }

    H5Gclose(dst_meta);
    H5Gclose(src_obs);
    return rc;
}

/* ── Convert h5seurat factors (levels/values) to h5ad categoricals ──────────── */

static int convert_factors_in_group(hid_t obs_grp) {
    hsize_t n_members;
    H5Gget_num_objs(obs_grp, &n_members);

    for (hsize_t i = 0; i < n_members; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(obs_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(obs_grp, i);

        if (obj_type != H5G_GROUP) continue;

        hid_t child = H5Gopen2(obs_grp, name, H5P_DEFAULT);
        if (child < 0) continue;

        /* Check if this is a h5seurat factor (has levels + values) */
        if (!sc_has_dataset(child, "levels") || !sc_has_dataset(child, "values")) {
            H5Gclose(child);
            continue;
        }

        /* Read levels (string array) */
        hid_t lev_dset = H5Dopen2(child, "levels", H5P_DEFAULT);
        hid_t lev_space = H5Dget_space(lev_dset);
        hsize_t n_levels;
        H5Sget_simple_extent_dims(lev_space, &n_levels, NULL);

        hid_t str_type = sc_create_vlen_str_type();
        char **levels = (char **)calloc(n_levels, sizeof(char *));
        if (!levels) {
            H5Tclose(str_type);
            H5Sclose(lev_space);
            H5Dclose(lev_dset);
            H5Gclose(child);
            continue;
        }
        H5Dread(lev_dset, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, levels);
        H5Sclose(lev_space);
        H5Dclose(lev_dset);

        /* Read values (1-based int array) */
        hid_t val_dset = H5Dopen2(child, "values", H5P_DEFAULT);
        hid_t val_space = H5Dget_space(val_dset);
        hsize_t n_vals;
        H5Sget_simple_extent_dims(val_space, &n_vals, NULL);

        int *vals = (int *)malloc(n_vals * sizeof(int));
        if (!vals) {
            for (hsize_t j = 0; j < n_levels; j++)
                if (levels[j]) H5free_memory(levels[j]);
            free(levels);
            H5Tclose(str_type);
            H5Sclose(val_space);
            H5Dclose(val_dset);
            H5Gclose(child);
            continue;
        }
        H5Dread(val_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
        H5Sclose(val_space);
        H5Dclose(val_dset);

        H5Gclose(child);

        /* Convert 1-based → 0-based */
        for (hsize_t j = 0; j < n_vals; j++)
            vals[j] -= 1;

        /* Delete the old factor group and recreate as categorical */
        H5Ldelete(obs_grp, name, H5P_DEFAULT);

        hid_t cat_grp = H5Gcreate2(obs_grp, name, H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);

        /* Write categories */
        hsize_t lev_dims = n_levels;
        hid_t lev_sp = H5Screate_simple(1, &lev_dims, NULL);
        hid_t cat_dset = H5Dcreate2(cat_grp, "categories", str_type, lev_sp,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(cat_dset, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, levels);
        H5Dclose(cat_dset);
        H5Sclose(lev_sp);

        /* Write codes */
        hsize_t val_dims = n_vals;
        hid_t val_sp = H5Screate_simple(1, &val_dims, NULL);
        hid_t code_dset = H5Dcreate2(cat_grp, "codes", H5T_NATIVE_INT, val_sp,
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(code_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals);
        H5Dclose(code_dset);
        H5Sclose(val_sp);

        /* Set categorical encoding attributes */
        sc_set_str_attr(cat_grp, "encoding-type", "categorical");
        sc_set_str_attr(cat_grp, "encoding-version", "0.2.0");

        /* Set ordered=False (required by anndata) */
        {
            hid_t bool_space = H5Screate(H5S_SCALAR);
            hid_t bool_attr = H5Acreate2(cat_grp, "ordered", H5T_NATIVE_HBOOL,
                                          bool_space, H5P_DEFAULT, H5P_DEFAULT);
            hbool_t ordered = 0;
            H5Awrite(bool_attr, H5T_NATIVE_HBOOL, &ordered);
            H5Aclose(bool_attr);
            H5Sclose(bool_space);
        }

        H5Gclose(cat_grp);

        /* Cleanup */
        for (hsize_t j = 0; j < n_levels; j++)
            if (levels[j]) H5free_memory(levels[j]);
        free(levels);
        free(vals);

        /* Deleting/creating a member invalidates iteration — restart */
        H5Tclose(str_type);
        return convert_factors_in_group(obs_grp);
    }

    return SC_OK;
}

/* ── Convert h5ad categoricals (categories/codes) to h5seurat factors ──────── */

static int convert_categoricals_to_factors(hid_t meta_grp) {
    hsize_t n_members;
    H5Gget_num_objs(meta_grp, &n_members);

    for (hsize_t i = 0; i < n_members; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(meta_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(meta_grp, i);

        if (obj_type != H5G_GROUP) continue;

        hid_t child = H5Gopen2(meta_grp, name, H5P_DEFAULT);
        if (child < 0) continue;

        /* Check if this is an h5ad categorical (has categories + codes) */
        if (!sc_has_dataset(child, "categories") || !sc_has_dataset(child, "codes")) {
            H5Gclose(child);
            continue;
        }

        /* Read categories (string array) */
        hid_t cat_dset = H5Dopen2(child, "categories", H5P_DEFAULT);
        hid_t cat_space = H5Dget_space(cat_dset);
        hsize_t n_cats;
        H5Sget_simple_extent_dims(cat_space, &n_cats, NULL);

        hid_t str_type = sc_create_vlen_str_type();
        char **cats = (char **)calloc(n_cats, sizeof(char *));
        if (!cats) {
            H5Tclose(str_type);
            H5Sclose(cat_space);
            H5Dclose(cat_dset);
            H5Gclose(child);
            continue;
        }
        H5Dread(cat_dset, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, cats);
        H5Sclose(cat_space);
        H5Dclose(cat_dset);

        /* Read codes (0-based int array) */
        hid_t code_dset = H5Dopen2(child, "codes", H5P_DEFAULT);
        hid_t code_space = H5Dget_space(code_dset);
        hsize_t n_codes;
        H5Sget_simple_extent_dims(code_space, &n_codes, NULL);

        int *codes = (int *)malloc(n_codes * sizeof(int));
        if (!codes) {
            for (hsize_t j = 0; j < n_cats; j++)
                if (cats[j]) H5free_memory(cats[j]);
            free(cats);
            H5Tclose(str_type);
            H5Sclose(code_space);
            H5Dclose(code_dset);
            H5Gclose(child);
            continue;
        }
        H5Dread(code_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, codes);
        H5Sclose(code_space);
        H5Dclose(code_dset);

        H5Gclose(child);

        /* Convert 0-based → 1-based */
        for (hsize_t j = 0; j < n_codes; j++)
            codes[j] += 1;

        /* Delete the old categorical group and recreate as factor */
        H5Ldelete(meta_grp, name, H5P_DEFAULT);

        hid_t fac_grp = H5Gcreate2(meta_grp, name, H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);

        /* Write levels */
        hsize_t lev_dims = n_cats;
        hid_t lev_sp = H5Screate_simple(1, &lev_dims, NULL);
        hid_t lev_dset = H5Dcreate2(fac_grp, "levels", str_type, lev_sp,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(lev_dset, str_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, cats);
        H5Dclose(lev_dset);
        H5Sclose(lev_sp);

        /* Write values */
        hsize_t val_dims = n_codes;
        hid_t val_sp = H5Screate_simple(1, &val_dims, NULL);
        hid_t val_dset = H5Dcreate2(fac_grp, "values", H5T_NATIVE_INT, val_sp,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(val_dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, codes);
        H5Dclose(val_dset);
        H5Sclose(val_sp);

        H5Gclose(fac_grp);

        /* Cleanup */
        for (hsize_t j = 0; j < n_cats; j++)
            if (cats[j]) H5free_memory(cats[j]);
        free(cats);
        free(codes);

        /* Deleting/creating a member invalidates iteration — restart */
        H5Tclose(str_type);
        return convert_categoricals_to_factors(meta_grp);
    }

    return SC_OK;
}

/* ── obs: h5seurat → h5ad ───────────────────────────────────────────────────── */

int sc_stream_obs_h5seurat_to_h5ad(hid_t src, hid_t dst, const char *assay) {
    (void)assay;

    /* Fast path: H5Ocopy copies meta.data → obs as a single operation,
     * preserving all datasets, groups, and attributes without
     * decompress/recompress cycles. */
    if (sc_has_group(src, "meta.data")) {
        int rc = sc_copy_group_h5ocopy(src, "meta.data", dst, "obs");
        if (rc != SC_OK) {
            /* Fallback: manual copy */
            hid_t dst_obs = sc_create_or_open_group(dst, "obs");
            hid_t src_meta = H5Gopen2(src, "meta.data", H5P_DEFAULT);
            rc = sc_stream_df_group(src_meta, dst_obs, SC_GZIP_LEVEL);
            H5Gclose(src_meta);
            if (rc != SC_OK) { H5Gclose(dst_obs); return rc; }
            H5Gclose(dst_obs);
        }
    }

    hid_t dst_obs = H5Gopen2(dst, "obs", H5P_DEFAULT);
    if (dst_obs < 0)
        dst_obs = sc_create_or_open_group(dst, "obs");

    /* Convert h5seurat factors (levels/values, 1-based) to
     * h5ad categoricals (categories/codes, 0-based) */
    convert_factors_in_group(dst_obs);

    /* Fix attributes: rename colnames → column-order for h5ad.
     * Use H5Arename to preserve on-disk type exactly (see note above
     * in sc_stream_obs_h5ad_to_h5seurat about Linux libhdf5 vlen quirks). */
    if (H5Aexists(dst_obs, "colnames") > 0) {
        if (H5Aexists(dst_obs, "column-order") > 0)
            H5Adelete(dst_obs, "column-order");
        H5Arename(dst_obs, "colnames", "column-order");
    }

    /* Ensure AnnData encoding attributes */
    sc_set_str_attr(dst_obs, "encoding-type", "dataframe");
    sc_set_str_attr(dst_obs, "encoding-version", "0.2.0");

    /* If _index dataset doesn't exist, create from cell.names */
    if (!sc_has_dataset(dst_obs, "_index")) {
        if (sc_has_dataset(src, "cell.names")) {
            hid_t src_cn = H5Dopen2(src, "cell.names", H5P_DEFAULT);
            sc_copy_dataset_chunked(src_cn, dst_obs, "_index", SC_GZIP_LEVEL);
            H5Dclose(src_cn);
        }
    }
    /* Always set _index as scalar string (h5seurat stores it as 1-elem array,
     * but anndata requires a scalar string attribute) */
    sc_set_str_attr(dst_obs, "_index", "_index");

    /* Set encoding attrs on _index dataset (anndata expects string-array encoding) */
    if (sc_has_dataset(dst_obs, "_index")) {
        hid_t dst_idx = H5Dopen2(dst_obs, "_index", H5P_DEFAULT);
        sc_set_str_attr(dst_idx, "encoding-type", "string-array");
        sc_set_str_attr(dst_idx, "encoding-version", "0.2.0");
        H5Dclose(dst_idx);
    }

    /* AnnData spec requires encoding-type on every column dataset.
     * Factors already have "categorical"; _index has "string-array".
     * Numeric/boolean columns need "array". Iterate all dataset members
     * and set "array" on any that don't already have encoding-type. */
    {
        hsize_t n_members;
        H5Gget_num_objs(dst_obs, &n_members);
        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(dst_obs, i, name, sizeof(name));
            int obj_type = H5Gget_objtype_by_idx(dst_obs, i);
            if (obj_type == H5G_DATASET) {
                hid_t dset = H5Dopen2(dst_obs, name, H5P_DEFAULT);
                if (H5Aexists(dset, "encoding-type") <= 0) {
                    /* Check if it's a string or numeric dataset */
                    hid_t dtype = H5Dget_type(dset);
                    H5T_class_t cls = H5Tget_class(dtype);
                    if (cls == H5T_STRING) {
                        sc_set_str_attr(dset, "encoding-type", "string-array");
                    } else {
                        sc_set_str_attr(dset, "encoding-type", "array");
                    }
                    sc_set_str_attr(dset, "encoding-version", "0.2.0");
                    H5Tclose(dtype);
                }
                H5Dclose(dset);
            }
        }
    }

    H5Gclose(dst_obs);
    return SC_OK;
}

/* ── var: h5ad → h5seurat ───────────────────────────────────────────────────── */

int sc_stream_var_h5ad_to_h5seurat(hid_t src, hid_t dst, const char *assay) {
    if (!sc_has_group(src, "var")) return SC_OK;

    hid_t src_var = H5Gopen2(src, "var", H5P_DEFAULT);

    /* Create assay group path: /assays/<assay>/ */
    hid_t assays = sc_create_or_open_group(dst, "assays");
    hid_t assay_grp = sc_create_or_open_group(assays, assay);

    /* Copy _index as features */
    if (sc_has_dataset(src_var, "_index")) {
        hid_t src_idx = H5Dopen2(src_var, "_index", H5P_DEFAULT);
        sc_copy_dataset_chunked(src_idx, assay_grp, "features", SC_GZIP_LEVEL);
        H5Dclose(src_idx);
    }

    /* Copy var group for metadata preservation */
    hid_t dst_var = sc_create_or_open_group(dst, "var");
    sc_stream_df_group(src_var, dst_var, SC_GZIP_LEVEL);
    H5Gclose(dst_var);

    H5Gclose(assay_grp);
    H5Gclose(assays);
    H5Gclose(src_var);
    return SC_OK;
}

/* ── var: h5seurat → h5ad ───────────────────────────────────────────────────── */

int sc_stream_var_h5seurat_to_h5ad(hid_t src, hid_t dst, const char *assay) {
    hid_t dst_var = sc_create_or_open_group(dst, "var");

    /* Set dataframe encoding */
    sc_set_str_attr(dst_var, "encoding-type", "dataframe");
    sc_set_str_attr(dst_var, "encoding-version", "0.2.0");

    /* Get feature names from /assays/<assay>/features */
    char feat_path[512];
    snprintf(feat_path, sizeof(feat_path), "assays/%s/features", assay);

    if (sc_has_dataset(src, feat_path)) {
        hid_t src_feat = H5Dopen2(src, feat_path, H5P_DEFAULT);
        sc_copy_dataset_chunked(src_feat, dst_var, "_index", SC_GZIP_LEVEL);

        /* Set encoding attrs on _index */
        hid_t dst_idx = H5Dopen2(dst_var, "_index", H5P_DEFAULT);
        sc_set_str_attr(dst_idx, "encoding-type", "string-array");
        sc_set_str_attr(dst_idx, "encoding-version", "0.2.0");
        H5Dclose(dst_idx);

        H5Dclose(src_feat);
    }

    sc_set_str_attr(dst_var, "_index", "_index");

    /* If there's a /var group in h5seurat, copy its columns too */
    if (sc_has_group(src, "var")) {
        hid_t src_var = H5Gopen2(src, "var", H5P_DEFAULT);
        /* Copy individual columns, skipping _index */
        hsize_t n_members;
        H5Gget_num_objs(src_var, &n_members);
        for (hsize_t i = 0; i < n_members; i++) {
            char name[SC_MAX_NAME_LEN];
            H5Gget_objname_by_idx(src_var, i, name, sizeof(name));
            if (strcmp(name, "_index") == 0) continue;
            if (strcmp(name, "features") == 0) continue;

            int obj_type = H5Gget_objtype_by_idx(src_var, i);
            if (obj_type == H5G_DATASET) {
                hid_t src_dset = H5Dopen2(src_var, name, H5P_DEFAULT);
                if (!sc_has_dataset(dst_var, name)) {
                    sc_copy_dataset_chunked(src_dset, dst_var, name, SC_GZIP_LEVEL);
                    hid_t dst_dset = H5Dopen2(dst_var, name, H5P_DEFAULT);
                    copy_all_attrs(src_dset, dst_dset);
                    H5Dclose(dst_dset);
                }
                H5Dclose(src_dset);
            }
        }
        H5Gclose(src_var);
    }

    /* Set column-order if not present */
    if (H5Aexists(dst_var, "column-order") <= 0) {
        const char *empty = NULL;
        hsize_t zero = 0;
        hid_t space = H5Screate_simple(1, &zero, NULL);
        hid_t tid = sc_create_vlen_str_type();
        hid_t attr = H5Acreate2(dst_var, "column-order", tid, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, tid, &empty);
        H5Aclose(attr);
        H5Tclose(tid);
        H5Sclose(space);
    }

    H5Gclose(dst_var);
    return SC_OK;
}
