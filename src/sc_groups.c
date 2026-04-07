/*
 * sc_groups.c — Streaming transfer of obsm, obsp, layers, uns, varm, varp
 *
 * These are mostly recursive group copies with format-specific handling:
 * - obsm: 2D embedding arrays (X_pca, X_umap, etc.)
 * - obsp: sparse graphs (connectivities, distances) as CSR groups
 * - layers: additional expression matrices as CSR groups
 * - uns: unstructured metadata (nested groups, scalars, arrays)
 */

#include "sc_convert.h"
#include <ctype.h>

/* ── H5Ocopy fast path for whole-group copy across files ────────────────────
 *
 * Copies an entire HDF5 object (group hierarchy, datasets, attributes) in a
 * single HDF5 call. Much faster than manual recursive iteration because HDF5
 * can optimize the internal transfer and skip decompress/recompress cycles.
 *
 * Use when NO structural transformation is needed (same layout in src & dst).
 * Falls back to sc_copy_group_recursive() if H5Ocopy fails.
 */
int sc_copy_group_h5ocopy(hid_t src_loc, const char *src_name,
                           hid_t dst_loc, const char *dst_name) {
    herr_t status = H5Ocopy(src_loc, src_name, dst_loc, dst_name,
                             H5P_DEFAULT, H5P_DEFAULT);
    return (status >= 0) ? SC_OK : SC_ERR_HDF;
}

/* ── Recursive group copy (fallback for custom transformations) ────────────── */

int sc_copy_group_recursive(hid_t src, hid_t dst, int gzip_level) {
    /* Copy group attributes */
    H5O_info2_t oinfo;
    H5Oget_info3(src, &oinfo, H5O_INFO_NUM_ATTRS);
    int num_attrs = (int)oinfo.num_attrs;
    for (int i = 0; i < num_attrs; i++) {
        hid_t attr = H5Aopen_by_idx(src, ".", H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, (hsize_t)i, H5P_DEFAULT, H5P_DEFAULT);
        char name[SC_MAX_NAME_LEN];
        H5Aget_name(attr, sizeof(name), name);

        hid_t space = H5Aget_space(attr);
        hid_t type = H5Aget_type(attr);

        if (H5Aexists(dst, name) > 0)
            H5Adelete(dst, name);

        hid_t dst_attr = H5Acreate2(dst, name, type, space,
                                     H5P_DEFAULT, H5P_DEFAULT);

        if (H5Tget_class(type) == H5T_STRING && H5Tis_variable_str(type) > 0) {
            int ndims = H5Sget_simple_extent_ndims(space);
            hid_t memtype = sc_create_vlen_str_type();
            if (ndims == 0) {
                char *str = NULL;
                H5Aread(attr, memtype, &str);
                H5Awrite(dst_attr, memtype, &str);
                if (str) H5free_memory(str);
            } else {
                hsize_t dims;
                H5Sget_simple_extent_dims(space, &dims, NULL);
                char **strs = (char **)calloc(dims, sizeof(char *));
                if (strs) {
                    H5Aread(attr, memtype, strs);
                    H5Awrite(dst_attr, memtype, strs);
                    for (hsize_t j = 0; j < dims; j++)
                        if (strs[j]) H5free_memory(strs[j]);
                    free(strs);
                }
            }
            H5Tclose(memtype);
        } else {
            size_t sz = H5Aget_storage_size(attr);
            if (sz > 0) {
                void *buf = malloc(sz);
                if (buf) {
                    H5Aread(attr, type, buf);
                    H5Awrite(dst_attr, type, buf);
                    free(buf);
                }
            }
        }

        H5Aclose(dst_attr);
        H5Tclose(type);
        H5Sclose(space);
        H5Aclose(attr);
    }

    /* Iterate group members */
    hsize_t n_members;
    H5Gget_num_objs(src, &n_members);

    for (hsize_t i = 0; i < n_members; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(src, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(src, i);

        if (obj_type == H5G_DATASET) {
            hid_t src_dset = H5Dopen2(src, name, H5P_DEFAULT);
            int rc = sc_copy_dataset_chunked(src_dset, dst, name, gzip_level);
            if (rc == SC_OK) {
                /* Copy dataset attributes */
                hid_t dst_dset = H5Dopen2(dst, name, H5P_DEFAULT);
                H5O_info2_t dinfo;
                H5Oget_info3(src_dset, &dinfo, H5O_INFO_NUM_ATTRS);
                int na = (int)dinfo.num_attrs;
                for (int a = 0; a < na; a++) {
                    hid_t sa = H5Aopen_by_idx(src_dset, ".", H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, (hsize_t)a, H5P_DEFAULT, H5P_DEFAULT);
                    char aname[SC_MAX_NAME_LEN];
                    H5Aget_name(sa, sizeof(aname), aname);

                    hid_t as = H5Aget_space(sa);
                    hid_t at = H5Aget_type(sa);

                    if (H5Aexists(dst_dset, aname) > 0)
                        H5Adelete(dst_dset, aname);

                    hid_t da = H5Acreate2(dst_dset, aname, at, as,
                                           H5P_DEFAULT, H5P_DEFAULT);

                    if (H5Tget_class(at) == H5T_STRING &&
                        H5Tis_variable_str(at) > 0) {
                        hid_t mt = sc_create_vlen_str_type();
                        int nd = H5Sget_simple_extent_ndims(as);
                        if (nd == 0) {
                            char *s = NULL;
                            H5Aread(sa, mt, &s);
                            H5Awrite(da, mt, &s);
                            if (s) H5free_memory(s);
                        } else {
                            hsize_t d;
                            H5Sget_simple_extent_dims(as, &d, NULL);
                            char **ss = calloc(d, sizeof(char*));
                            if (ss) {
                                H5Aread(sa, mt, ss);
                                H5Awrite(da, mt, ss);
                                for (hsize_t j = 0; j < d; j++)
                                    if (ss[j]) H5free_memory(ss[j]);
                                free(ss);
                            }
                        }
                        H5Tclose(mt);
                    } else {
                        size_t asz = H5Aget_storage_size(sa);
                        if (asz > 0) {
                            void *ab = malloc(asz);
                            if (ab) {
                                H5Aread(sa, at, ab);
                                H5Awrite(da, at, ab);
                                free(ab);
                            }
                        }
                    }

                    H5Aclose(da);
                    H5Tclose(at);
                    H5Sclose(as);
                    H5Aclose(sa);
                }
                H5Dclose(dst_dset);
            }
            H5Dclose(src_dset);
            if (rc != SC_OK) return rc;
        }
        else if (obj_type == H5G_GROUP) {
            hid_t src_child = H5Gopen2(src, name, H5P_DEFAULT);
            hid_t dst_child = H5Gcreate2(dst, name, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
            int rc = sc_copy_group_recursive(src_child, dst_child, gzip_level);
            H5Gclose(dst_child);
            H5Gclose(src_child);
            if (rc != SC_OK) return rc;
        }
    }

    return SC_OK;
}

/* ── obsm: embedding matrices ───────────────────────────────────────────────── */

int sc_stream_obsm(hid_t src, hid_t dst, sc_direction_t dir) {
    const char *src_path, *dst_path;

    if (dir == SC_H5AD_TO_H5SEURAT) {
        src_path = "obsm";
        dst_path = "reductions";
    } else {
        /* h5seurat→h5ad: reductions → obsm */
        src_path = "reductions";
        dst_path = "obsm";
    }

    if (!sc_has_group(src, src_path)) {
        /* Create empty target group with dict encoding for h5ad */
        if (dir == SC_H5SEURAT_TO_H5AD) {
            hid_t g = sc_create_or_open_group(dst, dst_path);
            sc_set_str_attr(g, "encoding-type", "dict");
            sc_set_str_attr(g, "encoding-version", "0.1.0");
            H5Gclose(g);
        }
        return SC_OK;
    }

    hid_t src_grp = H5Gopen2(src, src_path, H5P_DEFAULT);
    hid_t dst_grp = sc_create_or_open_group(dst, dst_path);

    if (dir == SC_H5SEURAT_TO_H5AD) {
        sc_set_str_attr(dst_grp, "encoding-type", "dict");
        sc_set_str_attr(dst_grp, "encoding-version", "0.1.0");
    }

    /* For h5ad→h5seurat: detect n_obs to decide if obsm needs transpose.
     * Standard h5ad: obsm dim0 == n_obs (C order).
     * R-generated h5ad: obsm dim1 == n_obs (Fortran order). */
    hsize_t n_obs = 0;
    if (dir == SC_H5AD_TO_H5SEURAT && sc_has_group(src, "obs")) {
        hid_t obs_grp = H5Gopen2(src, "obs", H5P_DEFAULT);
        if (sc_has_dataset(obs_grp, "_index")) {
            hid_t idx = H5Dopen2(obs_grp, "_index", H5P_DEFAULT);
            hid_t sp = H5Dget_space(idx);
            H5Sget_simple_extent_dims(sp, &n_obs, NULL);
            H5Sclose(sp);
            H5Dclose(idx);
        }
        H5Gclose(obs_grp);
    }

    /* Iterate and copy each member */
    hsize_t n;
    H5Gget_num_objs(src_grp, &n);

    for (hsize_t i = 0; i < n; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(src_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(src_grp, i);

        if (dir == SC_H5AD_TO_H5SEURAT) {
            /* h5ad obsm: datasets like X_pca → reductions/pca/cell.embeddings */
            if (obj_type == H5G_DATASET) {
                /* Strip X_ prefix for reduction name */
                const char *red_name = name;
                if (strncmp(name, "X_", 2) == 0)
                    red_name = name + 2;

                hid_t red_grp = sc_create_or_open_group(dst_grp, red_name);
                hid_t src_dset = H5Dopen2(src_grp, name, H5P_DEFAULT);

                /* Detect obsm dimension order and transpose if needed.
                 * Standard h5ad (Python): HDF5 dims [n_obs, n_comp] → transpose.
                 * R-generated h5ad: HDF5 dims [n_comp, n_obs] → direct copy.
                 * h5seurat needs: [n_comp, n_obs] (Fortran order). */
                hid_t emb_space = H5Dget_space(src_dset);
                int emb_ndims = H5Sget_simple_extent_ndims(emb_space);
                hsize_t emb_dims[2] = {0, 0};
                if (emb_ndims == 2)
                    H5Sget_simple_extent_dims(emb_space, emb_dims, NULL);
                H5Sclose(emb_space);

                if (emb_ndims == 2 && n_obs > 0 && emb_dims[0] == n_obs) {
                    /* Standard h5ad: dim0 is n_obs → transpose */
                    sc_copy_2d_transposed(src_dset, red_grp, "cell.embeddings",
                                           SC_GZIP_LEVEL);
                } else {
                    /* Already in Fortran order or can't detect → direct copy */
                    sc_copy_dataset_chunked(src_dset, red_grp, "cell.embeddings",
                                              SC_GZIP_LEVEL);
                }
                H5Dclose(src_dset);

                /* Set required h5seurat reduction attributes */
                sc_set_str_attr(red_grp, "active.assay", "RNA");
                {
                    /* Build key: e.g. "pca" → "PC_", "umap" → "UMAP_" */
                    char key[SC_MAX_NAME_LEN];
                    if (strcmp(red_name, "pca") == 0)
                        strcpy(key, "PC_");
                    else if (strcmp(red_name, "tsne") == 0)
                        strcpy(key, "tSNE_");
                    else {
                        /* Default: uppercase name + "_" */
                        size_t k;
                        for (k = 0; red_name[k] && k < sizeof(key) - 2; k++)
                            key[k] = (char)toupper((unsigned char)red_name[k]);
                        key[k] = '_';
                        key[k + 1] = '\0';
                    }
                    sc_set_str_attr(red_grp, "key", key);
                }
                /* global=0 means this reduction is not global */
                {
                    int global_val = 0;
                    hid_t s = H5Screate(H5S_SCALAR);
                    hid_t a = H5Acreate2(red_grp, "global", H5T_NATIVE_INT,
                                          s, H5P_DEFAULT, H5P_DEFAULT);
                    H5Awrite(a, H5T_NATIVE_INT, &global_val);
                    H5Aclose(a);
                    H5Sclose(s);
                }
                /* Create empty misc group (required by h5seurat) */
                {
                    hid_t misc = H5Gcreate2(red_grp, "misc", H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
                    H5Gclose(misc);
                }

                H5Gclose(red_grp);
            } else if (obj_type == H5G_GROUP) {
                /* Could be a sparse obsm entry — copy as-is */
                hid_t src_child = H5Gopen2(src_grp, name, H5P_DEFAULT);
                hid_t dst_child = sc_create_or_open_group(dst_grp, name);
                sc_copy_group_recursive(src_child, dst_child, SC_GZIP_LEVEL);
                H5Gclose(dst_child);
                H5Gclose(src_child);
            }
        } else {
            /* h5seurat reductions → h5ad obsm */
            if (obj_type == H5G_GROUP) {
                hid_t red_grp = H5Gopen2(src_grp, name, H5P_DEFAULT);

                /* Copy cell.embeddings as X_<name> */
                if (sc_has_dataset(red_grp, "cell.embeddings")) {
                    char dst_name[SC_MAX_NAME_LEN];
                    snprintf(dst_name, sizeof(dst_name), "X_%s", name);

                    hid_t src_emb = H5Dopen2(red_grp, "cell.embeddings",
                                              H5P_DEFAULT);
                    /* h5seurat cell.embeddings: HDF5 dims [n_comp, n_cells]
                     * (R Fortran order).
                     * h5ad obsm needs: [n_cells, n_comp] (C order).
                     * Must transpose [n_comp, n_cells] → [n_cells, n_comp]. */
                    sc_copy_2d_transposed(src_emb, dst_grp, dst_name,
                                            SC_GZIP_LEVEL);

                    /* Set array encoding on the obsm dataset */
                    hid_t dst_dset = H5Dopen2(dst_grp, dst_name, H5P_DEFAULT);
                    sc_set_str_attr(dst_dset, "encoding-type", "array");
                    sc_set_str_attr(dst_dset, "encoding-version", "0.2.0");
                    H5Dclose(dst_dset);

                    H5Dclose(src_emb);
                }
                H5Gclose(red_grp);
            }
        }
    }

    H5Gclose(dst_grp);
    H5Gclose(src_grp);
    return SC_OK;
}

/* ── obsp: graph matrices ───────────────────────────────────────────────────── */

int sc_stream_obsp(hid_t src, hid_t dst, sc_direction_t dir,
                    const char *assay) {
    const char *src_path, *dst_path;

    if (dir == SC_H5AD_TO_H5SEURAT) {
        src_path = "obsp";
        dst_path = "graphs";
    } else {
        src_path = "graphs";
        dst_path = "obsp";
    }

    if (!sc_has_group(src, src_path)) {
        /* Create empty group (no dict encoding for h5ad obsp) */
        hid_t g = sc_create_or_open_group(dst, dst_path);
        H5Gclose(g);
        return SC_OK;
    }

    hid_t src_grp = H5Gopen2(src, src_path, H5P_DEFAULT);
    hid_t dst_grp = sc_create_or_open_group(dst, dst_path);

    hsize_t n;
    H5Gget_num_objs(src_grp, &n);

    for (hsize_t i = 0; i < n; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(src_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(src_grp, i);

        if (obj_type == H5G_GROUP) {
            /* Fast path: H5Ocopy for sparse graphs (CSR structure preserved) */
            int copy_rc = sc_copy_group_h5ocopy(src_grp, name, dst_grp, name);
            if (copy_rc != SC_OK) {
                /* Fallback to manual copy */
                hid_t child = H5Gopen2(src_grp, name, H5P_DEFAULT);
                if (sc_has_dataset(child, "data") && sc_has_dataset(child, "indices")) {
                    hid_t dst_child = H5Gcreate2(dst_grp, name, H5P_DEFAULT,
                                                  H5P_DEFAULT, H5P_DEFAULT);
                    sc_stream_csr_copy(child, dst_child, SC_GZIP_LEVEL);
                    sc_copy_group_attrs(child, dst_child);
                    H5Gclose(dst_child);
                } else {
                    hid_t dst_child = sc_create_or_open_group(dst_grp, name);
                    sc_copy_group_recursive(child, dst_child, SC_GZIP_LEVEL);
                    H5Gclose(dst_child);
                }
                H5Gclose(child);
            }

            /* Fix graph matrix attributes for h5ad ↔ h5seurat direction */
            if (dir == SC_H5SEURAT_TO_H5AD) {
                hid_t dst_child = H5Gopen2(dst_grp, name, H5P_DEFAULT);
                /* Rename dims → shape */
                if (H5Aexists(dst_child, "dims") > 0) {
                    int64_t shape[2] = {0, 0};
                    hid_t dims_attr = H5Aopen(dst_child, "dims", H5P_DEFAULT);
                    hid_t dtype = H5Aget_type(dims_attr);
                    if (H5Tget_size(dtype) <= 4) {
                        int32_t s32[2];
                        H5Aread(dims_attr, H5T_NATIVE_INT32, s32);
                        shape[0] = s32[0]; shape[1] = s32[1];
                    } else {
                        H5Aread(dims_attr, H5T_NATIVE_INT64, shape);
                    }
                    H5Tclose(dtype);
                    H5Aclose(dims_attr);
                    H5Adelete(dst_child, "dims");
                    sc_set_int_array_attr(dst_child, "shape", shape, 2);
                }
                /* h5seurat stores graphs as CSC (R's dgCMatrix: column pointers),
                 * so label as csc_matrix — not csr_matrix — in h5ad obsp */
                sc_set_str_attr(dst_child, "encoding-type", "csc_matrix");
                sc_set_str_attr(dst_child, "encoding-version", "0.1.0");
                H5Gclose(dst_child);
            } else if (dir == SC_H5AD_TO_H5SEURAT) {
                hid_t dst_child = H5Gopen2(dst_grp, name, H5P_DEFAULT);
                /* Rename shape → dims */
                if (H5Aexists(dst_child, "shape") > 0) {
                    int64_t shape[2] = {0, 0};
                    hid_t shape_attr = H5Aopen(dst_child, "shape", H5P_DEFAULT);
                    hid_t dtype = H5Aget_type(shape_attr);
                    if (H5Tget_size(dtype) <= 4) {
                        int32_t s32[2];
                        H5Aread(shape_attr, H5T_NATIVE_INT32, s32);
                        shape[0] = s32[0]; shape[1] = s32[1];
                    } else {
                        H5Aread(shape_attr, H5T_NATIVE_INT64, shape);
                    }
                    H5Tclose(dtype);
                    H5Aclose(shape_attr);
                    H5Adelete(dst_child, "shape");
                    sc_set_int_array_attr(dst_child, "dims", shape, 2);
                }
                /* Remove h5ad-specific encoding attrs (not used by h5seurat) */
                if (H5Aexists(dst_child, "encoding-type") > 0)
                    H5Adelete(dst_child, "encoding-type");
                if (H5Aexists(dst_child, "encoding-version") > 0)
                    H5Adelete(dst_child, "encoding-version");
                /* Set assay.used so scLoadH5Seurat indexes this graph */
                if (assay != NULL)
                    sc_set_str_attr(dst_child, "assay.used", assay);
                H5Gclose(dst_child);
            }
        }
    }

    H5Gclose(dst_grp);
    H5Gclose(src_grp);
    return SC_OK;
}

/* ── varp: pairwise variable annotations (gene x gene) ─────────────────────── */
/*
 * h5seurat has no native varp slot. The scConvert R API stores varp under
 * misc/__varp__ (see tests/testthat/test-varp-roundtrip.R). This function
 * mirrors sc_stream_obsp but maps:
 *    h5ad:     /varp/{name}           (csc_matrix groups)
 *    h5seurat: /misc/__varp__/{name}  (dgCMatrix groups)
 *
 * varp matrices are gene x gene (symmetric in common use). No orientation
 * flip is needed — storage is identical to obsp (CSC, 0-based indices).
 */
int sc_stream_varp(hid_t src, hid_t dst, sc_direction_t dir) {
    const char *src_path, *dst_path;

    if (dir == SC_H5AD_TO_H5SEURAT) {
        src_path = "varp";
        dst_path = "misc/__varp__";
    } else {
        src_path = "misc/__varp__";
        dst_path = "varp";
    }

    if (!sc_has_group(src, src_path)) {
        /* Nothing to copy; ensure destination placeholder exists in h5ad. */
        if (dir == SC_H5SEURAT_TO_H5AD) {
            hid_t g = sc_create_or_open_group(dst, "varp");
            sc_set_str_attr(g, "encoding-type", "dict");
            sc_set_str_attr(g, "encoding-version", "0.1.0");
            H5Gclose(g);
        }
        return SC_OK;
    }

    /* For h5ad → h5seurat, parent misc/ must exist before creating misc/__varp__. */
    if (dir == SC_H5AD_TO_H5SEURAT) {
        if (!sc_has_group(dst, "misc")) {
            hid_t mg = H5Gcreate2(dst, "misc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (mg >= 0) H5Gclose(mg);
        }
    }

    hid_t src_grp = H5Gopen2(src, src_path, H5P_DEFAULT);
    hid_t dst_grp = sc_create_or_open_group(dst, dst_path);

    if (dir == SC_H5SEURAT_TO_H5AD) {
        sc_set_str_attr(dst_grp, "encoding-type", "dict");
        sc_set_str_attr(dst_grp, "encoding-version", "0.1.0");
    }

    hsize_t n;
    H5Gget_num_objs(src_grp, &n);

    for (hsize_t i = 0; i < n; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(src_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(src_grp, i);

        /* varp entries can be either groups (sparse CSR/CSC matrices) or
         * datasets (dense gene x gene arrays like attention weights).
         * Use H5Ocopy for both — it handles datasets and groups uniformly. */
        int copy_rc = sc_copy_group_h5ocopy(src_grp, name, dst_grp, name);
        if (copy_rc != SC_OK) {
            if (obj_type == H5G_GROUP) {
                hid_t child = H5Gopen2(src_grp, name, H5P_DEFAULT);
                if (sc_has_dataset(child, "data") && sc_has_dataset(child, "indices")) {
                    hid_t dst_child = H5Gcreate2(dst_grp, name, H5P_DEFAULT,
                                                  H5P_DEFAULT, H5P_DEFAULT);
                    sc_stream_csr_copy(child, dst_child, SC_GZIP_LEVEL);
                    sc_copy_group_attrs(child, dst_child);
                    H5Gclose(dst_child);
                } else {
                    hid_t dst_child = sc_create_or_open_group(dst_grp, name);
                    sc_copy_group_recursive(child, dst_child, SC_GZIP_LEVEL);
                    H5Gclose(dst_child);
                }
                H5Gclose(child);
            }
            /* For dataset children, if H5Ocopy failed there is nothing to
             * fall back to — skip and continue. */
            if (obj_type != H5G_GROUP) continue;
        }

        /* Attribute fixups apply only to sparse-matrix group children. For
         * dense-dataset children the source encoding-type is already "array"
         * and no dims/shape renaming is needed. */
        if (obj_type != H5G_GROUP) continue;

        /* Fix attributes for direction, mirroring obsp. */
        if (dir == SC_H5SEURAT_TO_H5AD) {
            hid_t dst_child = H5Gopen2(dst_grp, name, H5P_DEFAULT);
            if (H5Aexists(dst_child, "dims") > 0) {
                int64_t shape[2] = {0, 0};
                hid_t dims_attr = H5Aopen(dst_child, "dims", H5P_DEFAULT);
                hid_t dtype = H5Aget_type(dims_attr);
                if (H5Tget_size(dtype) <= 4) {
                    int32_t s32[2];
                    H5Aread(dims_attr, H5T_NATIVE_INT32, s32);
                    shape[0] = s32[0]; shape[1] = s32[1];
                } else {
                    H5Aread(dims_attr, H5T_NATIVE_INT64, shape);
                }
                H5Tclose(dtype);
                H5Aclose(dims_attr);
                H5Adelete(dst_child, "dims");
                sc_set_int_array_attr(dst_child, "shape", shape, 2);
            }
            /* Same convention as obsp: h5seurat dgCMatrix is column-major,
             * so label as csc_matrix in h5ad varp. */
            sc_set_str_attr(dst_child, "encoding-type", "csc_matrix");
            sc_set_str_attr(dst_child, "encoding-version", "0.1.0");
            H5Gclose(dst_child);
        } else if (dir == SC_H5AD_TO_H5SEURAT) {
            hid_t dst_child = H5Gopen2(dst_grp, name, H5P_DEFAULT);
            if (H5Aexists(dst_child, "shape") > 0) {
                int64_t shape[2] = {0, 0};
                hid_t shape_attr = H5Aopen(dst_child, "shape", H5P_DEFAULT);
                hid_t dtype = H5Aget_type(shape_attr);
                if (H5Tget_size(dtype) <= 4) {
                    int32_t s32[2];
                    H5Aread(shape_attr, H5T_NATIVE_INT32, s32);
                    shape[0] = s32[0]; shape[1] = s32[1];
                } else {
                    H5Aread(shape_attr, H5T_NATIVE_INT64, shape);
                }
                H5Tclose(dtype);
                H5Aclose(shape_attr);
                H5Adelete(dst_child, "shape");
                sc_set_int_array_attr(dst_child, "dims", shape, 2);
            }
            if (H5Aexists(dst_child, "encoding-type") > 0)
                H5Adelete(dst_child, "encoding-type");
            if (H5Aexists(dst_child, "encoding-version") > 0)
                H5Adelete(dst_child, "encoding-version");
            H5Gclose(dst_child);
        }
    }

    H5Gclose(dst_grp);
    H5Gclose(src_grp);
    return SC_OK;
}

/* ── layers: additional expression matrices ─────────────────────────────────── */

int sc_stream_layers(hid_t src, hid_t dst, sc_direction_t dir, int gzip_level) {
    if (!sc_has_group(src, "layers")) {
        if (dir == SC_H5SEURAT_TO_H5AD) {
            hid_t g = sc_create_or_open_group(dst, "layers");
            sc_set_str_attr(g, "encoding-type", "dict");
            sc_set_str_attr(g, "encoding-version", "0.1.0");
            H5Gclose(g);
        }
        return SC_OK;
    }

    hid_t src_grp = H5Gopen2(src, "layers", H5P_DEFAULT);
    hid_t dst_grp = sc_create_or_open_group(dst, "layers");

    if (dir == SC_H5SEURAT_TO_H5AD) {
        sc_set_str_attr(dst_grp, "encoding-type", "dict");
        sc_set_str_attr(dst_grp, "encoding-version", "0.1.0");
    }

    hsize_t n;
    H5Gget_num_objs(src_grp, &n);

    for (hsize_t i = 0; i < n; i++) {
        char name[SC_MAX_NAME_LEN];
        H5Gget_objname_by_idx(src_grp, i, name, sizeof(name));
        int obj_type = H5Gget_objtype_by_idx(src_grp, i);

        if (obj_type == H5G_GROUP) {
            /* Fast path: H5Ocopy for sparse layers (CSR structure preserved) */
            int copy_rc = sc_copy_group_h5ocopy(src_grp, name, dst_grp, name);
            if (copy_rc != SC_OK) {
                /* Fallback to manual copy */
                hid_t child = H5Gopen2(src_grp, name, H5P_DEFAULT);
                if (sc_has_dataset(child, "data") && sc_has_dataset(child, "indices")) {
                    hid_t dst_child = H5Gcreate2(dst_grp, name, H5P_DEFAULT,
                                                  H5P_DEFAULT, H5P_DEFAULT);
                    sc_stream_csr_copy(child, dst_child, gzip_level);
                    sc_copy_group_attrs(child, dst_child);
                    H5Gclose(dst_child);
                } else {
                    hid_t dst_child = sc_create_or_open_group(dst_grp, name);
                    sc_copy_group_recursive(child, dst_child, gzip_level);
                    H5Gclose(dst_child);
                }
                H5Gclose(child);
            }
        }
        else if (obj_type == H5G_DATASET) {
            hid_t src_dset = H5Dopen2(src_grp, name, H5P_DEFAULT);
            sc_copy_dataset_chunked(src_dset, dst_grp, name, gzip_level);
            H5Dclose(src_dset);
        }
    }

    H5Gclose(dst_grp);
    H5Gclose(src_grp);
    return SC_OK;
}

/* ── uns: unstructured metadata ─────────────────────────────────────────────── */

int sc_stream_uns(hid_t src, hid_t dst, sc_direction_t dir) {
    const char *src_path = "uns";
    const char *dst_path = "uns";

    /* For h5ad→h5seurat, uns goes to misc */
    if (dir == SC_H5AD_TO_H5SEURAT)
        dst_path = "misc";

    if (!sc_has_group(src, src_path)) {
        if (dir == SC_H5SEURAT_TO_H5AD) {
            hid_t g = sc_create_or_open_group(dst, dst_path);
            sc_set_str_attr(g, "encoding-type", "dict");
            sc_set_str_attr(g, "encoding-version", "0.1.0");
            H5Gclose(g);
        }
        return SC_OK;
    }

    /* Fast path: H5Ocopy copies the entire group hierarchy in one call,
     * avoiding manual iteration and decompress/recompress cycles */
    int rc = sc_copy_group_h5ocopy(src, src_path, dst, dst_path);

    if (rc == SC_OK && dir == SC_H5SEURAT_TO_H5AD) {
        /* Add dict encoding attributes required for h5ad */
        hid_t dst_grp = H5Gopen2(dst, dst_path, H5P_DEFAULT);
        sc_set_str_attr(dst_grp, "encoding-type", "dict");
        sc_set_str_attr(dst_grp, "encoding-version", "0.1.0");
        H5Gclose(dst_grp);
    }

    if (rc != SC_OK) {
        /* Fallback to recursive copy if H5Ocopy fails */
        hid_t src_grp = H5Gopen2(src, src_path, H5P_DEFAULT);
        hid_t dst_grp = sc_create_or_open_group(dst, dst_path);
        if (dir == SC_H5SEURAT_TO_H5AD) {
            sc_set_str_attr(dst_grp, "encoding-type", "dict");
            sc_set_str_attr(dst_grp, "encoding-version", "0.1.0");
        }
        rc = sc_copy_group_recursive(src_grp, dst_grp, SC_GZIP_LEVEL);
        H5Gclose(dst_grp);
        H5Gclose(src_grp);
    }

    return rc;
}

/* ── Ensure empty groups exist with dict encoding ───────────────────────────── */

int sc_ensure_empty_groups(hid_t file, sc_direction_t dir) {
    if (dir == SC_H5SEURAT_TO_H5AD) {
        const char *groups[] = {"obsm", "obsp", "varm", "varp", "layers", "uns"};
        for (int i = 0; i < 6; i++) {
            if (!sc_has_group(file, groups[i])) {
                hid_t g = sc_create_or_open_group(file, groups[i]);
                sc_set_str_attr(g, "encoding-type", "dict");
                sc_set_str_attr(g, "encoding-version", "0.1.0");
                H5Gclose(g);
            }
        }
    }
    return SC_OK;
}
