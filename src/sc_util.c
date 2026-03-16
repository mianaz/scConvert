/*
 * sc_util.c — HDF5 attribute helpers and utilities
 */

#include "sc_convert.h"

/* ── String attribute helpers ───────────────────────────────────────────────── */

hid_t sc_create_vlen_str_type(void) {
    hid_t tid = H5Tcopy(H5T_C_S1);
    H5Tset_size(tid, H5T_VARIABLE);
    H5Tset_cset(tid, H5T_CSET_UTF8);
    return tid;
}

int sc_set_str_attr(hid_t loc, const char *name, const char *value) {
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t tid = sc_create_vlen_str_type();
    hid_t attr;

    if (H5Aexists(loc, name) > 0)
        H5Adelete(loc, name);

    attr = H5Acreate2(loc, name, tid, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Tclose(tid); H5Sclose(space); return SC_ERR_HDF; }

    const char *ptr = value;
    herr_t err = H5Awrite(attr, tid, &ptr);

    H5Aclose(attr);
    H5Tclose(tid);
    H5Sclose(space);
    return (err < 0) ? SC_ERR_HDF : SC_OK;
}

int sc_set_str_array_attr(hid_t loc, const char *name,
                           const char **values, hsize_t n)
{
    hid_t space = H5Screate_simple(1, &n, NULL);
    hid_t tid = sc_create_vlen_str_type();
    hid_t attr;

    if (H5Aexists(loc, name) > 0)
        H5Adelete(loc, name);

    attr = H5Acreate2(loc, name, tid, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Tclose(tid); H5Sclose(space); return SC_ERR_HDF; }

    herr_t err = H5Awrite(attr, tid, values);

    H5Aclose(attr);
    H5Tclose(tid);
    H5Sclose(space);
    return (err < 0) ? SC_ERR_HDF : SC_OK;
}

int sc_set_int_array_attr(hid_t loc, const char *name,
                           const int64_t *values, hsize_t n)
{
    hid_t space = H5Screate_simple(1, &n, NULL);
    hid_t attr;

    if (H5Aexists(loc, name) > 0)
        H5Adelete(loc, name);

    /* Write as int32 — matches R/Python AnnData convention */
    int32_t *vals32 = (int32_t *)malloc(n * sizeof(int32_t));
    if (!vals32) { H5Sclose(space); return SC_ERR; }
    for (hsize_t i = 0; i < n; i++)
        vals32[i] = (int32_t)values[i];

    attr = H5Acreate2(loc, name, H5T_NATIVE_INT32, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { free(vals32); H5Sclose(space); return SC_ERR_HDF; }

    herr_t err = H5Awrite(attr, H5T_NATIVE_INT32, vals32);

    free(vals32);
    H5Aclose(attr);
    H5Sclose(space);
    return (err < 0) ? SC_ERR_HDF : SC_OK;
}

int sc_get_str_attr(hid_t loc, const char *name, char *buf, size_t buflen) {
    if (H5Aexists(loc, name) <= 0) {
        buf[0] = '\0';
        return SC_ERR;
    }

    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) return SC_ERR_HDF;

    hid_t tid = H5Aget_type(attr);
    char *tmp = NULL;

    if (H5Tis_variable_str(tid) > 0) {
        hid_t memtype = sc_create_vlen_str_type();
        H5Aread(attr, memtype, &tmp);
        H5Tclose(memtype);
        if (tmp) {
            strncpy(buf, tmp, buflen - 1);
            buf[buflen - 1] = '\0';
            H5free_memory(tmp);
        } else {
            buf[0] = '\0';
        }
    } else {
        size_t sz = H5Tget_size(tid);
        if (sz >= buflen) sz = buflen - 1;
        H5Aread(attr, tid, buf);
        buf[sz] = '\0';
    }

    H5Tclose(tid);
    H5Aclose(attr);
    return SC_OK;
}

int sc_get_str_array_attr(hid_t loc, const char *name,
                           char ***values, hsize_t *n)
{
    *values = NULL;
    *n = 0;

    if (H5Aexists(loc, name) <= 0)
        return SC_ERR;

    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) return SC_ERR_HDF;

    hid_t space = H5Aget_space(attr);
    int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims != 1) {
        H5Sclose(space);
        H5Aclose(attr);
        return SC_ERR;
    }

    hsize_t dims;
    H5Sget_simple_extent_dims(space, &dims, NULL);
    *n = dims;

    hid_t memtype = sc_create_vlen_str_type();
    char **raw = (char **)calloc(dims, sizeof(char *));
    if (!raw) {
        H5Tclose(memtype);
        H5Sclose(space);
        H5Aclose(attr);
        return SC_ERR;
    }
    H5Aread(attr, memtype, raw);

    /* Deep copy so we can free HDF5's memory */
    *values = (char **)calloc(dims, sizeof(char *));
    if (!*values) {
        for (hsize_t i = 0; i < dims; i++)
            if (raw[i]) H5free_memory(raw[i]);
        free(raw);
        H5Tclose(memtype);
        H5Sclose(space);
        H5Aclose(attr);
        return SC_ERR;
    }
    for (hsize_t i = 0; i < dims; i++) {
        if (raw[i]) {
            (*values)[i] = strdup(raw[i]);
            H5free_memory(raw[i]);
        } else {
            (*values)[i] = strdup("");
        }
    }
    free(raw);

    H5Tclose(memtype);
    H5Sclose(space);
    H5Aclose(attr);
    return SC_OK;
}

void sc_free_str_array(char **values, hsize_t n) {
    if (!values) return;
    for (hsize_t i = 0; i < n; i++)
        free(values[i]);
    free(values);
}

/* ── Group/dataset query helpers ────────────────────────────────────────────── */

int sc_has_group(hid_t loc, const char *name) {
    htri_t exists = H5Lexists(loc, name, H5P_DEFAULT);
    if (exists <= 0) return 0;

    H5O_info2_t info;
    if (H5Oget_info_by_name3(loc, name, &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0)
        return 0;
    return (info.type == H5O_TYPE_GROUP) ? 1 : 0;
}

int sc_has_dataset(hid_t loc, const char *name) {
    htri_t exists = H5Lexists(loc, name, H5P_DEFAULT);
    if (exists <= 0) return 0;

    H5O_info2_t info;
    if (H5Oget_info_by_name3(loc, name, &info, H5O_INFO_BASIC, H5P_DEFAULT) < 0)
        return 0;
    return (info.type == H5O_TYPE_DATASET) ? 1 : 0;
}

hid_t sc_create_or_open_group(hid_t loc, const char *name) {
    if (sc_has_group(loc, name))
        return H5Gopen2(loc, name, H5P_DEFAULT);
    return H5Gcreate2(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

int sc_get_encoding_type(hid_t loc, char *buf, size_t buflen) {
    return sc_get_str_attr(loc, "encoding-type", buf, buflen);
}

/* ── Copy all attributes from one HDF5 object to another ────────────────────── */

int sc_copy_group_attrs(hid_t src, hid_t dst) {
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
            hid_t memtype = sc_create_vlen_str_type();
            int ndims = H5Sget_simple_extent_ndims(space);
            if (ndims == 0) {
                char *str = NULL;
                H5Aread(attr, memtype, &str);
                H5Awrite(dst_attr, memtype, &str);
                if (str) H5free_memory(str);
            } else {
                hsize_t dims;
                H5Sget_simple_extent_dims(space, &dims, NULL);
                char **strs = (char **)calloc(dims, sizeof(char *));
                H5Aread(attr, memtype, strs);
                H5Awrite(dst_attr, memtype, strs);
                for (hsize_t j = 0; j < dims; j++)
                    if (strs[j]) H5free_memory(strs[j]);
                free(strs);
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
    return SC_OK;
}

