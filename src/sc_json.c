/*
 * sc_json.c — Minimal JSON parser for zarr v2 metadata
 *
 * Only handles the fixed structures found in .zarray, .zattrs, .zgroup.
 * Not a general-purpose JSON parser.
 */

#include "sc_convert.h"
#include <ctype.h>
#include <errno.h>

/* ── Read entire file into malloc'd NUL-terminated string ───────────────────── */

char *sc_json_read_file(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long len = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (len <= 0) { fclose(f); return NULL; }
    char *buf = (char *)malloc((size_t)len + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t nr = fread(buf, 1, (size_t)len, f);
    fclose(f);
    buf[nr] = '\0';
    return buf;
}

/* ── Skip whitespace ────────────────────────────────────────────────────────── */

static const char *skip_ws(const char *p) {
    while (*p && isspace((unsigned char)*p)) p++;
    return p;
}

/* ── Skip a JSON value (string, number, object, array, literal) ─────────── */

static const char *skip_value(const char *p) {
    p = skip_ws(p);
    if (*p == '"') {
        /* string */
        p++;
        while (*p && *p != '"') {
            if (*p == '\\') p++;
            p++;
        }
        if (*p == '"') p++;
    } else if (*p == '{') {
        /* object */
        int depth = 1;
        p++;
        while (*p && depth > 0) {
            if (*p == '{') depth++;
            else if (*p == '}') depth--;
            else if (*p == '"') {
                p++;
                while (*p && *p != '"') {
                    if (*p == '\\') p++;
                    p++;
                }
            }
            p++;
        }
    } else if (*p == '[') {
        /* array */
        int depth = 1;
        p++;
        while (*p && depth > 0) {
            if (*p == '[') depth++;
            else if (*p == ']') depth--;
            else if (*p == '"') {
                p++;
                while (*p && *p != '"') {
                    if (*p == '\\') p++;
                    p++;
                }
            }
            p++;
        }
    } else {
        /* number, bool, null */
        while (*p && *p != ',' && *p != '}' && *p != ']' && !isspace((unsigned char)*p))
            p++;
    }
    return p;
}

/* ── Find the value position for a top-level key ───────────────────────────── */

const char *sc_json_find_key(const char *json, const char *key) {
    if (!json || !key) return NULL;
    const char *p = skip_ws(json);
    if (*p != '{') return NULL;
    p++;

    while (*p) {
        p = skip_ws(p);
        if (*p == '}') return NULL;
        if (*p != '"') return NULL;

        /* Extract key name */
        p++;
        const char *kstart = p;
        while (*p && *p != '"') {
            if (*p == '\\') p++;
            p++;
        }
        size_t klen = (size_t)(p - kstart);
        if (*p == '"') p++;

        /* skip colon */
        p = skip_ws(p);
        if (*p == ':') p++;
        p = skip_ws(p);

        /* Check if this key matches */
        if (klen == strlen(key) && memcmp(kstart, key, klen) == 0) {
            return p; /* pointer to value */
        }

        /* Skip value */
        p = skip_value(p);
        p = skip_ws(p);
        if (*p == ',') p++;
    }
    return NULL;
}

/* ── Extract a string value ─────────────────────────────────────────────────── */

int sc_json_get_string(const char *json, const char *key, char *buf, size_t buflen) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return -1;
    v = skip_ws(v);
    if (*v != '"') return -1;
    v++;
    size_t i = 0;
    while (*v && *v != '"' && i < buflen - 1) {
        if (*v == '\\') {
            v++;
            if (!*v) break;
        }
        buf[i++] = *v++;
    }
    buf[i] = '\0';
    return 0;
}

/* ── Extract an integer value ───────────────────────────────────────────────── */

int sc_json_get_int(const char *json, const char *key, int64_t *val) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return -1;
    v = skip_ws(v);
    char *end;
    errno = 0;
    long long ll = strtoll(v, &end, 10);
    if (end == v || errno) return -1;
    *val = (int64_t)ll;
    return 0;
}

/* ── Extract a double value ─────────────────────────────────────────────────── */

int sc_json_get_double(const char *json, const char *key, double *val) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return -1;
    v = skip_ws(v);
    /* Handle string-encoded special values like "NaN", "Infinity" */
    if (*v == '"') {
        char buf[32];
        sc_json_get_string(json, key, buf, sizeof(buf));
        if (strcmp(buf, "NaN") == 0) { *val = 0.0 / 0.0; return 0; }
        if (strcmp(buf, "Infinity") == 0) { *val = 1.0 / 0.0; return 0; }
        if (strcmp(buf, "-Infinity") == 0) { *val = -1.0 / 0.0; return 0; }
        return -1;
    }
    char *end;
    errno = 0;
    double d = strtod(v, &end);
    if (end == v || errno) return -1;
    *val = d;
    return 0;
}

/* ── Extract a boolean value ────────────────────────────────────────────────── */

int sc_json_get_bool(const char *json, const char *key, int *val) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return -1;
    v = skip_ws(v);
    if (strncmp(v, "true", 4) == 0)  { *val = 1; return 0; }
    if (strncmp(v, "false", 5) == 0) { *val = 0; return 0; }
    return -1;
}

/* ── Check if a value is null ───────────────────────────────────────────────── */

int sc_json_is_null(const char *json, const char *key) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return 1;
    v = skip_ws(v);
    return (strncmp(v, "null", 4) == 0);
}

/* ── Extract array of integers ──────────────────────────────────────────────── */

int sc_json_get_int_array(const char *json, const char *key,
                           int64_t *arr, int max_n) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return -1;
    v = skip_ws(v);
    if (*v != '[') return -1;
    v++;
    int n = 0;
    while (*v && n < max_n) {
        v = skip_ws(v);
        if (*v == ']') break;
        char *end;
        long long ll = strtoll(v, &end, 10);
        if (end == v) break;
        arr[n++] = (int64_t)ll;
        v = end;
        v = skip_ws(v);
        if (*v == ',') v++;
    }
    return n;
}

/* ── Extract array of strings ───────────────────────────────────────────────── */

int sc_json_get_str_array(const char *json, const char *key,
                           char ***out_arr, int *out_n) {
    const char *v = sc_json_find_key(json, key);
    if (!v) { *out_arr = NULL; *out_n = 0; return -1; }
    v = skip_ws(v);
    if (*v != '[') { *out_arr = NULL; *out_n = 0; return -1; }
    v++;

    /* First pass: count */
    int count = 0;
    const char *scan = v;
    while (*scan) {
        scan = skip_ws(scan);
        if (*scan == ']') break;
        if (*scan == '"') {
            count++;
            scan++;
            while (*scan && *scan != '"') {
                if (*scan == '\\') scan++;
                scan++;
            }
            if (*scan == '"') scan++;
        }
        scan = skip_ws(scan);
        if (*scan == ',') scan++;
    }

    char **arr = (char **)calloc((size_t)count, sizeof(char *));
    if (!arr) { *out_arr = NULL; *out_n = 0; return -1; }

    /* Second pass: extract */
    const char *p = v;
    for (int i = 0; i < count; i++) {
        p = skip_ws(p);
        if (*p != '"') break;
        p++;
        const char *start = p;
        while (*p && *p != '"') {
            if (*p == '\\') p++;
            p++;
        }
        size_t slen = (size_t)(p - start);
        arr[i] = (char *)malloc(slen + 1);
        if (!arr[i]) {
            for (int j = 0; j < i; j++) free(arr[j]);
            free(arr);
            *out_arr = NULL;
            *out_n = 0;
            return -1;
        }
        memcpy(arr[i], start, slen);
        arr[i][slen] = '\0';
        if (*p == '"') p++;
        p = skip_ws(p);
        if (*p == ',') p++;
    }

    *out_arr = arr;
    *out_n = count;
    return 0;
}

void sc_json_free_str_array(char **arr, int n) {
    if (!arr) return;
    for (int i = 0; i < n; i++) free(arr[i]);
    free(arr);
}

/* ── Extract a sub-object as raw JSON string ────────────────────────────────── */

int sc_json_get_object(const char *json, const char *key, char *buf, size_t buflen) {
    const char *v = sc_json_find_key(json, key);
    if (!v) return -1;
    v = skip_ws(v);
    if (*v != '{') return -1;
    const char *end = skip_value(v);
    size_t len = (size_t)(end - v);
    if (len >= buflen) len = buflen - 1;
    memcpy(buf, v, len);
    buf[len] = '\0';
    return 0;
}

/* ── Parse .zarray JSON into metadata struct ────────────────────────────────── */

int sc_json_parse_zarray(const char *json, sc_zarr_meta_t *meta) {
    memset(meta, 0, sizeof(*meta));
    meta->order = 'C';
    meta->compressor_level = 4;

    sc_json_get_int(json, "zarr_format", &meta->zarr_format);
    sc_json_get_string(json, "dtype", meta->dtype, sizeof(meta->dtype));

    /* shape */
    int64_t shape[8] = {0};
    int ndim = sc_json_get_int_array(json, "shape", shape, 8);
    if (ndim > 0) {
        meta->ndim = ndim;
        meta->shape[0] = shape[0];
        if (ndim > 1) meta->shape[1] = shape[1];
    }

    /* chunks */
    int64_t chunks[8] = {0};
    int nc = sc_json_get_int_array(json, "chunks", chunks, 8);
    if (nc > 0) {
        meta->chunks[0] = chunks[0];
        if (nc > 1) meta->chunks[1] = chunks[1];
    }

    /* compressor */
    char comp[256] = {0};
    if (sc_json_get_object(json, "compressor", comp, sizeof(comp)) == 0) {
        sc_json_get_string(comp, "id", meta->compressor_id, sizeof(meta->compressor_id));
        int64_t lvl = 4;
        if (sc_json_get_int(comp, "level", &lvl) == 0)
            meta->compressor_level = (int)lvl;
    } else if (sc_json_is_null(json, "compressor")) {
        meta->compressor_id[0] = '\0';
    }

    /* order */
    char order[4] = "C";
    sc_json_get_string(json, "order", order, sizeof(order));
    meta->order = order[0];

    /* filters — check for vlen-utf8 */
    if (!sc_json_is_null(json, "filters")) {
        const char *fv = sc_json_find_key(json, "filters");
        if (fv && *fv == '[') {
            /* Scan for "vlen-utf8" in the filters array */
            if (strstr(fv, "vlen-utf8"))
                meta->has_vlen_utf8 = 1;
        }
    }

    /* fill_value */
    sc_json_get_double(json, "fill_value", &meta->fill_value);

    return 0;
}

/* ── Write .zgroup ──────────────────────────────────────────────────────────── */

int sc_json_write_zgroup(const char *dir) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/.zgroup", dir);
    FILE *f = fopen(path, "w");
    if (!f) return SC_ERR_IO;
    fprintf(f, "{\"zarr_format\":2}\n");
    fclose(f);
    return SC_OK;
}

/* ── Write .zarray ──────────────────────────────────────────────────────────── */

int sc_json_write_zarray(const char *dir, const sc_zarr_meta_t *meta) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/.zarray", dir);
    FILE *f = fopen(path, "w");
    if (!f) return SC_ERR_IO;

    /* shape */
    fprintf(f, "{\n  \"zarr_format\": 2,\n  \"shape\": [");
    for (int i = 0; i < meta->ndim; i++) {
        if (i > 0) fprintf(f, ", ");
        fprintf(f, "%lld", (long long)meta->shape[i]);
    }
    fprintf(f, "],\n  \"chunks\": [");
    for (int i = 0; i < meta->ndim; i++) {
        if (i > 0) fprintf(f, ", ");
        fprintf(f, "%lld", (long long)meta->chunks[i]);
    }
    fprintf(f, "],\n");
    fprintf(f, "  \"dtype\": \"%s\",\n", meta->dtype);

    /* compressor */
    if (meta->compressor_id[0]) {
        fprintf(f, "  \"compressor\": {\"id\": \"%s\", \"level\": %d},\n",
                meta->compressor_id, meta->compressor_level);
    } else {
        fprintf(f, "  \"compressor\": null,\n");
    }

    /* fill_value */
    if (strcmp(meta->dtype, "|b1") == 0) {
        fprintf(f, "  \"fill_value\": false,\n");
    } else if (meta->dtype[1] == 'f') {
        fprintf(f, "  \"fill_value\": 0.0,\n");
    } else {
        fprintf(f, "  \"fill_value\": 0,\n");
    }

    fprintf(f, "  \"order\": \"%c\",\n", meta->order);

    /* filters */
    if (meta->has_vlen_utf8) {
        fprintf(f, "  \"filters\": [{\"id\": \"vlen-utf8\"}]\n");
    } else {
        fprintf(f, "  \"filters\": null\n");
    }

    fprintf(f, "}\n");
    fclose(f);
    return SC_OK;
}

/* ── Write .zattrs with raw JSON string ─────────────────────────────────────── */

int sc_json_write_zattrs_str(const char *dir, const char *json_str) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/.zattrs", dir);
    FILE *f = fopen(path, "w");
    if (!f) return SC_ERR_IO;
    fprintf(f, "%s\n", json_str);
    fclose(f);
    return SC_OK;
}
