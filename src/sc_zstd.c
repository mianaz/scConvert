/*
 * sc_zstd.c — Zstandard codec bridge for the zarr reader/writer
 *
 * zarr v3 stores written by zarr-python 3 / anndata >= 0.12 default to the
 * `zstd` codec (usually wrapped in `sharding_indexed`). R has no CRAN Zstd
 * codec package, so when libzstd is available at build time (configure
 * defines HAVE_ZSTD) the package links it and exposes decompress/compress
 * entry points to R. Without libzstd the entry points report unavailability
 * and the R side falls back to the optional zstdlite package.
 *
 * Only built into the R shared library (BUILDING_R_PACKAGE); not part of the
 * standalone CLI binary.
 */
#ifdef BUILDING_R_PACKAGE

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_ZSTD
#include <zstd.h>
#endif

SEXP C_zstd_available(void) {
#ifdef HAVE_ZSTD
    return Rf_ScalarLogical(1);
#else
    return Rf_ScalarLogical(0);
#endif
}

SEXP C_zstd_decompress(SEXP raw_sexp) {
#ifdef HAVE_ZSTD
    if (TYPEOF(raw_sexp) != RAWSXP)
        Rf_error("zstd_decompress: input must be a raw vector");
    size_t clen = (size_t)XLENGTH(raw_sexp);
    const void *cbuf = RAW(raw_sexp);
    if (clen == 0) return Rf_allocVector(RAWSXP, 0);

    unsigned long long dlen = ZSTD_getFrameContentSize(cbuf, clen);
    if (dlen == ZSTD_CONTENTSIZE_ERROR)
        Rf_error("zstd_decompress: input is not a zstd frame");

    if (dlen == ZSTD_CONTENTSIZE_UNKNOWN) {
        /* Streaming decode with a growing buffer */
        ZSTD_DStream *ds = ZSTD_createDStream();
        if (!ds) Rf_error("zstd_decompress: cannot create stream");
        ZSTD_initDStream(ds);
        size_t cap = clen * 4 + 4096;
        unsigned char *out = (unsigned char *)malloc(cap);
        if (!out) { ZSTD_freeDStream(ds); Rf_error("zstd_decompress: out of memory"); }
        ZSTD_inBuffer in = { cbuf, clen, 0 };
        size_t pos = 0;
        while (in.pos < in.size) {
            if (cap - pos < ZSTD_DStreamOutSize()) {
                cap *= 2;
                unsigned char *tmp = (unsigned char *)realloc(out, cap);
                if (!tmp) { free(out); ZSTD_freeDStream(ds); Rf_error("zstd_decompress: out of memory"); }
                out = tmp;
            }
            ZSTD_outBuffer ob = { out, cap, pos };
            size_t rc = ZSTD_decompressStream(ds, &ob, &in);
            if (ZSTD_isError(rc)) {
                free(out); ZSTD_freeDStream(ds);
                Rf_error("zstd_decompress: %s", ZSTD_getErrorName(rc));
            }
            pos = ob.pos;
            if (rc == 0) break;
        }
        ZSTD_freeDStream(ds);
        SEXP res = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)pos));
        memcpy(RAW(res), out, pos);
        free(out);
        UNPROTECT(1);
        return res;
    }

    SEXP res = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)dlen));
    size_t rc = ZSTD_decompress(RAW(res), (size_t)dlen, cbuf, clen);
    if (ZSTD_isError(rc)) {
        UNPROTECT(1);
        Rf_error("zstd_decompress: %s", ZSTD_getErrorName(rc));
    }
    if (rc != (size_t)dlen) {
        SEXP trimmed = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)rc));
        memcpy(RAW(trimmed), RAW(res), rc);
        UNPROTECT(2);
        return trimmed;
    }
    UNPROTECT(1);
    return res;
#else
    (void)raw_sexp;
    Rf_error("scConvert was built without libzstd; install libzstd-dev "
             "(Debian/Ubuntu), zstd (Homebrew) or zstd (conda) and reinstall, "
             "or install the zstdlite R package");
    return R_NilValue;
#endif
}

SEXP C_zstd_compress(SEXP raw_sexp, SEXP level_sexp) {
#ifdef HAVE_ZSTD
    if (TYPEOF(raw_sexp) != RAWSXP)
        Rf_error("zstd_compress: input must be a raw vector");
    int level = Rf_asInteger(level_sexp);
    if (level == NA_INTEGER) level = 3;
    size_t dlen = (size_t)XLENGTH(raw_sexp);
    size_t bound = ZSTD_compressBound(dlen);
    unsigned char *buf = (unsigned char *)malloc(bound ? bound : 1);
    if (!buf) Rf_error("zstd_compress: out of memory");
    size_t rc = ZSTD_compress(buf, bound, RAW(raw_sexp), dlen, level);
    if (ZSTD_isError(rc)) {
        free(buf);
        Rf_error("zstd_compress: %s", ZSTD_getErrorName(rc));
    }
    SEXP res = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)rc));
    memcpy(RAW(res), buf, rc);
    free(buf);
    UNPROTECT(1);
    return res;
#else
    (void)raw_sexp; (void)level_sexp;
    Rf_error("scConvert was built without libzstd");
    return R_NilValue;
#endif
}

#endif /* BUILDING_R_PACKAGE */
