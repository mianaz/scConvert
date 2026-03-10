/*
 * sc_modality.c — Modality ↔ assay name mapping for h5mu support
 *
 * Provides bidirectional mapping between Seurat assay names (RNA, ADT, ATAC)
 * and h5mu modality names (rna, prot, atac). Unmapped names are lowercased
 * for assay→modality, or preserved as-is for modality→assay.
 */

#include "sc_convert.h"
#include <ctype.h>

/* ── Static mapping tables ─────────────────────────────────────────────────── */

typedef struct {
    const char *from;
    const char *to;
} sc_name_map_t;

/* Assay name → modality name */
static const sc_name_map_t assay_to_mod[] = {
    {"RNA",     "rna"},
    {"ADT",     "prot"},
    {"ATAC",    "atac"},
    {"Spatial", "spatial"},
    {"HTO",     "hto"},
    {"SCT",     "sct"},
    {"Protein", "prot"},
    {"Peaks",   "atac"},
    {NULL, NULL}
};

/* Modality name → assay name (includes common aliases) */
static const sc_name_map_t mod_to_assay[] = {
    {"rna",     "RNA"},
    {"prot",    "ADT"},
    {"atac",    "ATAC"},
    {"spatial", "Spatial"},
    {"hto",     "HTO"},
    {"sct",     "SCT"},
    {"adt",     "ADT"},
    {"cite",    "ADT"},
    {"peaks",   "ATAC"},
    {NULL, NULL}
};

/* ── Thread-local buffer for fallback lowercase conversion ─────────────────── */

static _Thread_local char sc_mod_buf[SC_MAX_NAME_LEN];

/* ── Public API ────────────────────────────────────────────────────────────── */

/*
 * sc_modality_to_assay — Map h5mu modality name to Seurat assay name.
 *
 * Returns mapped name from the table, or the original modality name
 * (pointer to input) if not found in the table.
 */
const char *sc_modality_to_assay(const char *modality) {
    if (!modality) return "RNA";
    for (const sc_name_map_t *m = mod_to_assay; m->from; m++) {
        if (strcasecmp(modality, m->from) == 0)
            return m->to;
    }
    return modality;
}

/*
 * sc_assay_to_modality — Map Seurat assay name to h5mu modality name.
 *
 * Returns mapped name from the table, or a lowercased copy in a
 * thread-local buffer if not found.
 */
const char *sc_assay_to_modality(const char *assay) {
    if (!assay) return "rna";
    for (const sc_name_map_t *m = assay_to_mod; m->from; m++) {
        if (strcmp(assay, m->from) == 0)
            return m->to;
    }
    /* Fallback: lowercase the assay name */
    size_t len = strlen(assay);
    if (len >= sizeof(sc_mod_buf)) len = sizeof(sc_mod_buf) - 1;
    for (size_t i = 0; i < len; i++)
        sc_mod_buf[i] = (char)tolower((unsigned char)assay[i]);
    sc_mod_buf[len] = '\0';
    return sc_mod_buf;
}
