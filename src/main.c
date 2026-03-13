/*
 * main.c — scConvert CLI entry point
 *
 * Usage:
 *   scconvert <input> <output> [options]
 *
 * Detects conversion direction from file extensions:
 *   .h5ad  → .h5seurat   (h5ad to h5seurat)
 *   .h5seurat → .h5ad    (h5seurat to h5ad)
 *
 * Options:
 *   --assay <name>   Assay name (default: RNA)
 *   --gzip <level>   Compression level 0-9 (default: 1)
 *   --overwrite      Overwrite existing output file
 *   --quiet          Suppress progress messages
 *   --version        Print version and exit
 *   --help           Print help and exit
 */

#include "sc_convert.h"
#include <unistd.h>
#include <sys/stat.h>

static void print_usage(const char *prog) {
    fprintf(stderr,
        "scConvert CLI v%s — Streaming single-cell format converter\n"
        "\n"
        "Usage: %s <input> <output> [options]\n"
        "\n"
        "Supported format pairs (auto-detected from extensions):\n"
        "  .h5ad     ↔ .h5seurat\n"
        "  .h5mu     ↔ .h5seurat   (multi-assay)\n"
        "  .h5mu     ↔ .h5ad       (single modality)\n"
        "  .loom     ↔ .h5seurat\n"
        "  .loom     ↔ .h5ad\n"
        "  .loom     ↔ .h5mu\n"
        "  .zarr     ↔ .h5seurat\n"
        "  .zarr     ↔ .h5ad\n"
        "  .zarr     ↔ .h5mu\n"
        "  .zarr     ↔ .loom\n"
        "\n"
        "For .rds, .soma use the R backend:\n"
        "  Rscript -e 'scConvert::scConvert_cli(\"in.h5ad\", \"out.rds\")'\n"
        "\n"
        "Options:\n"
        "  --assay <name>   Assay/modality name (default: RNA)\n"
        "  --gzip <level>   Compression level 0-9 (default: 1)\n"
        "  --overwrite      Overwrite existing output\n"
        "  --quiet          Suppress progress messages\n"
        "  --version        Print version\n"
        "  --help           Print this help\n"
        "\n"
        "Examples:\n"
        "  %s data.h5ad data.h5seurat\n"
        "  %s data.h5seurat data.h5ad --assay RNA --gzip 6\n"
        "  %s data.h5mu data.h5seurat\n"
        "  %s data.h5mu data.h5ad --assay RNA\n"
        "  %s data.loom data.h5seurat\n"
        "  %s data.h5ad data.loom\n",
        SC_VERSION_STRING, prog, prog, prog, prog, prog, prog, prog);
}

static int has_extension(const char *path, const char *ext) {
    size_t plen = strlen(path);
    size_t elen = strlen(ext);
    if (plen < elen) return 0;
    return strcasecmp(path + plen - elen, ext) == 0;
}

static int is_supported_format(const char *path) {
    return has_extension(path, ".h5ad") ||
           has_extension(path, ".h5seurat") ||
           has_extension(path, ".h5mu") ||
           has_extension(path, ".loom") ||
           has_extension(path, ".zarr");
}

static sc_direction_t detect_direction(const char *input, const char *output) {
    /* h5ad ↔ h5seurat */
    if (has_extension(input, ".h5ad") && has_extension(output, ".h5seurat"))
        return SC_H5AD_TO_H5SEURAT;
    if (has_extension(input, ".h5seurat") && has_extension(output, ".h5ad"))
        return SC_H5SEURAT_TO_H5AD;
    /* h5mu ↔ h5seurat */
    if (has_extension(input, ".h5mu") && has_extension(output, ".h5seurat"))
        return SC_H5MU_TO_H5SEURAT;
    if (has_extension(input, ".h5seurat") && has_extension(output, ".h5mu"))
        return SC_H5SEURAT_TO_H5MU;
    /* h5mu ↔ h5ad */
    if (has_extension(input, ".h5mu") && has_extension(output, ".h5ad"))
        return SC_H5MU_TO_H5AD;
    if (has_extension(input, ".h5ad") && has_extension(output, ".h5mu"))
        return SC_H5AD_TO_H5MU;
    /* loom ↔ h5seurat */
    if (has_extension(input, ".loom") && has_extension(output, ".h5seurat"))
        return SC_LOOM_TO_H5SEURAT;
    if (has_extension(input, ".h5seurat") && has_extension(output, ".loom"))
        return SC_H5SEURAT_TO_LOOM;
    /* loom ↔ h5ad */
    if (has_extension(input, ".loom") && has_extension(output, ".h5ad"))
        return SC_LOOM_TO_H5AD;
    if (has_extension(input, ".h5ad") && has_extension(output, ".loom"))
        return SC_H5AD_TO_LOOM;
    /* loom ↔ h5mu */
    if (has_extension(input, ".loom") && has_extension(output, ".h5mu"))
        return SC_LOOM_TO_H5MU;
    if (has_extension(input, ".h5mu") && has_extension(output, ".loom"))
        return SC_H5MU_TO_LOOM;
    /* zarr ↔ h5seurat */
    if (has_extension(input, ".zarr") && has_extension(output, ".h5seurat"))
        return SC_ZARR_TO_H5SEURAT;
    if (has_extension(input, ".h5seurat") && has_extension(output, ".zarr"))
        return SC_H5SEURAT_TO_ZARR;
    /* zarr ↔ h5ad */
    if (has_extension(input, ".zarr") && has_extension(output, ".h5ad"))
        return SC_ZARR_TO_H5AD;
    if (has_extension(input, ".h5ad") && has_extension(output, ".zarr"))
        return SC_H5AD_TO_ZARR;
    /* zarr ↔ h5mu */
    if (has_extension(input, ".zarr") && has_extension(output, ".h5mu"))
        return SC_ZARR_TO_H5MU;
    if (has_extension(input, ".h5mu") && has_extension(output, ".zarr"))
        return SC_H5MU_TO_ZARR;
    /* zarr ↔ loom */
    if (has_extension(input, ".zarr") && has_extension(output, ".loom"))
        return SC_ZARR_TO_LOOM;
    if (has_extension(input, ".loom") && has_extension(output, ".zarr"))
        return SC_LOOM_TO_ZARR;
    /* Unsupported pairs require the R backend */
    if (!is_supported_format(input) || !is_supported_format(output))
        return SC_DIRECTION_UNKNOWN;
    /* Default: try to infer from input */
    if (has_extension(input, ".h5ad"))
        return SC_H5AD_TO_H5SEURAT;
    if (has_extension(input, ".h5mu"))
        return SC_H5MU_TO_H5SEURAT;
    return SC_H5SEURAT_TO_H5AD;
}

int main(int argc, char **argv) {
    sc_opts_t opts = {
        .input_path  = NULL,
        .output_path = NULL,
        .assay_name  = "RNA",
        .gzip_level  = SC_GZIP_LEVEL,
        .verbose     = 1,
        .overwrite   = 0
    };

    /* Parse arguments */
    int positional = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
        if (strcmp(argv[i], "--version") == 0 || strcmp(argv[i], "-V") == 0) {
            printf("scConvert CLI v%s\n", SC_VERSION_STRING);
            return 0;
        }
        if (strcmp(argv[i], "--quiet") == 0 || strcmp(argv[i], "-q") == 0) {
            opts.verbose = 0;
            continue;
        }
        if (strcmp(argv[i], "--overwrite") == 0) {
            opts.overwrite = 1;
            continue;
        }
        if (strcmp(argv[i], "--assay") == 0 && i + 1 < argc) {
            opts.assay_name = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "--gzip") == 0 && i + 1 < argc) {
            opts.gzip_level = atoi(argv[++i]);
            if (opts.gzip_level < 0) opts.gzip_level = 0;
            if (opts.gzip_level > 9) opts.gzip_level = 9;
            continue;
        }
        if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
        /* Positional arguments */
        if (positional == 0)
            opts.input_path = argv[i];
        else if (positional == 1)
            opts.output_path = argv[i];
        else {
            fprintf(stderr, "Too many arguments\n");
            print_usage(argv[0]);
            return 1;
        }
        positional++;
    }

    if (!opts.input_path || !opts.output_path) {
        fprintf(stderr, "Error: input and output paths required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    /* Check input exists (file or directory for zarr) */
    struct stat in_stat;
    if (stat(opts.input_path, &in_stat) != 0) {
        fprintf(stderr, "Error: cannot read input: %s\n", opts.input_path);
        return 1;
    }

    /* Check output doesn't exist (unless --overwrite) */
    struct stat out_stat;
    if (!opts.overwrite && stat(opts.output_path, &out_stat) == 0) {
        fprintf(stderr, "Error: output exists: %s (use --overwrite)\n",
                opts.output_path);
        return 1;
    }
    /* If overwrite and output is a directory (zarr), remove it */
    if (opts.overwrite && stat(opts.output_path, &out_stat) == 0) {
        if (S_ISDIR(out_stat.st_mode))
            sc_rmdir_recursive(opts.output_path);
        else
            unlink(opts.output_path);
    }

    /* Detect direction */
    opts.direction = detect_direction(opts.input_path, opts.output_path);

    /* Disable HDF5 error printing (we handle errors ourselves) */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    int rc;
    switch (opts.direction) {
    case SC_H5AD_TO_H5SEURAT:
        rc = sc_h5ad_to_h5seurat(&opts);
        break;
    case SC_H5SEURAT_TO_H5AD:
        rc = sc_h5seurat_to_h5ad(&opts);
        break;
    case SC_H5MU_TO_H5SEURAT:
        rc = sc_h5mu_to_h5seurat(&opts);
        break;
    case SC_H5SEURAT_TO_H5MU:
        rc = sc_h5seurat_to_h5mu(&opts);
        break;
    case SC_H5MU_TO_H5AD:
        rc = sc_h5mu_to_h5ad(&opts);
        break;
    case SC_H5AD_TO_H5MU:
        rc = sc_h5ad_to_h5mu(&opts);
        break;
    case SC_LOOM_TO_H5SEURAT:
        rc = sc_loom_to_h5seurat(&opts);
        break;
    case SC_H5SEURAT_TO_LOOM:
        rc = sc_h5seurat_to_loom(&opts);
        break;
    case SC_LOOM_TO_H5AD:
        rc = sc_loom_to_h5ad(&opts);
        break;
    case SC_H5AD_TO_LOOM:
        rc = sc_h5ad_to_loom(&opts);
        break;
    case SC_LOOM_TO_H5MU:
        rc = sc_loom_to_h5mu(&opts);
        break;
    case SC_H5MU_TO_LOOM:
        rc = sc_h5mu_to_loom(&opts);
        break;
    case SC_ZARR_TO_H5SEURAT:
        rc = sc_zarr_to_h5seurat(&opts);
        break;
    case SC_H5SEURAT_TO_ZARR:
        rc = sc_h5seurat_to_zarr(&opts);
        break;
    case SC_ZARR_TO_H5AD:
        rc = sc_zarr_to_h5ad(&opts);
        break;
    case SC_H5AD_TO_ZARR:
        rc = sc_h5ad_to_zarr(&opts);
        break;
    case SC_ZARR_TO_H5MU:
        rc = sc_zarr_to_h5mu(&opts);
        break;
    case SC_H5MU_TO_ZARR:
        rc = sc_h5mu_to_zarr(&opts);
        break;
    case SC_ZARR_TO_LOOM:
        rc = sc_zarr_to_loom(&opts);
        break;
    case SC_LOOM_TO_ZARR:
        rc = sc_loom_to_zarr(&opts);
        break;
    case SC_DIRECTION_UNKNOWN:
        fprintf(stderr,
            "Error: this format pair requires the R backend.\n"
            "Use the R wrapper instead:\n"
            "  Rscript -e 'scConvert::sc_cli_convert(\"%s\", \"%s\")'\n"
            "\n"
            "The C binary supports: .h5ad, .h5seurat, .h5mu, .loom, .zarr\n"
            "The R backend adds:    .rds, .soma, .spatialdata.zarr\n",
            opts.input_path, opts.output_path);
        rc = SC_ERR_ARG;
        break;
    default:
        fprintf(stderr, "Error: unknown conversion direction\n");
        rc = SC_ERR_ARG;
        break;
    }

    if (rc != SC_OK) {
        fprintf(stderr, "Conversion failed (error code %d)\n", rc);
        /* Clean up partial output (file or directory) */
        struct stat out_st;
        if (stat(opts.output_path, &out_st) == 0) {
            if (S_ISDIR(out_st.st_mode))
                sc_rmdir_recursive(opts.output_path);
            else
                unlink(opts.output_path);
        }
        return 1;
    }

    return 0;
}
