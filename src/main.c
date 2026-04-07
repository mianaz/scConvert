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
#include <sys/wait.h>
#include <dirent.h>
#include <limits.h>
#include <errno.h>

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

/* Detect a CosMx SMI flat-file bundle. Accepts any directory that contains
 * at least 2 of the canonical NanoString CosMx CSV filenames. */
static int is_cosmx_dir(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)) return 0;
    DIR *d = opendir(path);
    if (!d) return 0;
    int hits_expr = 0, hits_meta = 0, hits_fov = 0, hits_tx = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (strstr(ent->d_name, "exprMat_file"))     hits_expr = 1;
        else if (strstr(ent->d_name, "metadata_file"))     hits_meta = 1;
        else if (strstr(ent->d_name, "fov_positions_file")) hits_fov  = 1;
        else if (strstr(ent->d_name, "tx_file"))           hits_tx   = 1;
    }
    closedir(d);
    return (hits_expr + hits_meta + hits_fov + hits_tx) >= 2;
}

/* Detect vendor raw formats that must be processed via the R backend. */
static int needs_r_delegation(const char *input) {
    if (has_extension(input, ".gef")) return 1;
    if (has_extension(input, ".cellbin.gef")) return 1;
    if (is_cosmx_dir(input)) return 1;
    return 0;
}

/* Check that Rscript is available on PATH without invoking a shell. */
static int rscript_available(void) {
    const char *path_env = getenv("PATH");
    if (!path_env) return 0;
    char *path_copy = strdup(path_env);
    if (!path_copy) return 0;
    int found = 0;
    char *saveptr = NULL;
    for (char *tok = strtok_r(path_copy, ":", &saveptr);
         tok != NULL;
         tok = strtok_r(NULL, ":", &saveptr)) {
        char candidate[PATH_MAX];
        int n = snprintf(candidate, sizeof(candidate), "%s/Rscript", tok);
        if (n > 0 && (size_t)n < sizeof(candidate) && access(candidate, X_OK) == 0) {
            found = 1;
            break;
        }
    }
    free(path_copy);
    return found;
}

/* Delegate conversion of vendor raw formats to the R backend via Rscript.
 *
 * We use execvp + fork rather than system() so filenames are passed as argv
 * strings that the shell never sees — no escaping, no injection surface.
 * stdout and stderr are inherited from the parent so users see R messages
 * directly.
 */
static int delegate_to_r(const sc_opts_t *opts) {
    if (!rscript_available()) {
        fprintf(stderr,
            "Error: '%s' requires the R backend but Rscript was not found on PATH.\n"
            "Install R and the scConvert package, then retry.\n",
            opts->input_path);
        return SC_ERR;
    }

    /* Absolutise input and output paths so the child does not depend on cwd.
     * realpath() fails for non-existent output files — absolutise the parent
     * directory and rejoin the basename instead. */
    char abs_input[PATH_MAX];
    if (realpath(opts->input_path, abs_input) == NULL) {
        fprintf(stderr, "Error: cannot resolve input path '%s'\n", opts->input_path);
        return SC_ERR;
    }

    char abs_output[PATH_MAX];
    {
        char out_copy[PATH_MAX];
        if (strlen(opts->output_path) >= sizeof(out_copy)) {
            fprintf(stderr, "Error: output path too long\n");
            return SC_ERR;
        }
        strncpy(out_copy, opts->output_path, sizeof(out_copy) - 1);
        out_copy[sizeof(out_copy) - 1] = '\0';

        /* Split into parent dir + basename */
        char *slash = strrchr(out_copy, '/');
        char parent_abs[PATH_MAX];
        const char *base;
        if (slash == NULL) {
            if (realpath(".", parent_abs) == NULL) {
                fprintf(stderr, "Error: cannot resolve cwd\n");
                return SC_ERR;
            }
            base = out_copy;
        } else {
            *slash = '\0';
            const char *parent = (out_copy[0] == '\0') ? "/" : out_copy;
            if (realpath(parent, parent_abs) == NULL) {
                fprintf(stderr, "Error: cannot resolve output parent directory '%s'\n", parent);
                return SC_ERR;
            }
            base = slash + 1;
        }
        int nw = snprintf(abs_output, sizeof(abs_output), "%s/%s", parent_abs, base);
        if (nw < 0 || (size_t)nw >= sizeof(abs_output)) {
            fprintf(stderr, "Error: output path too long after absolutisation\n");
            return SC_ERR;
        }
    }

    /* Build an R expression that calls scConvert() with explicit boolean args.
     * Filenames are NEVER interpolated into the R string — they come in as
     * separate argv entries after '--args', which R exposes via commandArgs(). */
    const char *expr =
        "args <- commandArgs(trailingOnly = TRUE); "
        "suppressPackageStartupMessages(library(scConvert)); "
        "scConvert::scConvert("
            "source = args[[1]], dest = args[[2]], "
            "overwrite = as.logical(args[[3]]), "
            "verbose = as.logical(args[[4]]))";

    if (opts->verbose) {
        fprintf(stderr,
                "[scConvert] Delegating to R backend: %s -> %s\n",
                abs_input, abs_output);
    }

    pid_t pid = fork();
    if (pid < 0) {
        fprintf(stderr, "Error: fork failed\n");
        return SC_ERR;
    }
    if (pid == 0) {
        /* Child: exec Rscript with a fixed argv. With `Rscript -e EXPR
         * arg1 arg2 ...` every following positional reaches the R script
         * via commandArgs(trailingOnly=TRUE) verbatim. NOT prepending an
         * `--args` marker — under -e mode Rscript would forward the
         * literal "--args" as the first positional. */
        char *child_argv[] = {
            (char *)"Rscript",
            (char *)"--vanilla",
            (char *)"-e",
            (char *)expr,
            abs_input,
            abs_output,
            (char *)(opts->overwrite ? "TRUE" : "FALSE"),
            (char *)(opts->verbose   ? "TRUE" : "FALSE"),
            NULL
        };
        execvp("Rscript", child_argv);
        /* execvp returned -> failure */
        fprintf(stderr, "Error: failed to exec Rscript: %s\n", strerror(errno));
        _exit(127);
    }

    /* Parent: wait for child */
    int status;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        fprintf(stderr, "Error: waitpid failed\n");
        return SC_ERR;
    }
    if (!WIFEXITED(status)) {
        fprintf(stderr, "Error: R backend terminated abnormally\n");
        return SC_ERR;
    }
    int code = WEXITSTATUS(status);
    if (code != 0) {
        fprintf(stderr, "Error: R backend exited with status %d\n", code);
        return SC_ERR;
    }
    return SC_OK;
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

    /* Vendor raw formats (.gef, .cellbin.gef, CosMx flat-file bundles) have
     * no C streaming path. Delegate to the R backend via Rscript. */
    if (needs_r_delegation(opts.input_path)) {
        int drc = delegate_to_r(&opts);
        return drc == SC_OK ? 0 : 1;
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
