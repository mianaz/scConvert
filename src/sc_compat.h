/*
 * sc_compat.h — R/CLI compatibility layer
 *
 * The C source files in src/ serve dual purpose:
 *   1. Compiled into the R package shared library (.so) via Makevars
 *   2. Compiled into a standalone CLI binary via Makefile
 *
 * For the R .so, CRAN requires REprintf() instead of fprintf(stderr, ...)
 * and snprintf() instead of sprintf(). This header provides a SC_MSG macro
 * that resolves to the correct function depending on the build target.
 *
 * When building for R (Makevars defines -DBUILDING_R_PACKAGE), SC_MSG uses
 * REprintf. When building the standalone CLI (Makefile), SC_MSG uses fprintf.
 */

#ifndef SC_COMPAT_H
#define SC_COMPAT_H

#ifdef BUILDING_R_PACKAGE
#include <R.h>
#include <Rinternals.h>
#define SC_MSG(...) REprintf(__VA_ARGS__)
#else
#include <stdio.h>
#define SC_MSG(...) fprintf(stderr, __VA_ARGS__)
#endif

#endif /* SC_COMPAT_H */
