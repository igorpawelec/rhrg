/* Registration of the compiled routines. Copyright (C) 2025 Igor Pawelec. GPLv3. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_watershed(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rag(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_watershed", (DL_FUNC) &C_watershed, 5},
    {"C_rag",       (DL_FUNC) &C_rag,       5},
    {NULL, NULL, 0}
};

void attribute_visible R_init_rhrg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
