/* Per-region statistics and region adjacency, in one pass over the pixels.
 *
 * Mirrors pyhrg.hrg._build_adjacency_and_stats, including the row-major
 * traversal order: floating-point addition is not associative, so summing
 * the same heights in a different order gives a slightly different variance,
 * which could tip a threshold comparison. Same order, same bits.
 *
 * Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.
 */
#include <R.h>
#include <Rinternals.h>

/* labels, chm: row-major, nr*nc. Returns list(n, mean, var, cy, cx, ea, eb). */
SEXP C_rag(SEXP s_labels, SEXP s_chm, SEXP s_nr, SEXP s_nc, SEXP s_nreg) {
    const int nr = asInteger(s_nr), nc = asInteger(s_nc);
    const int nreg = asInteger(s_nreg);
    const R_xlen_t n = (R_xlen_t) nr * nc;
    const int *lab = INTEGER(s_labels);
    const double *chm = REAL(s_chm);

    SEXP s_n    = PROTECT(allocVector(REALSXP, nreg + 1));
    SEXP s_mean = PROTECT(allocVector(REALSXP, nreg + 1));
    SEXP s_var  = PROTECT(allocVector(REALSXP, nreg + 1));
    SEXP s_cy   = PROTECT(allocVector(REALSXP, nreg + 1));
    SEXP s_cx   = PROTECT(allocVector(REALSXP, nreg + 1));
    double *rn = REAL(s_n), *rmean = REAL(s_mean), *rvar = REAL(s_var);
    double *rcy = REAL(s_cy), *rcx = REAL(s_cx);

    double *sum  = (double *) R_alloc(nreg + 1, sizeof(double));
    double *sum2 = (double *) R_alloc(nreg + 1, sizeof(double));
    double *sumr = (double *) R_alloc(nreg + 1, sizeof(double));
    double *sumc = (double *) R_alloc(nreg + 1, sizeof(double));
    for (int i = 0; i <= nreg; i++)
        rn[i] = sum[i] = sum2[i] = sumr[i] = sumc[i] = 0.0;

    for (int r = 0; r < nr; r++) {
        for (int c = 0; c < nc; c++) {
            const int L = lab[(R_xlen_t) r * nc + c];
            if (L <= 0) continue;
            const double v = chm[(R_xlen_t) r * nc + c];
            rn[L]   += 1.0;
            sum[L]  += v;
            sum2[L] += v * v;
            sumr[L] += r;
            sumc[L] += c;
        }
    }
    for (int i = 1; i <= nreg; i++) {
        if (rn[i] > 0) {
            rmean[i] = sum[i] / rn[i];
            rvar[i]  = sum2[i] / rn[i] - rmean[i] * rmean[i];
            if (rvar[i] < 0.0) rvar[i] = 0.0;   /* clamp fp noise */
            rcy[i] = sumr[i] / rn[i];
            rcx[i] = sumc[i] / rn[i];
        } else {
            rmean[i] = rvar[i] = rcy[i] = rcx[i] = 0.0;
        }
    }

    /* edges: right and down only, so each boundary is seen once per pixel
     * pair. Duplicates are kept: their count is the shared border length. */
    R_xlen_t cap = 1024, ne = 0;
    int *ea = (int *) malloc(cap * sizeof(int));
    int *eb = (int *) malloc(cap * sizeof(int));
    if (!ea || !eb) { free(ea); free(eb); error("rag: out of memory"); }

    for (int r = 0; r < nr; r++) {
        for (int c = 0; c < nc; c++) {
            const int L = lab[(R_xlen_t) r * nc + c];
            if (L <= 0) continue;
            for (int k = 0; k < 2; k++) {
                const int rr = r + (k == 1), cc = c + (k == 0);
                if (rr >= nr || cc >= nc) continue;
                const int M = lab[(R_xlen_t) rr * nc + cc];
                if (M <= 0 || M == L) continue;
                if (ne == cap) {
                    cap *= 2;
                    int *ta = (int *) realloc(ea, cap * sizeof(int));
                    int *tb = (int *) realloc(eb, cap * sizeof(int));
                    if (!ta || !tb) { free(ta ? ta : ea); free(tb ? tb : eb);
                                      error("rag: out of memory"); }
                    ea = ta; eb = tb;
                }
                ea[ne] = L < M ? L : M;
                eb[ne] = L < M ? M : L;
                ne++;
            }
        }
    }

    SEXP s_ea = PROTECT(allocVector(INTSXP, ne));
    SEXP s_eb = PROTECT(allocVector(INTSXP, ne));
    for (R_xlen_t i = 0; i < ne; i++) { INTEGER(s_ea)[i] = ea[i]; INTEGER(s_eb)[i] = eb[i]; }
    free(ea); free(eb);

    SEXP out = PROTECT(allocVector(VECSXP, 7));
    SET_VECTOR_ELT(out, 0, s_n);   SET_VECTOR_ELT(out, 1, s_mean);
    SET_VECTOR_ELT(out, 2, s_var); SET_VECTOR_ELT(out, 3, s_cy);
    SET_VECTOR_ELT(out, 4, s_cx);  SET_VECTOR_ELT(out, 5, s_ea);
    SET_VECTOR_ELT(out, 6, s_eb);
    SEXP nm = PROTECT(allocVector(STRSXP, 7));
    const char *nms[] = {"n","mean","var","cy","cx","ea","eb"};
    for (int i = 0; i < 7; i++) SET_STRING_ELT(nm, i, mkChar(nms[i]));
    setAttrib(out, R_NamesSymbol, nm);

    UNPROTECT(9);
    return out;
}
