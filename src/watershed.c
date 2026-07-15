/* Marker-based watershed by priority flooding.
 *
 * Reproduces skimage.segmentation.watershed (connectivity=1, compactness=0)
 * exactly, so that rHRG and pyHRG label the same pixels the same way. The
 * details that matter for that:
 *
 *   - the queue is ordered by (value, age), age being a strict insertion
 *     counter, so ties flood first-in-first-out rather than arbitrarily;
 *   - a pixel is labelled when it is pushed, not when it is popped;
 *   - neighbours are visited in raveled order (up, left, right, down),
 *     which is what decides plateaus.
 *
 * Data arrives row-major: R passes t(matrix), whose column-major storage
 * is byte-identical to numpy's row-major layout. Without that the neighbour
 * order would silently transpose and plateaus would resolve differently.
 *
 * Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.
 */
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

typedef struct { double value; R_xlen_t age; R_xlen_t index; } Elem;

/* binary min-heap ordered by (value, age) */
static int elem_lt(const Elem *a, const Elem *b) {
    if (a->value != b->value) return a->value < b->value;
    return a->age < b->age;
}

static void heap_push(Elem *h, R_xlen_t *n, Elem e) {
    R_xlen_t i = (*n)++;
    h[i] = e;
    while (i > 0) {
        R_xlen_t p = (i - 1) / 2;
        if (elem_lt(&h[i], &h[p])) { Elem t = h[i]; h[i] = h[p]; h[p] = t; i = p; }
        else break;
    }
}

static Elem heap_pop(Elem *h, R_xlen_t *n) {
    Elem top = h[0];
    h[0] = h[--(*n)];
    R_xlen_t i = 0;
    for (;;) {
        R_xlen_t l = 2*i + 1, r = l + 1, m = i;
        if (l < *n && elem_lt(&h[l], &h[m])) m = l;
        if (r < *n && elem_lt(&h[r], &h[m])) m = r;
        if (m == i) break;
        Elem t = h[i]; h[i] = h[m]; h[m] = t; i = m;
    }
    return top;
}

/* image, markers, mask: row-major, length nr*nc. Returns labels row-major. */
SEXP C_watershed(SEXP s_image, SEXP s_markers, SEXP s_mask,
                 SEXP s_nr, SEXP s_nc) {
    const int nr = asInteger(s_nr), nc = asInteger(s_nc);
    const R_xlen_t n = (R_xlen_t) nr * nc;
    const double *img = REAL(s_image);
    const int *mrk = INTEGER(s_markers);
    const int *msk = INTEGER(s_mask);

    SEXP s_out = PROTECT(allocVector(INTSXP, n));
    int *out = INTEGER(s_out);
    for (R_xlen_t i = 0; i < n; i++) out[i] = 0;

    Elem *heap = (Elem *) R_alloc(n + 1, sizeof(Elem));
    R_xlen_t hn = 0, age = 0;

    /* seed: markers enter in raveled order, as skimage's flatnonzero gives */
    for (R_xlen_t i = 0; i < n; i++) {
        if (mrk[i] != 0 && msk[i]) {
            out[i] = mrk[i];
            Elem e = { img[i], age++, i };
            heap_push(heap, &hn, e);
        }
    }

    const int dr[4] = { -1,  0,  0,  1 };   /* up, left, right, down: */
    const int dc[4] = {  0, -1,  1,  0 };   /* raveled order, as skimage */

    while (hn > 0) {
        Elem e = heap_pop(heap, &hn);
        const int r = (int) (e.index / nc), c = (int) (e.index % nc);
        const int lab = out[e.index];

        for (int k = 0; k < 4; k++) {
            const int rr = r + dr[k], cc = c + dc[k];
            if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
            const R_xlen_t j = (R_xlen_t) rr * nc + cc;
            if (!msk[j] || out[j] != 0) continue;
            out[j] = lab;                       /* label on push */
            Elem ne = { img[j], age++, j };
            heap_push(heap, &hn, ne);
        }
    }

    UNPROTECT(1);
    return s_out;
}
