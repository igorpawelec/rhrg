# Prove rHRG and pyHRG produce identical crowns.
#
# This is the claim the two packages rest on, so it is checked rather than
# asserted. It lives in tools/ and not tests/ because it needs Python and
# pyHRG installed, which CRAN machines do not have.
#
#   python3 tools/generate_pyhrg_reference.py     # writes tools/reference/
#   Rscript tools/cross_validate_against_pyhrg.R
#
# Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.

library(rhrg)
dir <- file.path("tools", "reference")
if (!dir.exists(dir))
  stop("run tools/generate_pyhrg_reference.py first", call. = FALSE)

meta <- utils::read.csv(file.path(dir, "meta.csv"), header = FALSE,
                        col.names = c("k", "vt", "rule", "n_contested"))
bad <- 0L; diff_px <- 0L; tot_px <- 0L; contested_bad <- 0L

for (i in seq_len(nrow(meta))) {
  k <- meta$k[i]
  rd <- function(f) as.matrix(utils::read.csv(
    file.path(dir, sprintf("%s_%d.csv", f, k)), header = FALSE))
  sm <- rd("sm"); ref <- rd("out")
  tops <- rd("tops"); if (ncol(tops) == 1L) tops <- t(tops)

  got <- grow_crowns(sm, tops, variance_thresh = meta$vt[i], mask_thresh = 1,
                     conflict_rule = as.character(meta$rule[i]))
  d <- sum(got != ref)
  tot_px <- tot_px + length(ref); diff_px <- diff_px + d
  if (d > 0L) {
    bad <- bad + 1L
    cat(sprintf("  case %d differs in %d px (vt=%s, rule=%s)\n",
                k, d, meta$vt[i], meta$rule[i]))
  }
  if (attr(got, "n_contested") != meta$n_contested[i]) contested_bad <- contested_bad + 1L
}

cat(sprintf("\ncases                : %d\n", nrow(meta)))
cat(sprintf("cases differing      : %d\n", bad))
cat(sprintf("contested mismatches : %d\n", contested_bad))
cat(sprintf("pixels differing     : %d of %s\n", diff_px, format(tot_px, big.mark = ",")))
if (bad == 0L && contested_bad == 0L) {
  cat("\nrHRG == pyHRG, bit for bit.\n")
} else {
  stop("rHRG and pyHRG disagree", call. = FALSE)
}
