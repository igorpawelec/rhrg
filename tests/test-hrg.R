library(rhrg)
ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

g <- function(n, cy, cx, h, s) {
  yy <- row(matrix(0, n, n)); xx <- col(matrix(0, n, n))
  h * exp(-((yy - cy)^2 + (xx - cx)^2) / (2 * s^2))
}

# a stand of well-separated trees
n <- 90
chm <- matrix(0, n, n); tops <- NULL
for (r in seq(15, 75, by = 30)) for (cc in seq(15, 75, by = 30)) {
  chm <- pmax(chm, g(n, r, cc, 20 + (r + cc) %% 8, 6))
  tops <- rbind(tops, c(r, cc))
}
tops <- matrix(as.integer(tops), ncol = 2)

cat("growing\n")
out <- grow_crowns(chm, tops, variance_thresh = 6, mask_thresh = 1)
ok(identical(dim(out), dim(chm)), "crown raster matches the CHM")
ok(is.integer(out), "crown ids are integers")
ok(all(out[chm < 1] == 0), "background stays 0")
ok(max(out) <= nrow(tops), "no id exceeds the tree count")

cat("reproducibility\n")
a <- grow_crowns(chm, tops, variance_thresh = 6, mask_thresh = 1)
b <- grow_crowns(chm, tops, variance_thresh = 6, mask_thresh = 1)
ok(identical(a, b), "repeated runs are identical")

cat("conflict rules\n")
for (rule in CONFLICT_RULES) {
  o <- grow_crowns(chm, tops, variance_thresh = 6, mask_thresh = 1,
                   conflict_rule = rule)
  ok(is.integer(o), paste(rule, "runs"))
}
ok(inherits(try(grow_crowns(chm, tops, conflict_rule = "bogus"), silent = TRUE),
            "try-error"), "bad conflict_rule rejected")

cat("one tree, two detected tops\n")
# a broad crown carrying two maxima: the classic over-detection
wide <- g(60, 30, 30, 24, 11) + g(60, 30, 24, 1.5, 3) + g(60, 30, 36, 2, 3)
two <- matrix(c(30L, 30L, 24L, 36L), ncol = 2)
m <- grow_crowns(wide, two, variance_thresh = 60, mask_thresh = 1,
                 conflict_rule = "height")
ok(length(setdiff(unique(as.vector(m)), 0)) == 1, "the two tops merge into one crown")
ok(attr(m, "n_contested") > 0, "the merge is recorded as contested")

p <- grow_crowns(wide, two, variance_thresh = 60, mask_thresh = 1,
                 protect_seeds = TRUE)
ok(length(setdiff(unique(as.vector(p)), 0)) == 2, "protect_seeds keeps both trees")

cat("seed order must not decide the winner\n")
fwd <- grow_crowns(wide, two, variance_thresh = 60, mask_thresh = 1)
rev <- grow_crowns(wide, two[2:1, , drop = FALSE], variance_thresh = 60,
                   mask_thresh = 1)
ok(identical(fwd == 2L, rev == 1L), "reversing the tops only swaps the ids")

cat("edge cases\n")
z <- grow_crowns(matrix(1, 20, 20), matrix(integer(0), 0, 2), mask_thresh = 0.5)
ok(max(z) == 0, "no seeds gives no crowns")
ok(inherits(try(grow_crowns(chm, matrix(c(999L, 999L), ncol = 2)), silent = TRUE),
            "try-error"), "out-of-raster seed rejected")
ok(is.integer(grow_crowns(chm, tops, variance_thresh = 6, mask_thresh = 1,
                          retry_rejected = TRUE)), "retry_rejected runs")
ok(is.integer(grow_crowns(chm, tops, variance_thresh = 6, mask_thresh = 1,
                          morpho_radius = 2L)), "morpho_radius runs")

cat("pipeline\n")
r <- delineate(chm, hmin = 5, detect_ws = 5, variance_thresh = 6,
               mask_thresh = 1, quiet = TRUE)
ok(is.list(r) && all(c("crowns", "tops", "n_contested") %in% names(r)),
   "delineate returns crowns, tops, n_contested")
ok(identical(dim(r$crowns), dim(chm)), "pipeline raster matches the CHM")
r2 <- delineate(chm, hmin = 5, detect_ws = 5, merge_distance = 4,
                screen_hmin = 8, variance_thresh = 6, mask_thresh = 1,
                quiet = TRUE)
ok(nrow(r2$tops) > 0, "optional merge/screen steps run")


# ── max_iters ─────────────────────────────────────────────────────────
#
# It defaulted to 200 until 0.3.0 and the bound bit: on chm_33_2012.tif 332
# of 492 crowns stopped there with candidates still queued, and the crown
# count read 132 against 63 once lifted. Boundaries barely moved (2.9% of
# the partition) because a truncated grow blocks merges rather than
# misplacing pixels, so the damage landed on the headline number.
cat("max_iters
")

set.seed(5); nn <- 90
yy <- row(matrix(0, nn, nn)); xx <- col(matrix(0, nn, nn))
mchm <- matrix(0, nn, nn)
for (i in 1:28) {
  r <- sample(6:(nn - 6), 1); c <- sample(6:(nn - 6), 1)
  mchm <- pmax(mchm, runif(1, 12, 26) *
                 exp(-((yy - r)^2 + (xx - c)^2) / (2 * runif(1, 3, 6)^2)))
}
msm <- smooth_chm(mchm, ws = 3, method = "median")
mtops <- as_pixels(merge_tops(detect_tops(msm, hmin = 5, ws = 5)))

ok(is.null(formals(grow_crowns)$max_iters),
   "max_iters defaults to NULL, meaning natural termination")

n_full <- length(setdiff(as.vector(
  grow_crowns(msm, mtops, variance_thresh = 20, mask_thresh = 1)), 0))

warned <- FALSE
cut <- withCallingHandlers(
  grow_crowns(msm, mtops, variance_thresh = 20, mask_thresh = 1, max_iters = 2),
  warning = function(w) {
    warned <<- grepl("max_iters", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
ok(warned, "a binding cap warns rather than truncating in silence")

n_cut <- length(setdiff(as.vector(cut), 0))
ok(n_cut > n_full,
   sprintf("a truncated grow leaves more, smaller crowns (%d against %d)",
           n_cut, n_full))

quiet <- TRUE
withCallingHandlers(
  invisible(grow_crowns(msm, mtops, variance_thresh = 20, mask_thresh = 1,
                        max_iters = 1e6)),
  warning = function(w) { quiet <<- FALSE; invokeRestart("muffleWarning") })
ok(quiet, "a cap that does not bind stays silent")

cat("\nall hrg tests passed\n")
