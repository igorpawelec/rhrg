# Plain-R tests: any error fails R CMD check. No testthat dependency, so
# these run anywhere the package installs.
library(rhrg)

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

cat("smoothing\n")
m <- matrix(1, 5, 5); m[3, 3] <- 99
ok(smooth_chm(m, 3, "median")[3, 3] == 1, "median kills a spike")
ok(identical(dim(smooth_chm(m, 3)), dim(m)), "shape preserved")
ok(all(smooth_chm(m, 1) == m), "ws=1 is a no-op")
for (meth in c("median", "mean", "gaussian", "maximum"))
  ok(all(is.finite(smooth_chm(m, 3, meth))), paste(meth, "runs"))
ok(inherits(try(smooth_chm(m, 0), silent = TRUE), "try-error"), "ws=0 rejected")
ok(inherits(try(smooth_chm(m, 3, "nope"), silent = TRUE), "try-error"),
   "bad method rejected")

cat("tree tops\n")
chm <- matrix(0, 30, 30); chm[10, 10] <- 25; chm[22, 20] <- 18
sm <- smooth_chm(chm, 3, "maximum")
t <- detect_tops(sm, hmin = 5, ws = 5)
ok(nrow(t) == 2, "two peaks found")
ok(ncol(t) == 2, "tops are (row, col)")
ok(is.matrix(t) && is.numeric(t), "tops are a numeric matrix")

e <- detect_tops(matrix(0, 10, 10), hmin = 5)
ok(identical(dim(e), c(0L, 2L)), "empty detection is 0x2, not 0")
invisible(e[, 1])   # must not error

ok(nrow(detect_tops(sm, hmin = 100, ws = 5)) == 0, "hmin filters everything")

cat("merging\n")
tp <- rbind(c(10, 10), c(10, 12), c(50, 50))
ok(nrow(merge_tops(tp, 5)) == 2, "close pair merges")
ok(nrow(merge_tops(tp, 0.5)) == 3, "distant tops kept")
ok(all(merge_tops(rbind(c(10, 10), c(10, 20)), 20)[1, ] == c(10, 15)),
   "merged top is the centroid")
ok(identical(dim(merge_tops(matrix(numeric(0), 0, 2))), c(0L, 2L)), "empty merge")

cat("screening\n")
ok(nrow(screen_tops(sm, t, 100)) == 0, "screen drops short tops")
ok(nrow(screen_tops(sm, t, 0)) == 2, "screen keeps tall tops")
ok(nrow(screen_tops(sm, rbind(c(-5, -5), c(500, 500)), 0)) == 0,
   "out-of-raster tops dropped, not wrapped")

cat("as_pixels\n")
ok(all(as_pixels(rbind(c(3.7, 4.2))) == c(3L, 4L)), "floors to integers")
ok(is.integer(as_pixels(rbind(c(3.7, 4.2)))), "returns integers")

cat("\nall chm/treetops tests passed\n")
