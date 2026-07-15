#' Delineate tree crowns from a canopy height model
#'
#' Runs the whole pipeline: smooth, detect tops, optionally merge and screen
#' them, then grow crowns. Each stage is also available on its own
#' ([smooth_chm()], [detect_tops()], [merge_tops()], [screen_tops()],
#' [grow_crowns()]) if you want to inspect or replace one.
#'
#' @param chm Numeric matrix, or a `SpatRaster`, or a path to a raster.
#'   Rasters need the \pkg{terra} package.
#' @param smooth_ws,smooth_method Passed to [smooth_chm()].
#' @param hmin,detect_ws Passed to [detect_tops()].
#' @param merge_distance Numeric or `NULL`. If given, merge tops within this
#'   radius.
#' @param screen_hmin Numeric or `NULL`. If given, drop tops below this height.
#' @param quiet Logical. Suppress progress messages.
#' @param ... Passed to [grow_crowns()], e.g. `variance_thresh`,
#'   `conflict_rule`, `mask_thresh`.
#'
#' @return A list with:
#'   \describe{
#'     \item{crowns}{Integer matrix of crown ids, 0 for background. If `chm`
#'       was a `SpatRaster` or a path, a `SpatRaster` carrying the same CRS
#'       and extent.}
#'     \item{tops}{Two-column matrix of the tree tops actually used.}
#'     \item{n_contested}{Regions claimed by more than one crown.}
#'   }
#'
#' @examples
#' yy <- row(matrix(0, 60, 60)); xx <- col(matrix(0, 60, 60))
#' chm <- pmax(22 * exp(-((yy - 20)^2 + (xx - 20)^2) / 60),
#'             26 * exp(-((yy - 40)^2 + (xx - 40)^2) / 60))
#' res <- delineate(chm, hmin = 5, detect_ws = 5, variance_thresh = 6,
#'                  mask_thresh = 1, quiet = TRUE)
#' length(setdiff(unique(as.vector(res$crowns)), 0))
#' @export
delineate <- function(chm,
                      smooth_ws = 3L,
                      smooth_method = "median",
                      hmin = 2,
                      detect_ws = 3L,
                      merge_distance = NULL,
                      screen_hmin = NULL,
                      quiet = FALSE,
                      ...) {
  tmpl <- NULL
  if (inherits(chm, "SpatRaster") || is.character(chm)) {
    r <- read_chm(chm)
    tmpl <- r$template
    m <- r$chm
  } else {
    m <- as.matrix(chm)
  }
  say <- function(...) if (!quiet) message("  ", ...)

  say(sprintf("CHM %dx%d, range %.1f-%.1f m", nrow(m), ncol(m),
              min(m, na.rm = TRUE), max(m, na.rm = TRUE)))

  sm <- smooth_chm(m, ws = smooth_ws, method = smooth_method)
  say(sprintf("smooth: %s ws=%d", smooth_method, smooth_ws))

  tops <- detect_tops(sm, hmin = hmin, ws = detect_ws)
  say(sprintf("detect: %d tops (hmin=%s, ws=%d)", nrow(tops), hmin, detect_ws))

  if (!is.null(merge_distance)) {
    n0 <- nrow(tops); tops <- merge_tops(tops, merge_distance)
    say(sprintf("merge: %d -> %d (distance=%s)", n0, nrow(tops), merge_distance))
  }
  if (!is.null(screen_hmin)) {
    n0 <- nrow(tops); tops <- screen_tops(sm, tops, screen_hmin)
    say(sprintf("screen: %d -> %d (>=%s m)", n0, nrow(tops), screen_hmin))
  }

  crowns <- grow_crowns(sm, as_pixels(tops), ...)
  nc <- length(setdiff(unique(as.vector(crowns)), 0))
  merged <- nrow(tops) - nc
  say(sprintf("delineate: %d crowns%s, %d contested", nc,
              if (merged > 0) sprintf(", %d merged", merged) else "",
              attr(crowns, "n_contested")))

  out <- crowns
  if (!is.null(tmpl)) out <- .to_raster(crowns, tmpl)
  list(crowns = out, tops = tops, n_contested = attr(crowns, "n_contested"))
}
