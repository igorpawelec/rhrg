#' Smooth a canopy height model
#'
#' A raw CHM carries pit noise and a ragged crown surface, both of which
#' [detect_tops()] would read as extra trees. Smoothing trades a little
#' height accuracy for far fewer false tops.
#'
#' @param chm Numeric matrix. Canopy height model.
#' @param ws Integer. Window size in pixels. Larger windows suppress more
#'   false tops but also merge genuinely adjacent crowns; keep it below the
#'   narrowest crown you want to preserve.
#' @param method One of `"median"`, `"mean"`, `"gaussian"`, `"maximum"`.
#'   `"median"` removes pits and spikes while keeping crown edges sharp and
#'   is the usual choice. `"maximum"` flattens apexes into plateaus and is
#'   rarely what you want before tree detection.
#'
#' @return A numeric matrix, same dimensions as `chm`.
#'
#' @details
#' Edges are handled by reflection, matching `scipy.ndimage`, so results
#' agree with pyHRG.
#'
#' `NA` is not treated specially and will spread by `ws %/% 2`. Fill nodata
#' before smoothing if that matters.
#'
#' @examples
#' chm <- matrix(runif(400, 0, 30), 20, 20)
#' sm <- smooth_chm(chm, ws = 3)
#' @export
smooth_chm <- function(chm, ws = 3L, method = c("median", "mean", "gaussian", "maximum")) {
  method <- match.arg(method)
  chm <- as.matrix(chm)
  if (!is.numeric(chm)) stop("chm must be numeric", call. = FALSE)
  if (length(dim(chm)) != 2L) stop("chm must be a 2-D matrix", call. = FALSE)
  ws <- as.integer(ws)
  if (is.na(ws) || ws < 1L) stop("ws must be >= 1", call. = FALSE)
  # See detect_tops(): an even window sits half a pixel off, so the result
  # depends on the raster's orientation. Smoothing a 40x55 scene and its
  # mirror image differed by up to 10.5 m at ws=4. "gaussian" is exempt
  # because ws only sets sigma there and the kernel stays symmetric.
  if (method != "gaussian" && ws %% 2L == 0L)
    stop("ws must be odd for method = '", method, "', got ", ws,
         ". An even window has no centre pixel, so it sits half a pixel off ",
         "and the result depends on the raster's orientation. Use ", ws - 1L,
         " or ", ws + 1L, ".", call. = FALSE)
  if (ws == 1L) return(chm)

  # NA is skipped, not propagated. A window holding nothing but NA has no
  # statistic and stays NA.
  #
  # This matches pyHRG 0.5.0, which had to stop handing NaN to
  # scipy.ndimage: those filters are undefined on NaN rather than NaN-aware.
  # With the values 1..9 and one NaN, median_filter returned 9, 5 or 4
  # depending on where the NaN sat, and maximum_filter returned nan, 9 or 8.
  # Before that, R propagated NA instead and lost 3% of a real raster per
  # smoothing pass; neither behaviour was defensible.
  if (method == "gaussian") return(.gaussian_filter(chm, sigma = ws / 3))
  fun <- switch(method,
                median  = function(v) if (all(is.na(v))) NA_real_ else
                            stats::median(v, na.rm = TRUE),
                mean    = function(v) if (all(is.na(v))) NA_real_ else
                            mean(v, na.rm = TRUE),
                maximum = function(v) if (all(is.na(v))) NA_real_ else
                            max(v, na.rm = TRUE))
  .focal(chm, ws, fun)
}

# Reflected padding, as scipy.ndimage uses by default.
.pad_reflect <- function(m, p) {
  nr <- nrow(m); nc <- ncol(m)
  ri <- c(rev(seq_len(p)), seq_len(nr), nr - seq_len(p) + 1L)
  ci <- c(rev(seq_len(p)), seq_len(nc), nc - seq_len(p) + 1L)
  m[ri, ci, drop = FALSE]
}

.focal <- function(m, ws, fun) {
  p <- ws %/% 2L
  pm <- .pad_reflect(m, p)
  nr <- nrow(m); nc <- ncol(m)
  # stack every window offset, then reduce across the stack: ws^2 vectorised
  # passes instead of nr*nc calls into fun().
  cube <- array(NA_real_, c(nr, nc, ws * ws))
  k <- 1L
  for (dr in seq_len(ws) - 1L) {
    for (dc in seq_len(ws) - 1L) {
      cube[, , k] <- pm[seq_len(nr) + dr, seq_len(nc) + dc, drop = FALSE]
      k <- k + 1L
    }
  }
  out <- apply(cube, c(1, 2), fun)
  matrix(as.double(out), nr, nc)
}

.gaussian_filter <- function(m, sigma) {
  # Normalised convolution when NA is present: weight the kernel by which
  # cells exist, so a hole neither drags neighbours towards zero nor
  # spreads. Matches pyHRG's gaussian_filter(filled)/gaussian_filter(mask).
  if (anyNA(m)) {
    valid <- !is.na(m)
    filled <- m; filled[!valid] <- 0
    num <- .gaussian_filter(filled, sigma)
    den <- .gaussian_filter(matrix(as.double(valid), nrow(m), ncol(m)), sigma)
    out <- num / den
    out[den == 0] <- NA_real_
    return(out)
  }
  # separable, truncated at 4 sigma, as scipy.ndimage.gaussian_filter
  r <- as.integer(4 * sigma + 0.5)
  x <- -r:r
  k <- exp(-(x^2) / (2 * sigma^2)); k <- k / sum(k)
  pm <- .pad_reflect(m, r)
  nr <- nrow(m); nc <- ncol(m)
  tmp <- matrix(0, nr, nrow(t(pm)) - 0)
  # rows then columns
  acc <- matrix(0, nr, ncol(pm))
  for (i in seq_along(x)) acc <- acc + k[i] * pm[seq_len(nr) + r + x[i], , drop = FALSE]
  out <- matrix(0, nr, nc)
  for (i in seq_along(x)) out <- out + k[i] * acc[, seq_len(nc) + r + x[i], drop = FALSE]
  out
}
