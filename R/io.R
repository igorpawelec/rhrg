#' Read a CHM raster
#'
#' Thin wrapper over \pkg{terra}, kept out of the algorithm so that
#' [grow_crowns()] and friends work on plain matrices with no spatial
#' dependency at all.
#'
#' @param x A `SpatRaster` or a path to a raster file.
#' @param band Integer. Layer to read.
#'
#' @return A list with `chm` (numeric matrix) and `template` (a single-layer
#'   `SpatRaster` carrying the CRS, extent and resolution).
#'
#' @details
#' Nodata is left as-is rather than converted to `NA`: growing works on raw
#' values and a nodata of -9999 falls below any sensible `mask_thresh`.
#' Convert explicitly if your nodata is a plausible height.
#' @export
read_chm <- function(x, band = 1L) {
  .need_terra()
  r <- if (inherits(x, "SpatRaster")) x else terra::rast(x)
  r <- r[[band]]
  list(chm = terra::as.matrix(r, wide = TRUE), template = r)
}

#' Write a crown raster
#'
#' @param crowns Integer matrix of crown ids, or a `SpatRaster`.
#' @param filename Output path.
#' @param template A `SpatRaster` giving CRS and extent. Required when
#'   `crowns` is a matrix.
#' @param ... Passed to `terra::writeRaster`.
#' @return `filename`, invisibly.
#' @export
write_crowns <- function(crowns, filename, template = NULL, ...) {
  .need_terra()
  r <- if (inherits(crowns, "SpatRaster")) crowns else {
    if (is.null(template)) stop("template is required for a matrix", call. = FALSE)
    .to_raster(crowns, template)
  }
  terra::writeRaster(r, filename, overwrite = TRUE, ...)
  invisible(filename)
}

#' Convert crowns to polygons
#'
#' @param crowns Integer matrix or `SpatRaster` of crown ids.
#' @param template A `SpatRaster` giving CRS and extent, when `crowns` is a
#'   matrix.
#' @return A `SpatVector` of polygons, one per crown; background dropped.
#' @export
crowns_to_polygons <- function(crowns, template = NULL) {
  .need_terra()
  r <- if (inherits(crowns, "SpatRaster")) crowns else {
    if (is.null(template)) stop("template is required for a matrix", call. = FALSE)
    .to_raster(crowns, template)
  }
  r[r == 0] <- NA
  terra::as.polygons(r, dissolve = TRUE)
}

.to_raster <- function(m, template) {
  r <- terra::rast(template)
  terra::values(r) <- as.integer(t(m))   # terra fills row-major
  names(r) <- "crown"
  r
}

.need_terra <- function() {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("this needs the terra package:\n  install.packages(\"terra\")",
         call. = FALSE)
}
