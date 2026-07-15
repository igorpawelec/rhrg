% Generated from roxygen comments in R/io.R; regenerate with roxygen2::roxygenise()
\name{crowns_to_polygons}
\alias{crowns_to_polygons}
\title{Convert crowns to polygons}
\usage{
crowns_to_polygons(crowns, template = NULL)
}
\arguments{
\item{crowns}{Integer matrix or \code{SpatRaster} of crown ids.}
\item{template}{A \code{SpatRaster} giving CRS and extent, when \code{crowns}
is a matrix.}
}
\value{A \code{SpatVector} of polygons, one per crown; background dropped.}
\description{
Convert crowns to polygons.
}
