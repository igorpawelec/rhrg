% Generated from roxygen comments in R/io.R; regenerate with roxygen2::roxygenise()
\name{write_crowns}
\alias{write_crowns}
\title{Write a crown raster}
\usage{
write_crowns(crowns, filename, template = NULL, ...)
}
\arguments{
\item{crowns}{Integer matrix of crown ids, or a \code{SpatRaster}.}
\item{filename}{Output path.}
\item{template}{A \code{SpatRaster} giving CRS and extent. Required when
\code{crowns} is a matrix.}
\item{...}{Passed to \code{terra::writeRaster}.}
}
\value{\code{filename}, invisibly.}
\description{
Write a crown raster.
}
