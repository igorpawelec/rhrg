% Generated from roxygen comments in R/treetops.R; regenerate with roxygen2::roxygenise()
\name{as_pixels}
\alias{as_pixels}
\title{Convert subpixel tops to integer pixel positions}
\usage{
as_pixels(tops)
}
\arguments{
\item{tops}{Two-column matrix of (row, col) positions.}
}
\value{An integer two-column matrix suitable for \code{\link{grow_crowns}}.}
\description{
Convert subpixel tops to integer pixel positions.
}
\examples{
as_pixels(rbind(c(3.7, 4.2)))
}
