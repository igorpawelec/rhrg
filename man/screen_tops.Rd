% Generated from roxygen comments in R/treetops.R; regenerate with roxygen2::roxygenise()
\name{screen_tops}
\alias{screen_tops}
\title{Drop tree tops standing on low CHM values}
\usage{
screen_tops(chm, tops, hmin)
}
\arguments{
\item{chm}{Numeric matrix, the CHM the tops were detected on.}
\item{tops}{Two-column matrix of (row, col) positions.}
\item{hmin}{Numeric. Minimum height to keep.}
}
\value{
A two-column matrix. Tops outside the raster are dropped rather than wrapping
around it.
}
\description{
Useful as a second, stricter pass: detect at 2 m so the growing can see the
understorey, then keep only tops above 10 m.
}
\examples{
chm <- matrix(0, 20, 20); chm[10, 10] <- 25
screen_tops(chm, rbind(c(10, 10)), hmin = 10)
}
