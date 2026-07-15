% Generated from roxygen comments in R/treetops.R; regenerate with roxygen2::roxygenise()
\name{detect_tops}
\alias{detect_tops}
\title{Detect tree tops as local maxima of a CHM}
\usage{
detect_tops(chm, hmin = 2, ws = 3L)
}
\arguments{
\item{chm}{Numeric matrix, ideally already smoothed by
\code{\link{smooth_chm}}. On a raw CHM this over-detects badly.}
\item{hmin}{Numeric. Minimum height for a pixel to be considered, in the units
of the CHM. Sets the floor between canopy and understorey.}
\item{ws}{Integer. Neighbourhood size, which acts as the minimum spacing
between tops. Over-detection is the safer error: \code{\link{delineate}}
merges surplus tops back together, but a top never detected is a tree lost.}
}
\value{
A two-column numeric matrix of (row, col) positions, 1-based. Zero detections
give a 0 x 2 matrix, not \code{NULL}, so callers can index \code{[, 1]}
unconditionally.
}
\description{
A pixel is a candidate when it equals the maximum of its \code{ws} x \code{ws}
neighbourhood and exceeds \code{hmin}. Adjacent candidates, which a flat apex
produces, are grouped into connected components and each component is reduced
to its height-weighted centre of mass, giving subpixel positions.
}
\references{
Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest.
\emph{Photogrammetric Engineering & Remote Sensing} 70(5), 589-604.
}
\examples{
chm <- matrix(0, 20, 20); chm[10, 10] <- 25
detect_tops(smooth_chm(chm, 3, "maximum"), hmin = 5, ws = 3)
}
