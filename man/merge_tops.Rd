% Generated from roxygen comments in R/treetops.R; regenerate with roxygen2::roxygenise()
\name{merge_tops}
\alias{merge_tops}
\title{Merge tree tops that lie close together}
\usage{
merge_tops(tops, distance = 5)
}
\arguments{
\item{tops}{Two-column matrix of (row, col) positions.}
\item{distance}{Numeric. Merge radius in pixels.}
}
\value{A two-column matrix with at most as many rows as \code{tops}.}
\description{
Tops closer than \code{distance} are almost always one tree seen through
several local maxima. Neighbours within the radius are grouped transitively
and each group is replaced by its centroid.
}
\details{
Grouping is transitive: a chain of tops each within \code{distance} of the next
collapses into one even if the endpoints are far apart. Keep the threshold well
below the typical crown diameter.
}
\examples{
merge_tops(rbind(c(10, 10), c(10, 12), c(50, 50)), distance = 5)
}
