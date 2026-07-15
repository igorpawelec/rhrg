% Generated from roxygen comments in R/chm.R; regenerate with roxygen2::roxygenise()
\name{smooth_chm}
\alias{smooth_chm}
\title{Smooth a canopy height model}
\usage{
smooth_chm(chm, ws = 3L, method = c("median", "mean", "gaussian", "maximum"))
}
\arguments{
\item{chm}{Numeric matrix. Canopy height model.}
\item{ws}{Integer. Window size in pixels. Larger windows suppress more false
tops but also merge genuinely adjacent crowns; keep it below the narrowest
crown you want to preserve.}
\item{method}{One of \code{"median"}, \code{"mean"}, \code{"gaussian"},
\code{"maximum"}. \code{"median"} removes pits and spikes while keeping crown
edges sharp and is the usual choice. \code{"maximum"} flattens apexes into
plateaus and is rarely what you want before tree detection.}
}
\value{A numeric matrix, same dimensions as \code{chm}.}
\description{
A raw CHM carries pit noise and a ragged crown surface, both of which
\code{\link{detect_tops}} would read as extra trees. Smoothing trades a little
height accuracy for far fewer false tops.
}
\details{
Edges are handled by reflection, matching \code{scipy.ndimage}, so results
agree with pyHRG.

\code{NA} is not treated specially and will spread by \code{ws \%/\% 2}. Fill
nodata before smoothing if that matters.
}
\examples{
chm <- matrix(runif(400, 0, 30), 20, 20)
sm <- smooth_chm(chm, ws = 3)
}
