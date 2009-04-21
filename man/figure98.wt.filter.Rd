\name{figure98.wt.filter}
\alias{figure98.wt.filter}
\title{
  Plot a DWT Wavelet or Scaling Filter for Specific Levels.
}
\description{
  Plots a DWT Wavelet or Scaling Filter for specific levels similar to
  Figure 98 in \emph{Wavelet Methods for Time Series Analysis} by
  Percival and Walden (2000).
}
\usage{
figure98.wt.filter(filter, levels = NULL, wavelet = TRUE, y.normalize =
TRUE)
}
\arguments{
  \item{filter}{Either a 'wt.filter' object, a character string
    specifying a wavelet filter, or a numeric vector of wavelet
    coefficients. See 'help(wt.filter)' for acceptable filter 
          names.}
  \item{levels}{Number or vector indicating levels of filter to
    plot. See Details.}
  \item{wavelet}{A logical flag indicating whether to plot the wavelet
    (high pass) or scaling (low pass) filter.}
  \item{y.normalize}{A logical flag indicating whether to vertically
    normalize each level of the filter to the plotting space available
    for each level of the filter.}
}
\details{
    If a single number is specified for \code{levels}, then the filter
    of levels 1 through \code{levels} will be plotted.
    Otherwise, a vector will specify which levels of the wavelet
    coefficients will be plotted. If \code{levels} is not defined, it
    will default to the vector 1:7.

    The plotting space available for each level of the filter is
    dictated by the value of greatest magnitude of all the levels of
    the filter plotted. The vertical plotting space for each level will
    then be 2 times the absolute value of this magnitude.

    The label 'L sub j' on the horizontal axis varies from level to
    level, but 'L sub j' indicates the length of a filter at level
    j. The filter vector at a given level j is indexed from 0 to 'L sub
    j' - 1. 'L sub j' is defined by formula 96a in \emph{Wavelet Methods
      for Time Series Analysis} by Percival and Walden (2000).
}
\author{
    Kelvin Ma, kkym@u.washington.edu
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
  Series Analysis}, Cambridge University Press.  
}
\seealso{
  \code{\link{wt.filter}}
}
\examples{
# Plotting LA8 Wavelet Filter Coefficients Levels 1 through 7.
filter <- wt.filter()
figure98.wt.filter(filter)

# Alternatively
figure98.wt.filter("la8")

# Plotting D4 Scaling Filter Coefficients Levels 1, 3, and 5 and not
# vertically normalizing each level to its plotting region.
figure98.wt.filter("d4", levels = c(1,3,5), wavelet = FALSE, y.normalize
= FALSE)
}
\keyword{ts}
