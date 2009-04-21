\name{plot.dwt.multiple}
\alias{plot.dwt.multiple}
\title{
  Plot Multiple DWT Objects
}
\description{
  Plot wavelet and scaling coefficients of multiple DWT objects.  
}
\usage{
\method{plot.dwt}{multiple}(x, levels = NULL, ylim = NULL, draw.dashed.lines =
TRUE, draw.level.labels = TRUE, col = c("red","blue"), ...)
}
\arguments{
  \item{x}{A list of \code{dwt} objects.}
  \item{levels}{Number, vector, or list of two vectors indicating range
    of levels to plot.  See details.}
  \item{ylim}{Vector specifying the lower and upper limits of the
    vertical y-axis of each plot.}
  \item{draw.dashed.lines}{Boolean indicating whether dashed lines
    should be drawn between levels of wavelet and scaling coefficients.}
  \item{draw.level.labels}{Boolean indicating whether the labels for the
     levels of wavelet and scaling coefficients should be drawn.}  
  \item{col}{Vector of length 2 for the alternating colors of the
    coefficients to be drawn.  The colors alternate by level of the
    wavelet or scaling coefficients drawn.  The first element of the
    vector is the color of the coefficients of the first level drawn.}
  \item{...}{Additional paramters that are acceptable arguments to the
    generic \code{plot} function}
}
\details{
    If a single number is specified for \code{levels}, then the wavelet
    coefficients of levels 1 through \code{levels} will be plotted.
    Otherwise, a vector or the first element of a list will specify
    which levels of the wavelet coefficients will be plotted.
    Unless specified in the second element of a list, only one level
    of scaling coefficients will be plotted and this level is equal
    to the highest level of the wavelet coefficients plotted. If a DWT
    object is defined for multiple time series, only the data pertaining
    to the first time series of the DWT object is plotted. Thus, only
    the wavelet coefficients and scaling coefficients of the first time
    series of the DWT objects will be plotted.

    For each \code{dwt}object in the list of \code{x},
    \code{plot.dwt.multiple} takes the coefficients of the \code{dwt}
    object and concatenates wavelet coefficients and scaling
    coefficients by levels specified in \code{levels}.  The wavelet
    coefficients will always be plotted preceding the scaling
    coefficients.

    This function allows users to visually examine differences in the
    DWT transform of a time series using different filters (different
    \code{dwt} objects).    

    For an example, see Figure 126 of \emph{Wavelet Methods for Time
      Series Analysis} by Percival and Walden (2000).
}
\author{
    Kelvin Ma, kkym@u.washington.edu
}  
\references{
    Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for
      Time Series Analysis}, Cambridge University Press.
}
\seealso{
  \code{\link{plot.dwt}}
}
\examples{
X <- rnorm(2048)
# Create DWT Object of X with the "la8" filter.
dwtobj1 <- dwt(X, filter = "la8")
# Create DWT Object of X with the "d4" filter.
dwtobj2 <- dwt(X, filter = "d4")
# Create DWT Object of X with the "haar" filter
dwtobj3 <- dwt(X, filter = "haar")
# Create DWT Object of X with the "c6" filter
dwtobj4 <- dwt(X, filter = "c6")

#Create list of dwt objects
dwtlist <- list(dwtobj1, dwtobj2, dwtobj3, dwtobj4)

# Plot the dwt objects and the wavelet coefficients of level 1 through 6
# and the scaling coefficients of level 6.  The first level drawn will
# be purple and the next level drawn will be gold.
plot.dwt.multiple(dwtlist, levels = 6, col = c("purple", "gold"))

# Plot the dwt objects and the wavelet coefficients of level 1, 3, and 5
# and scaling coefficients of level 2, and 4.
plot.dwt.multiple(dwtlist, levels = list(c(1,3,5), c(2,4)))

}
\keyword{ts}
