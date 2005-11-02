\name{squaredgain.wt.filter}
\alias{squaredgain.wt.filter}
\title{
  Plot Squared Gain Function for Wavelet Filter
}
\description{
  Plots squared gain function for a wavelet filter.
}
\usage{
squaredgain.wt.filter(filter, level = 1, N = NULL, draw.bands = TRUE,
wavelet = TRUE)
}
\arguments{
  \item{filter}{Either a 'wt.filter' object, a character string
    specifying a wavelet filter, or a numeric vector of wavelet
    coefficients. See 'help(wt.filter)' for acceptable filter names.}
  \item{level}{Level of wavelet filter. Applicable only if a character
    string is supplied in \code{filter}}
  \item{N}{Length of filter vector after padding the wavelet filter with
  zeroes. The zeroes are inserted after the filter values.}
  \item{draw.bands}{Draws dashed lines marking the lower and upper ends
    of the nominal pass-band for the wavelet filter.}
  \item{wavelet}{A logical flag indicating whether to plot the squared
    gain function for the wavelet (high pass) or scaling (low pass)
    filter.}
}
\details{
  If \code{N} is not specified, and the filter supplied is of length
  less than 1024, then the filter vector will be padded to a length of
  1024.  Otherwise, if the filter supplied is of length greater than
  1024, then the filter vector will be padded to the first 'power of 2'
  that is greater than the length of the filter supplied.
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
# Plotting the squared gain function of LA8 Wavelet Filter Coefficients.
filter <- wt.filter("la8")
squaredgain.wt.filter(filter)

# Plotting the squared gain function of LA8 Scaling Filter Coefficients.
squaredgain.wt.filter(filter, wavelet = FALSE)

# Plotting the squared gain function of Haar Wavelet Filter Coefficients
# without supplying a filter object.
squaredgain.wt.filter("haar")
}
\keyword{ts}
