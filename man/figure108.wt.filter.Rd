\name{figure108.wt.filter}
\alias{figure108.wt.filter}
\title{
  Plot Multiple DWT Wavelet or Scaling Filters
}
\description{
  Plots multiple DWT Wavelet or Scaling Filters similar to
  Figure 108 in \emph{Wavelet Methods for Time Series Analysis} by
  Percival and Walden (2000).
}
\usage{
figure108.wt.filter(filter.objects, level = 1, l = NULL, wavelet = TRUE)
}
\arguments{
  \item{filter.objects}{List containing 'wt.filter' objects, character
    strings specifying a wavelet filter, or numeric vectors of wavelet
    coefficients. The list can contain a combination of 'wt.filter'
    objects, character strings, and numeric vectors. If only one filter
    is to be plotted, a single 'wt.filter' object, character string, or
    numeric vector may be supplied. See 'help(wt.filter)' for acceptable
    filter names.}
  \item{level}{If \code{filter.object} is only a single 'wt.filter'
    object, character string, or numeric vector, then the level of the
    filter can be specified in \code{level}.  Defaulted to 1.}
  \item{l}{Single integer representing the right hand limit of the
    horizontal axis. If unspecified, it will default to the length
    of the filter of greatest length given in
    \code{filter.object}. Cannot be less than the length of the filter
    of greatest length.}
  \item{wavelet}{A logical flag indicating whether to plot the wavelet
    (high pass) or scaling (low pass) filter.}
  }
\details{
    The plotting space available for each filter is dictated by the
    value of greatest magnitude of all the filters plotted. The vertical
    plotting space for each level will then be 2 times the absolute
    value of this magnitude.

    The filters are successively plotted in the order given in
    \code{filter.object}, where the first filter in \code{filter.object}
    is drawn at the top of the plot region, and the successive filters
    are plotted below.  
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
# Plotting the LA8 Wavelet Filter
filter <- wt.filter()
figure108.wt.filter(filter)

# Alternatively
figure108.wt.filter("la8")

# Plotting the Haar, D4, D6 Wavelet Filters
figure108.wt.filter(list("haar", "d4", "d6"))

# Plotting the Haar, D4, D6 Scaling Filters
figure108.wt.filter(list("haar", "d4", "d6"), wavelet = FALSE)

# Alternatively
haar <- wt.filter("haar")
d6 <- wt.filter("d6")
figure108.wt.filter(list(haar, "d4", d6), wavelet = FALSE)

# Adding an "made up" filter (represented by c(1,-1,1,-1)
figure108.wt.filter(list(haar, "d4", d6, c(1,-1,1,-1)), wavelet = FALSE)
}
\keyword{ts}
