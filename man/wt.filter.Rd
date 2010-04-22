\name{wt.filter}
\alias{wt.filter}
\alias{print.wt.filter}
\alias{summary.wt.filter}
\title{
  Wavelet Transform Filter
}
\description{
  Generates a wavelet transform filter.
}
\usage{
wt.filter(filter="la8", modwt=FALSE, level=1)
}
\arguments{
  \item{filter}{A character string indicating which wavelet transform
    filter to compute or a numeric vector of wavelet (high pass) filter
    coefficients (not scaling (low pass) coefficients). If a numeric
    vector is supplied, the length must be even.}
  \item{modwt}{A logical value indicating whether to compute the
    maximal overlap discrete wavelet transform filter.}
  \item{level}{An integer value indicating the level of the wavelet
    filter to compute.}
}
\details{
  The character strings currently supported are derived from one of four
  classes of wavelet transform filters: Daubechies, Least Asymetric,
  Best Localized and Coiflet. The prefixes for filters of these classes
  are \code{d}, \code{la}, \code{bl} and \code{c},
  respectively. Following the prefix, the filter name consists of an
  integer indicating length.  Supported lengths are as follows:
  \describe{
    \item{Daubechies}{2,4,6,8,10,12,14,16,18,20.}
    \item{Least Asymetric}{8,10,12,14,16,18,20.}
    \item{Best Localized}{14,18,20.}
    \item{Coiflet}{6,12,18,24,30.}
    }
  Thus, to obtain the Daubechies wavelet transform filter of length 4,
  the character string \code{"d4"} can be passed to
  \code{wt.filter}.
  \cr
  \cr
  This naming convention has one exception: the Daubechies wavelet
  transform filter of length 2 is denoted by \code{haar} instead of
  \code{d2}.
}
\value{
  Returns an object of class \code{wt.filter}, which is an S4 object
  with slots 
  \item{L}{An integer representing the length of the wavelet and scaling
    filters.}
  \item{h}{A numeric vector of wavelet filter coefficients.}
  \item{g}{A numeric vector of scaling filter coefficients.}
  \item{wt.class}{A character string indicating the class of the wavelet
    transform filter. Possible values are \code{"Daubechies"},
    \code{"Least Asymetric"}, \code{"Best Localized"}, and
    \code{"Coiflet"}. If the \code{wt.filter} object is generated from a
    numeric vector of wavelet coefficients, \code{wt.class} is
    \code{"none"}.}
  \item{wt.name}{A character string indicating the name of the wavlet
    filter as listed in the Details section, above. If the
    \code{wt.filter} object is generated from a numeric vector of
    wavelet coefficients, \code{wt.name} is \code{"none"}.}
  \item{transform}{A character string indicating whether the resulting
    wavelet transform object contains DWT or MODWT coefficients. Possible
    values are \code{"dwt"} and \code{"modwt"}.}
}
\note{
  The notation \code{h} and \code{g} for wavelet and scaling
  coefficients, respectively, follows Percival and Walden (2000).  In
  other texts and articles the reverse notation is often adopted.
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.
}
\seealso{
  \code{\link{wt.filter.qmf}},
  \code{\link{dwt}},
  \code{\link{modwt}}.
}
\examples{
wt.filter("la14")

wt.filter(1:10, modwt=TRUE)
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
