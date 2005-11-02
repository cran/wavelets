\name{wt.filter.shift}
\alias{wt.filter.shift}
\title{
  Wavelet Transform Filter Phase Shift
}
\description{
  Computes phase shifts for wavelet transform coefficients corresponding
  to a particular filter.
}
\usage{
wt.filter.shift(filter, J, wavelet=TRUE, coe=FALSE, modwt=FALSE)
}
\arguments{
  \item{filter}{A \code{wt.filter} object, a character string indicating
    which wavelet transform filter to compute, or a numeric vector of
    wavelet (high pass) filter coefficients (not scaling (low pass)
    coefficients). If a numeric vector is supplied, the length must be
    even.}
  \item{J}{A vector of positive integers indicating levels for which
    to compute phase shifts.}
  \item{wavelet}{A logical flag indicating whether to compute the
    wavelet (high pass) or scaling (low pass) phase shift(s).}
  \item{coe}{A logical value indicating whether to use the center of
    energy method (see Percival and Walden 2000, page 118) to compute
    the phase shift(s).}
  \item{modwt}{A logical value indicating whether to compute MODWT phase
    shift(s).}
}
\details{
  For wavelet filters of class 'Least Assymetric' or 'Coiflet', the
  default method for calculating phase shifts is outlined on pages
  112-114 and page 124 of Percival and Walden 2000.  For wavelet filters of
  class 'Best Localized', the default shifts are given on page 119 of
  Percival and Walden 2000. For the Haar filter, both the level \eqn{j}
  wavelet and scaling phase shifts are half the length of the level
  \eqn{j} wavelet and scaling filters and the phase shifts for the D(4)
  filter are determined by specifying \eqn{\nu = -1} in equations (114a)
  and (114b) of Percival and Walden 2000.

  For all other filters, the center of energy method is the default
  method for computing phase shifts (thus rendering the \code{coe}
  argument irrelevant). If \code{coe = TRUE}, then the center of energy
  method is used regardless of filter class.

  By default the DWT phase shifts are computed, using the MODWT phase
  shifts and the methods outlined on pages 115-116 of Percival and
  Walden 2000. If \code{modwt = TRUE}, then only the MODWT phase shifts
  are computed.
}
\value{
  A vector of phase shifts, expressed in absolute value.  Each element
  corresponds to the equivalent element in \code{J}.
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.
}
\seealso{
  \code{\link{wt.filter}},
  \code{\link{align}},
}
\examples{
wt.filter.shift("la14", J=1:6)
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
