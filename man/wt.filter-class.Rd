\name{wt.filter-class}
\Rdversion{1.1}
\docType{class}
\alias{wt.filter-class}
\title{
  Wavelet Transform Filter Object
}
\description{
  An S4 object containing wavelet transform filter information.
}
\section{Slots}{
  \describe{
  \item{L}{An integer representing the length of the wavelet and scaling
    filters.}
  \item{level}{An integer representing the level of the wavelet filter.}
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
}
\note{
  The notation \code{h} and \code{g} for wavelet and scaling
  coefficients, respectively, follows Percival and Walden (2000).  In
  other texts and articles the reverse notation is often adopted.
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods
    for Time Series Analysis}, Cambridge University Press.
}
\seealso{
  \code{\link{wt.filter.qmf}},
  \code{\link{dwt}},
  \code{\link{modwt}}.
}

\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
\keyword{classes}