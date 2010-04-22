\name{mra-class}
\Rdversion{1.1}
\docType{class}
\alias{mra-class}
\title{
  Multiresolution Analysis Object
}
\description{
  An S4 object containing multiresolution analysis information.
}
\section{Slots}{
  \describe{
  \item{D}{A list with element \eqn{i} comprised of a matrix containing
    the \eqn{i}th level wavelet detail.}
  \item{S}{A list with element \eqn{i} comprised of a matrix containing
    the \eqn{i}th level wavelet smooths.}
  \item{filter}{A \code{wt.filter} object containing information for
    the filter used in the decomposition. See \code{help(wt.filter)} for
    details.}
  \item{level}{An integer value representing the level of wavelet
    decomposition.}
  \item{boundary}{A character string indicating the boundary method used
    in the wavelet decomposition. Valid values are "periodic" or
    "reflection".}
  \item{series}{The original time series, \code{X}, in matrix format.}
  \item{class.X}{A character string indicating the class of the input
    series.  Possible values are \code{"ts"}, \code{"mts"},
    \code{"numeric"}, \code{"matrix"}, or \code{"data.frame"}.}
  \item{attr.X}{A list containing the attributes information of the
    original time series, \code{X}.  This is useful if \code{X} is an
    object of class \code{ts} or \code{mts} and it is desired to retain
    relevant time information. If the original time series, \code{X}, is
    a matrix or has no attributes, then \code{attr.X} is an empty list.}
  \item{method}{A character string indicating which type of wavelet
    decomposition was performed (either "dwt" or "modwt").}
  }
}
\references{Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods
    for Time Series Analysis}, Cambridge University Press.}
\seealso{
  \code{\link{dwt}},
  \code{\link{modwt}},
  \code{\link{wt.filter}}.
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
\keyword{classes}