\name{dwt-class}
\Rdversion{1.1}
\docType{class}
\alias{dwt-class}
\title{
  Discrete Wavelet Transform Object
}
\description{
  An S4 object containing discrete wavelet transform information.
}
\section{Slots}{
  \describe{
  \item{W}{A list with element \eqn{i} comprised of a matrix containing
    the \eqn{i}th level wavelet coefficients.}
  \item{V}{A list with element \eqn{i} comprised of a matrix containing
    the \eqn{i}th level scaling coefficients.}
  \item{filter}{A \code{wt.filter} object containing information for
    the filter used in the decomposition. See \code{help(wt.filter)} for
    details.}
  \item{level}{An integer value representing the level of wavelet
    decomposition.}
  \item{n.boundary}{A numeric vector indicating the number of boundary
    coefficients at each level of the decomposition.}
  \item{boundary}{A character string indicating the boundary method used
    in the decomposition. Valid values are "periodic" or "reflection".}
  \item{series}{The original time series, \code{X}, in matrix format.}
  \item{class.X}{A character string indicating the class of the input
    series.  Possible values are \code{"ts"}, \code{"mts"},
    \code{"numeric"}, \code{"matrix"}, or \code{"data.frame"}.}
  \item{attr.X}{A list containing the attributes information of the
    original time series, \code{X}.  This is useful if \code{X} is an
    object of class \code{ts} or \code{mts} and it is desired to retain
    relevant time information. If the original time series, \code{X}, is
    a matrix or has no attributes, then \code{attr.X} is an empty list.}
  \item{aligned}{A logical value indicating whether the wavelet and
    scaling coefficients have been phase shifted so as to be aligned
    with relevant time information from the original series. The value
    of this slot is initially FALSE and can only be changed to TRUE via
    the \code{align} function, with the \code{dwt} object as input.}
  \item{coe}{A logical value indicating whether the center of energy
    method was used in phase alignement of the wavelet and scaling
    coefficients. By default, this value is FALSE (and will always be
    FALSE when \code{aligned} is FALSE) and will be set to true if the
    \code{dwt} object is phase shifted via the \code{align} function and
    center of energy method.}
  }
}
\references{Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods
    for Time Series Analysis}, Cambridge University Press.}
\seealso{
  \code{\link{dwt}},
  \code{\link{modwt}},
  \code{\link{modwt-class}},
  \code{\link{wt.filter}}.
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
\keyword{classes}