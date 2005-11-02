\name{mra}
\alias{mra}
\title{
  Multiresolution Analysis
}
\description{
  Computes the multiresolution analysis for a univariate or multivariate
  time series.
}
\usage{
mra(X, filter="la8", n.levels, boundary="periodic", fast=TRUE, method="dwt")
}
\arguments{
  \item{X}{A univariate or multivariate time series. Numeric vectors,
    matrices and data frames are also accepted.}
  \item{filter}{Either a \code{wt.filter} object, a character string
    indicating which wavelet filter to use in the decomposition, or a
    numeric vector of wavelet coefficients (not scaling coefficients). See
    \code{help(wt.filter)} for acceptable filter names.}
  \item{n.levels}{An integer specifying the level of the decomposition.
    By default this is the value J such that the length of \eqn{X} is at
    least as great as the length of the level \eqn{J} wavelet filter,
    but less than the length of the level \eqn{J+1} wavelet
    filter. Thus, \eqn{J \le \log{(\frac{N-1}{L-1}+1)}}{j <=
      log((N-1)/(L-1)+1)}, where \eqn{N} is the length of \eqn{X}.}
  \item{boundary}{A character string indicating which boundary method to
    use. \code{boundary = "periodic"} and \code{boundary = "reflection"}
    are the only supported methods at this time.}
  \item{fast}{A logical flag which, if true, indicates that the pyramid
    algorithm is computed with an internal C function.  Otherwise, only
    R code is used in all computations.}
  \item{method}{A character string, taking values "dwt" or "modwt",
    that indicates which type of transform to use when computing the
    MRA.}
}
\value{
  Returns an object of class \code{mra}, which is an S4 object with
  slots 
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
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.}
\seealso{
  \code{\link{dwt}},
  \code{\link{modwt}},
  \code{\link{wt.filter}}.
}
\examples{
# obtain the two series listed in Percival and Walden (2000), page 42
X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)

# combine them and compute MRA
newX <- cbind(X1,X2)
mra.out <- mra(newX, n.levels=3, boundary="reflection")
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
