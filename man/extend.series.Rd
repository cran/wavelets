\name{extend.series}
\alias{extend.series}
\title{
  Extend a Time Series
}
\description{
  Extends a univariate or multivariate time series beyond its original
  length.
}
\usage{
extend.series(X, method="reflection", length="double", n, j)
}
\arguments{
  \item{X}{A univariate or multivariate time series. Numeric vectors,
    matrices and data frames are accepted.}
  \item{method}{A character string indicating which extension method to
    use. Possible values are \code{"periodic"}, \code{"reflection"},
    \code{"zero"}, \code{"mean"}, and \code{"reflection.inverse"}.}
  \item{length}{A character string indicating how to determine the
    length of the extended series. Possible values are
    \code{"aribitrary"}, \code{"powerof2"}, and \code{"double"}.}
  \item{n}{An integer value specifying the length of the extended
    series. This argument is only relevant when \code{length =
      "arbitrary"}.}
  \item{j}{An integer value specifying the power of two of which the
    length of the extended series should be a multiple. This argument is
    only relevant for \code{length = "powerof2"}.}
}
\details{The original time series, \code{X}, is extended to a new length
  as determined by \code{length}: if \code{length = "arbitrary"}, a value
  for the argument \code{n} must be specified in order to extend the
  series to length \code{n}; if \code{length = "powerof2"} the series
  will be extended to the nearst multiple of \eqn{2^j}, where \code{j}
  must be specified as an argument to the function;  if \code{length =
  "double"} the series will be extended to twice its original length.

  Once the length of the extended series is determined, the actual
  extension is performed using one of five methods: if \code{method =
  "periodic"} a periodic extension is made by concatenating the series
  with itself (without changing the order of the values); if
  \code{method = "reflection"} an extension is made by concatenating
  the series with a reflection (reverse ordering) of itself; if
  \code{method = "zero"} the series is padded with zeros; if
  \code{method = "mean"} the series is padded with its mean; if
  \code{method = "reflection.inverse"} an extension is made by
  concatenating the series with a reflection (reverse ordering) of
  itself that is further reflected of the horizontal axis \eqn{y = x_N}
  where \eqn{x_N} is the last value of the observed series.

  In the case that the series is extended beyond twice it's original
  length, the \code{"periodic"} and \code{"reflection"} methods are
  repetetively applied so that every subseries of length \eqn{N} (the
  length of the original series) is either a periodic or reflected
  extension of any other \eqn{N} adjacent values. Similarly, the
  \code{"reflection.inverse"} method is repetetively applied so that
  every \eqn{2N} values are a periodic extension of any adjacent
  \eqn{2N} values.
}
\value{
  The extended time series, in its original class.
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.}
\seealso{
  \code{\link{dwt}},
  \code{\link{modwt}}.
}
\examples{
# obtain the two series listed in Percival and Walden (2000), page 42
X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)

# combine them and extend using reflection method
newX <- cbind(X1,X2)
newX.ext1 <- extend.series(newX, length="arbitrary", n=64)
plot.ts(newX.ext1)

# apply another extension method
newX.ext2 <- extend.series(newX, method="reflection.inverse",
                           length="powerof2", j=6)
plot.ts(newX.ext2) 
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
