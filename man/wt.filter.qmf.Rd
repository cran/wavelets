\name{wt.filter.qmf}
\alias{wt.filter.qmf}
\title{
  Quadrature Mirror Filter
}
\description{
  Computes the quadrature mirror filter of a series of even length.
}
\usage{
wt.filter.qmf(x, inverse=FALSE)
}
\arguments{
  \item{x}{A numeric vector of even length}
  \item{inverse}{A logical flag indicating whether to compute the
    inverse quadrature mirror filter.}
}
\details{The quadrature mirror filter is computed as outlined in
  Percival and Walden, page 75.  Specifically, the forward QMF makes use
  of a wavelet filter to obtain a scaling filter and the inverse QMF
  does the opposite.}
\value{
  A numeric vector which contains the values of \code{x} in reverse
  order, with the values of the resulting vector at even indices being
  negated in the case that \code{inverse=FALSE} and the values at odd
  indices being negated in the case that \code{inverse=TRUE}.
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
\examples{
# obtain both the la8 wavelet and scaling filters
la8 <- wt.filter("la8")

# now compare output of qmf with values in wt.filter object
wt.filter.qmf(la8@h)
la8@g
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
