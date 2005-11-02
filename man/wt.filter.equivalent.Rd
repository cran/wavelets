\name{wt.filter.equivalent}
\alias{wt.filter.equivalent}
\title{
  Equivalent Wavelet Transform Filter
}
\description{
  Generates an equivalent filter for a specified wavelet transform
  filter and level..
}
\usage{
wt.filter.equivalent(wt.filter, J)
}
\arguments{
  \item{wt.filter}{A \code{wt.filter} object.}
  \item{J}{Level of equivalent filter to compute.}
}
\details{
  A wavelet transform of level \eqn{j} can be computed by iterating the
  pyramid algorithm \eqn{j} times with a level 1 wavelet filter, or by
  simply implementing the pyramid once with the level \eqn{j} equivalent
  filter.  The equivalent filter is obtained by iteratively filtering
  the level 1 wavelet filter and scaling filters with upsampled versions
  of themselves.  For details regarding the computation, see Percival
  and Walden (2000), eqs. 95e and 96e.
}
\value{
  Returns the original \code{wt.filter} object with modifications made
  to the wavelet and scaling filter coefficients, filter length and
  filter level.
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.
}
\seealso{
  \code{\link{wt.filter}}.
}
\examples{
wf <- wt.filter("la14")
wt.filter.equivalent(wf, 3)
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
