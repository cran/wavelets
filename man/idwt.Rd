\name{idwt}
\alias{idwt}
\title{
  Inverse Discrete Wavelet Transform
}
\description{
  Computes the inverse discrete wavelet transform for a discrete wavelet
  transform that was obtained from a univariate or multivariate time
  series.
}
\usage{
idwt(wt, fast=TRUE)
}
\arguments{
  \item{wt}{A \code{dwt} object.}
  \item{fast}{A logical flag which, if true, indicates that the inverse
    pyramid algorithm is computed with an internal C function.
    Otherwise, only R code is used in all computations.}
}
\details{The inverse discrete wavelet transform is computed via the
  inverse pyramid algorithm, using pseudocode written by Percival and
  Walden (2000), p. 101.}
\value{An object with class and attributes equivalent to the original
  series that was used to compute the DWT.  In general, the output will
  be equivalent to the original series (i.e. X = idwt(dwt(X))), however
  when thresholding or shrinkage methods are used on the \code{dwt}
  object, the output of \code{idwt} may differ from the original series.
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.}
\seealso{
  \code{\link{dwt}},
  \code{\link{modwt}},
  \code{\link{imodwt}}.
}
\examples{
# obtain the two series listed in Percival and Walden (2000), page 42
X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)

# combine them and compute DWT
newX <- cbind(X1,X2)
wt <- dwt(newX, n.levels=3, boundary="reflection")

# compute the inverse DWT
invX <- idwt(wt)

# compare
newX
invX
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
