\name{align}
\alias{align}
\title{
  Align Wavelet Transform Coefficients
}
\description{
  Aligns both wavelet (high pass) and scaling (low pass) coefficients of
  objects of class \code{dwt}, \code{modwt}, \code{dwpt}, and
  \code{modwpt} using phase shift values computed by
  \code{wt.filter.shift}.
}
\usage{
align(wt, coe=FALSE, inverse=FALSE)
}
\arguments{
  \item{wt}{An object of class \code{dwt}, \code{modwt}, \code{dwpt}, or
    \code{modwpt}.}
  \item{coe}{Logical value indicating whether to use center of energy
    method in computing phase shifts.}
  \item{inverse}{Logical value indicating whether to shift wavelet and
    scaling coefficients of an aligned object back to their original
    positions.}
}
\details{
  An object of class \code{dwt}, \code{modwt}, is characterized as
  'unaligned' if the value in the \code{aligned} slot is \code{FALSE}.
  Similarly, these objects are classified as 'aligned' if the value in
  the \code{aligned} slot is \code{TRUE}. Thus, \code{align} will only
  operate on an 'unaligned' wavelet transform object if \code{inverse =
    FALSE} and on an 'aligned' wavelet transform object if \code{inverse
    = TRUE}.

  The argument \code{coe} is passed to the \code{wt.filter.shift}
  function to determine what method to use for computing phase shifts
  (see documentation for \code{wt.filter.shift}).
}
\value{
  Either an aligned or unaligned object of the same class as \code{wt}
  (see Details above). 
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
    Series Analysis}, Cambridge University Press.
}
\seealso{
  \code{\link{wt.filter.shift}},
  \code{\link{dwt}},
  \code{\link{modwt}}.
}
\examples{
# obtain the two series listed in Percival and Walden (2000), page 42
X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)

# combine them and compute DWT
newX <- cbind(X1,X2)
wt <- dwt(newX, n.levels=3, boundary="reflection", fast=FALSE)

# align
wt.aligned <- align(wt)
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
