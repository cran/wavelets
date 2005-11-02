\name{scalingshift.dwt}
\alias{scalingshift.dwt}
\title{
  Find Circular Shift of DWT Scaling Coefficients
}
\description{
  Finds circular shift of DWT scaling coefficients at a specified level,
  for a specified filter length.
}
\usage{
scalingshift.dwt(L, j, N = NULL)
}
\arguments{
  \item{L}{Length of wavelet transform filter used.}
  \item{j}{Level of DWT scaling coefficients to be shifted.}
  \item{N}{Length of original time series sample.}
}
\details{
  This function computes the circular shift associated with a wavelet
  transform filter of length \code{L} and level \code{j}. \code{L} must
  be of even length.

  Each coefficient's index will be circularly shifted forward by the
  value outputted. For example, if we have a vector of DWT scaling
  coefficients (with NAs inserted) with indices 1,2,3,4,5,6,7,8 and the
  value outputted from the function \code{scalingshift.dwt} is 2, the
  new order of the vector of aligned DWT wavelet coefficients would be
  7,8,1,2,3,4,5,6.

  If \code{N} is specified, the function will output the "shift" Modulo
  N. Because the shift is circular, the value outputted when N is
  specified is equivalent to the case when N is not specified.

  This function calculates the absolute value of vjH in formulas (114c)
  in \emph{Wavelet Methods for Time Series Analysis} by Percival and
  Walden (2000). (And optionally vjH Modulo N)
}
\value{
  \item{shift}{Circular shift for specified level of DWT scaling
    coefficients for a specified filter length.}
}
\author{
    Kelvin Ma, kkym@u.washington.edu
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
  Series Analysis}, Cambridge University Press.  
}
\seealso{
  \code{\link{waveletshift.dwt}} and \code{\link{plot.dwt}}
}
\examples{
# Finding the circular shift for a wavelet transform filter of length 8,
# pertaining to DWT Scaling Coefficients of level 5.
scalingshift.dwt(8, 5)

# If the sample size of the original time series is of length 1024.
scalingshift.dwt(8, 5, N = 1024)
}
\keyword{ts}
