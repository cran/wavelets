\name{dwt.forward}
\alias{dwt.forward}
\alias{dwt.backward}
\alias{modwt.forward}
\alias{modwt.backward}
\title{
  Discrete Wavelet Transform and Maximal Overlap Discrete Wavelet
  Tranform Forward and Backward Pyramid Algorithm
}
\description{
  Implementation of DWT and MODWT forward and backward pyramid
  algorithms.
}
\usage{
dwt.forward(V, filter)
dwt.backward(W, V, filter)
modwt.forward(V, filter, j)
modwt.backward(W, V, filter, j)
}
\arguments{
  \item{W}{A vector of wavelet coefficients.}
  \item{V}{A vector of scaling coefficients.}
  \item{filter}{A \code{wt.filter} object.}
  \item{j}{The level of wavelet and scaling coefficients to compute (for
    forward algorithm) or the level of wavelet and scaling coefficient
    inputs (for inverse algorithm).}
}
\details{An implementation of the DWT and MODWT forward and backward
  pyramid algorithms using pseudocode written by Percival and Walden 
  (2000), pp. 100-101, 177-178. These functions are intended primarily
  as helper functions for the \code{dwt}, \code{modwt}, \code{idwt} and
  \code{imodwt} functions.}
\value{
  \code{dwt.forward} and \code{modwt.forward} return a list of two
  elements containing vectors of wavelet and scaling coefficients for
  the subsequent level of analysis. \code{dwt.backward} and
  \code{modwt.backward} return a vector of scaling coefficients for the
  previous level of analysis.
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

# compute the LA8 wavelet filter for both DWT and MODWT
la8.dwt <- wt.filter()
la8.modwt <- wt.filter(modwt=TRUE)

# compute the DWT and MODWT level one wavelet and scaling coefficients
wt.dwt <- dwt.forward(X1, la8.dwt)
wt.modwt <- modwt.forward(X2, la8.modwt, 1)

# compute the original series with the level one coefficients
newX.dwt <- dwt.backward(wt.dwt$W, wt.dwt$V, la8.dwt)
newX.modwt <- modwt.backward(wt.modwt$W, wt.modwt$V, la8.modwt, 1)
}
\author{Eric Aldrich. ealdrich@gmail.com.}
\keyword{ts}
