\name{plot.modwt}
\alias{plot.modwt}
\title{
  Plot MODWT Coefficients
}
\description{
  Plot MODWT wavelet and scaling coefficients.
}
\usage{
\method{plot}{modwt}(x, levels = NULL, draw.boundary = FALSE, type = "stack",
col.plot = "black", col.boundary = "red", X.xtick.at = NULL, X.ytick.at
= NULL, Stack.xtick.at = NULL, Stack.ytick.at = NULL, X.xlab = "t",
y.rlabs = TRUE, plot.X = TRUE, plot.W = TRUE, plot.V = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{dwt}.}
  \item{levels}{Number, vector, or list of two vectors indicating range
    of levels to plot.  See details.}
  \item{draw.boundary}{Logical value indicating whether to draw boundary
    coefficients.}
  \item{type}{Type of plot to draw.  Currently only "stack" is
    implemented.}
  \item{col.plot}{Color of wavelet and scaling coefficients.}
  \item{col.boundary}{Color of boundary coefficient lines.}
  \item{X.xtick.at}{Vector specifying the extreme tick mark locations
    and the number of intervals between those extreme tick marks on the
    horizontal axis of the plot of the original times series. This
    vector takes the form similar to par("xaxp"), and is defaulted to
    par("xaxp") when \code{X.xtick.at} is not specified.}
  \item{X.ytick.at}{Vector specifying the extreme tick mark locations
    and the number of intervals between those extreme tick marks on the
    vertical axis of the plot of the original times series. This vector
    takes the form similar to par("yaxp"), and is defaulted to
    par("yaxp") when \code{X.ytick.at} is not specified.}
  \item{Stack.xtick.at}{Vector of form similar to \code{X.xtick.at}
    specifying the tick mark locations on the horizontal
    axis of the stacked plot. This is applicable only if \code{plot.X}
    is FALSE.}
  \item{Stack.ytick.at}{Vector of form similar to \code{X.ytick.at}
    specifying the tick mark locations on the vertical
    axis of the stacked plot. This is applicable only if \code{plot.X}
    is FALSE.}
  \item{X.xlab}{String specifying the label of the horizontal axis of
    the plot of the original time series.}
  \item{y.rlabs}{Logical flag indicating whether to draw the vertical
    labels on the right vertical axis of the stacked plot.  These labels
    indicate the number of boundary coefficients to the right of the
    right boundary line, and the level of decomposition of the wavelet
    (or scaling) coefficients.}
  \item{plot.X}{Logical flag indicating whether to draw original time
    series.}
  \item{plot.W}{Logical flag indicating whether to draw the wavelet
    coefficients in the stacked plot.}
  \item{plot.V}{Logical flag indicating whether to draw the scaling
    coefficients in the stacked plot.}
  \item{...}{Additional paramters that are acceptable arguments to the
    generic \code{plot} function}
}
\details{
  \code{plot.modwt} plots the MODWT wavelet and scaling coefficients.
  The \code{type} parameter species the type of plot, which is currently
  only stack. If the MODWT object is defined for multiple time series,
  only the data pertaining to the first time series of the MODWT object
  is plotted. Thus, only the wavelet coefficients and scaling
  coefficients of the first time series of the MODWT object will be
  plotted. If a single number is specified for \code{levels}, then the
  wavelet coefficients of levels 1 through \code{levels} will be
  plotted. Otherwise, a vector or the first element of a list will
  specify which levels of the wavelet coefficients will be
  plotted. Unless specified in the second element of a list, only one
  level of scaling coefficients will be plotted and this level is equal
  to the highest level of the wavelet coefficients plotted.

  For each respective axis, the distance between a tick mark on the plot
  of the original time series is equivalent to the distance between a
  tick mark on the stacked plot. Thus, when altering the relative
  spacing of the tick marks on the plot of the original time series
  using \code{X.xtick.at} or \code{X.ytick.at}, the tick marks of the
  stacked plot are automatically adjusted. If the plot of the original
  time series is not drawn, then the user can alter the spacing of the
  tick marks using \code{Stack.xtick.at} and \code{Stack.ytick.at} for
  the horizontal and vertical axes, respectively.
  
  One of \code{plot.W} or \code{plot.V} must be TRUE.
}
\author{
    Kelvin Ma, kkym@u.washington.edu
}
\references{
  Percival, D. B. and A. T. Walden (2000) \emph{Wavelet Methods for Time
  Series Analysis}, Cambridge University Press.  
}
\seealso{
  \code{\link{stackplot}}
}
\examples{
X <- rnorm(2048)
modwtobj <- modwt(X)

# Plotting wavelet coefficients at levels 1 through 6 and scaling
# coefficients at level 6.
plot.modwt(modwtobj, levels = 6)

# Plotting wavelet coefficients of levels 1, 3, 5, and scaling
# coefficients of levels 4 and 5, and green boundary coefficients
# lines.
plot.modwt(modwtobj, levels = list(c(1,3,5),c(4,5)), draw.boundary =
TRUE, col.boundary = "green")

# Plotting wavelet coefficients at level 1 through 6 and not plotting
# any scaling coefficients.
plot.modwt(modwtobj, levels = 6, plot.V = FALSE)
}
\keyword{ts}
