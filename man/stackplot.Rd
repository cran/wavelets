\name{stackplot}
\alias{stackplot}
\title{
  Stacked Plotting
}
\description{
  Plot an object as a stacked set of plots, plotting each column using
  the same x-axis and plotting the y-values at the same scale with an
  offset equal to the magnitude of each column vector.
}
\usage{
stackplot(x, y = NULL, type = "h", axes.labels = FALSE,
xlab = "", ylab = "", y.llabs = NULL, y.rlabs = NULL, draw.divides =
TRUE, xtick.at = NULL, ytick.at = NULL, col = "black", main = "")
}
\arguments{
  \item{x}{Specifies the coordinates of the points to be
    plotted. Alternatively, \code{x} can be a list containing matrices
    representing the values to be plotted, or a time series object. The
    x-values will be specified by the rownames of each matrix. The
    y-values, or offsets of each level, will be represented by the
    values of the elements in the columns of each matrix.}
  \item{y}{If \code{x} specifies the coordinates of the points to be
    plotted, then \code{y} can be a list containing matrices that 
    represent the values to be plotted.  The y-values, or offsets of
    each level will be represented by the values of the elements in the
    columns of each matrix.}
  \item{type}{Type of plots to be drawn.  Currently only \code{h}, or
    histogram, and \code{l}, or lines, are implemented.}
  \item{axes.labels}{Boolean representing whether or not to display the
    numeric labels of the horizontal and vertical axes.}
  \item{xlab}{String indicating the label for the horizontal x-axis.}
  \item{ylab}{String indicating the label for the vertical y-axis.}
  \item{y.llabs}{Vector or list of labels to be placed on left side of
    plot.}
  \item{y.rlabs}{Vector or list of labels to be placed on right side of
    plot.}
  \item{draw.divides}{Logical value indicating whether to draw lines to
    separate each subplot in the stackplot.}
  \item{xtick.at}{Vector specifying the tick mark locations on the
    horizontal axis of the plot. Defaulted to par("xaxp").}
  \item{ytick.at}{Vector specifying the tick mark locations on the
    vertical axis of the plot. Defaulted to par("yaxp").}
  \item{col}{String representing color of plot.}
  \item{main}{String representing the title of the plot.}
}
\details{
  Each object specified by \code{x} or \code{y} will be
  coerced into a list, and each object of the list will be coerced into
  a matrix using the function \code{as.matrix}. Each column of a matrix
  will represent the vertical y-values for each individual plot in the
  stacked plot. Horizontal x-values are discussed below.
  
  \code{stackplot} will plot a matrix of values as a series of
  \code{stacked} plots, where the individual columns are plotted using
  the same x and y scales, but offset by the magnitude (= abs(max-min))
  of the column.  

  If the object to be plotted is a Time Series object, the \code{start}
  and \code{end} values of the object will represent the first and last
  horizontal x-values of the points to be plotted. The Time Series
  will be coerced into matrix.  

  The limits of the horizontal x-axis and vertical y-axis will be
  dictated by the data in the first object of the list.

  If the object to be plotted is a list of matrices, and \code{x} is not
  specified, the horizontal x-values of the points to be plotted for a
  given matrix will default to the rownames of that matrix. The rownames
  can be either \code{numeric} or \code{character}. If the first matrix
  of the list does not have rownames, then the horizontal x-values for
  each matrix (and its corresponding plot) will default to range from 1
  to the length of the columns of the respective matrix.

  If the first matrix of the list does have specified rownames, then the
  horizontal x-values assume the values of the rownames of the
  respective matrix. However, if another matrix of the list does not
  have specified rownames, then its horizontal x-values will default to
  range from 1 to the length of the columns of that respective matrix.

  The user can specify the type of plot (currently only "h" and "l") to
  be drawn for every element of the list in \code{x} or \code{y} using
  the vector \code{type}. If the Nth element (N>1) in \code{type} is NA
  and the Nth element of the list in \code{x} or \code{y} exists, then
  the type of plot for the Nth element of the list will be equal to the
  type of plot for the N-1th element of the list. Similarly for
  \code{col}, which specifies the colors for each element of the list in
  \code{x} or \code{y}.
 
  The labels \code{y.llabs} and \code{y.rlabs} will be placed at the
  "zero" lines of each individual plot starting from the bottom working
  up to the top.

}
\author{
    Kelvin Ma, kkym@u.washington.edu
}
\seealso{
  \code{\link{plot.dwt}} and \code{\link{plot.modwt}}
}
\examples{
A <- array(c(1,2,3,4), c(3,4))

# Plotting the vectors c(1,2,3), c(4,1,2), c(3,4,1), and c(2,3,4) in a
# stacked plot with 4 separate plots).  The x-values will be 1:3.
stackplot(A)

#Plotting A with numeric labels on the axis and making the plot "red".
#Also labeling the x-axis with "X-Label", and the y-axis with "Y-Label".
stackplot(A, axes.labels = TRUE, xlab = "X-Label", ylab = "Y-Label")

timeSeries <- ts(A)
# Plotting the Time Series created by coercing the array A into a Time
# Series.  Observe that this plot similar to the plot in the first example.
stackplot(timeSeries)

lfig <- 2:11
hfig <- 1:10
# Plotting a stacked plot of two types: histogram and lines.  Where the
# line joins the points (1,2) and (10,11), and the histogram is similar to
# the histogram formed by entering plot(1:10, type = "h").  The line will
# be red and the histogram will be yellow.
stackplot(list(lfig, hfig), type = c("l", "h"), col = c("red", "yellow"))

}
\keyword{ts}
