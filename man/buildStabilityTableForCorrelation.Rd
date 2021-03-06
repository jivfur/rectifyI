% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rectifiedI.R
\name{buildStabilityTableForCorrelation}
\alias{buildStabilityTableForCorrelation}
\title{Finds how many iterations are necessary to achieve stability in resampling method for rectifying I through pearson corrrelation.}
\usage{
buildStabilityTableForCorrelation(data, times = 10, samples = 100,
  plots = TRUE)
}
\arguments{
\item{data}{data structure after loading the file using \code{loadFile} function}

\item{times}{the number of times \code{rescaleI} will be executed. The default value is 100.}

\item{samples}{size of the resampling method. The default value is 1000}

\item{plots}{to draw the significance plot}
}
\value{
A vector with the average \eqn{\log(samples)} averages I
}
\description{
\code{buildStabilityTableForCorrelation}  finds how many iterations are necessary to achieve stability in resampling method, plotting in a log scale.
}
\examples{
fileInput <- system.file("testdata", "chen.csv", package="Irescale")
data <- loadFile(fileInput)
resultsChen<-buildStabilityTableForCorrelation(data=data,times=10,samples=100,plots=TRUE)
}
