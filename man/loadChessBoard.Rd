% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rectifiedI.R
\name{loadChessBoard}
\alias{loadChessBoard}
\title{Loads a chessboard or matrix alike input file.}
\usage{
loadChessBoard(fileName)
}
\arguments{
\item{fileName}{the path and file's name to load.}
}
\value{
data frame with two variables, the first variable is a vector with coordinate x (latitude) and y (longitude), the second variable contains the values of the variable of interest.
}
\description{
\code{loadChessBoard} is used when the input file has a 2D shape, this is a board shape, and it is only one variable of interest.
For example:
\tabular{cccccc}{
1\tab 1\tab 1\tab 1\tab 1\tab 1 \cr
2\tab 2\tab 2\tab 2\tab 2\tab \cr
3\tab 3\tab 3\tab 3\tab 3\tab \cr
4\tab 4\tab 4\tab 4\tab 4\tab \cr
}
}
\examples{
fileInput <- system.file("testdata", "chessboard.csv", package="Irescale")
data<-loadChessBoard(fileInput)
}
