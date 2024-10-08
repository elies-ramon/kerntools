% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_study.R
\name{histK}
\alias{histK}
\title{Kernel matrix histogram}
\usage{
histK(K, main = "Histogram of K", vn = FALSE, ...)
}
\arguments{
\item{K}{Kernel matrix (class "matrix").}

\item{main}{Plot title.}

\item{vn}{If TRUE, the value of the von Neumann entropy is shown in the plot.
(Defaults: FALSE).}

\item{...}{further arguments and graphical parameters passed to `plot.histogram`.}
}
\value{
An object of class "histogram".
}
\description{
`histK()` plots the histogram of a kernel matrix.
}
\details{
Information about the von Neumann entropy can be found at '?vonNeumann()'.
}
\examples{
data <- matrix(rnorm(150),ncol=50,nrow=30)
K <- RBF(data,g=0.01)
histK(K)
}
