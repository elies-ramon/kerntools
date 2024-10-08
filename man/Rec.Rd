% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performance_measures.R
\name{Rec}
\alias{Rec}
\title{Recall or Sensitivity or TPR}
\usage{
Rec(ct, multi.class = "macro")
}
\arguments{
\item{ct}{Confusion Matrix.}

\item{multi.class}{Should the results of each class be aggregated, and how?
Options: "none", "macro", "micro". (Defaults: "macro").}
}
\value{
TPR (a single value).
}
\description{
`Rec()` computes the Recall, also known as Sensitivity or TPR (True Positive Rate),
 between the output of a classification model and the actual values of the target.
}
\examples{
y <- c(rep("a",3),rep("b",2))
y_pred <- c(rep("a",2),rep("b",3))
ct <- table(y,y_pred)
Rec(ct)
}
