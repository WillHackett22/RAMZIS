% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestSimilarityFunction.R
\name{TestSimilarityFunction}
\alias{TestSimilarityFunction}
\title{TestSimilarityFunction}
\usage{
TestSimilarityFunction(
  df1,
  combo1,
  gl1,
  df2,
  combo2,
  gl2,
  MVCorrection,
  mn = FALSE
)
}
\arguments{
\item{df1}{Dataset 1}

\item{combo1}{Sampling of dataset 1}

\item{gl1}{Glycopeptide list ofdataset 1}

\item{df2}{Dataset 2}

\item{combo2}{Sampling of dataset 2}

\item{gl2}{glycopeptide list of dataset 2}

\item{MVCorrection}{If False, it will ignore missing values rather than count towards the overall average. Default=TRUE}

\item{mn}{The scaling factor. Default will use 1+mean(presence) of a glycopeptide/identification. Setting this to a number will override that process.}
}
\value{
Similarity between two matrices
}
\description{
TestSimilarityFunction
}
\examples{
#
}
