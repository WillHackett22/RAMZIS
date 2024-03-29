% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OverlapCaclulator.R
\name{OverlapCalculator}
\alias{OverlapCalculator}
\title{Overlap Calculator}
\usage{
OverlapCalculator(
  ReferenceDis,
  TestDis,
  zerohandling = T,
  Contribution = F,
  lb = -0.1,
  ub = 1.1,
  na.rm = T
)
}
\arguments{
\item{ReferenceDis}{reference distribution, generally higher}

\item{TestDis}{comparison distribution}

\item{zerohandling}{Default=T, if T, when the Test max is zero it splits the distribution across the two points bordering. As normal, all values between -0.001 and 0.001 are 1}

\item{Contribution}{Default=F, Make T when doing contribution overlaps for higher resolution}

\item{na.rm}{Default=T. Removes NA values from ReferenceDis and TestDis}
}
\value{
Overlap statistics
}
\description{
Overlap Calculator
}
\examples{
#
}
