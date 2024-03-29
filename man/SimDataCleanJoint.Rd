% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimDataCleanJoint.R
\name{SimDataCleanJoint}
\alias{SimDataCleanJoint}
\title{SimDataCleanJoint normalizes and standardizes both data sets. By default it standardizes to the joint max signal.}
\usage{
SimDataCleanJoint(
  filename1,
  filename2,
  kmin = 2,
  rel = "Joint",
  normvector = list("None", "None"),
  logoption = TRUE,
  kmin_sub = list(0, 0),
  rel_force = F,
  MVCorrection = T
)
}
\arguments{
\item{filename1}{First file or dataframe: sample 1}

\item{filename2}{Second file or dataframe: sample 2}

\item{kmin}{Minimum number of identifications needed to be considered real. Default=2}

\item{rel}{Choice of Standardization method. 'Within' scales compared to largest signal size in each sample. 'Joint' scales to largest signal value between both files. 'AsIs' does no scaling. NUMERIC scales to a specified number. Default='Joint'}

\item{normvector}{List of Normalization vectors to be multiplied against. Vector length should equal sample size. Default=list('None','None')}

\item{logoption}{Boolean indicating use of log transformation. Default=TRUE}

\item{kmin_sub}{Minimum number of identifications needed to be seen in a given file. Default=list(0,0)}

\item{rel_force}{Forces relativization even if normalization vector is given. Default=F}

\item{MVCorrection}{default=TRUE. Disable for elimination of Presence Terms in multiplicaiton}
}
\value{
Normalized and standardized cleaned data
}
\description{
SimDataCleanJoint normalizes and standardizes both data sets. By default it standardizes to the joint max signal.
}
\examples{
#
}
