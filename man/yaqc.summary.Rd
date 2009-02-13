\name{yaqc.summary}
\alias{yaqc.summary}
\alias{summary,YAQCStats}
\alias{summary,YAQCStats-method}
\alias{summary,YAQCStats,missing-method}


\title{Summarizes YAQCStats object outliers}
\description{
  Generates a data frame of QC outliers.
}
\usage{yaqc.summary(YAQCStatsObject,latex)}
\arguments{
    \item{YAQCStatsObject}{an object of class \code{"YAQCStats"}.}
    \item{latex}{logical value. If 'TRUE', the data frame is returned as a latex table.}

}
\details{
  The outliers of the \code{YAQCStatsObject} quality control metrics are
  summerized and returned as a data frame. The internal probe
  (bioB, bioC and bioD) and poly-A controls
  (dap, thr, phe and lys) are
  summerized in two rows.

  If \code{latex} is set to 'TRUE', the data frame is returned as
  a latex table (requires the \code{xtable} library).
  
  Note: not tested if the YAQCStatsObject has been generated with
  an Expression Set objects.
}

\examples{
  ## load data
  library(affydata)
  data(Dilution)
  ## create the yaqc object
  ## and plot it
  qobj <- yaqc(Dilution)
  summary(qobj)
}

\author{Laurent Gatto}
 
\keyword{misc}
