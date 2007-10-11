\name{getOutliers}
\alias{getOutliers}

\title{ Get outliers for the different YAQCStatsObject slots }
\description{
    This function retrives the outliers fot the different quality control metrics
    stored in a YAQCStatsObject. Outliers are defined as being outside of the 
    mean +/- 2 stdev range or mean/2, mean*1.5 for the scale factor.
}
\usage{
  getOutliers(YAQCStatsObject,slot)
}

\arguments{
   \item{YAQCStatsObject}{ an object of type \code{"YAQCStats"} }
   \item{slot}{ an object of type \code{string} describing the slot for which 
                the outliers should be retrieved (see details for possible 
                slot strings)
              }
}

\value{
    An object of type \code{"numeric"} giving the outliers names and values
}

\details{
   The slot strings that can be used are:
    \begin{description}
     \item[scale factor] "sfs"
     \item[average background] "avbg"
     \item[average noise] "avns"
     \item[percentage present] "pp"
     \item[$\beta$-actin 3'/5' ratio] "actin"
     \item[GAPDH 3'/5' ratio] "gapdh"
     \item[internal bioB control] "biob"
     \item[internal bioC control] "bioc"
     \item[internal bioD control] "biod"
     \item[Dap spike control] "dap"
     \item[Thr spike control] "thr"
     \item[Phe spike control] "phe"
     \item[Lys spike control] "lys"
     \end{description}
}

\examples{
  ## load data
  library(MAQCsubsetAFX)
  data(refA)
  ## create the yaqc object
  qobj <- yaqc(refA[,3:5])
  ## get outliers for the scale factor
  getOutliers(qobj,"sfs")
}

\author{ Laurent Gatto }

\keyword{ misc }

