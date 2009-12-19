\name{getAllInt}
\alias{getAllInt}

\title{ Get the summerized MAS5 values for a given spike probe }
\description{
    This function retrieves the expression intensities prvided by the object of
    class \code{"YAQCStats"} for the probe which names are compatible with the 
    given pattern and return their mean value.
}
\usage{
  getAllInt(YAQCStatsObject,pattern)
}

\arguments{
   \item{YAQCStatsObject}{ an object of type \code{"YAQCStats"} }
   \item{pattern}{ a pattern used to select the probe names to be used }
}

\value{
    An object of type \code{"numeric"} in which the mean expression intensities 
    of each array are given.
}


\examples{
    ## load a dataset
    library(affydata)
    data(Dilution)
    ## perform quality control
    qc <- yaqc(Dilution)
    ## get intensities for the biob
    ## spikes probe sets
    getAllInt(qc,"biob")
    ## or
    getAllInt(qc,"b[3|5|m]")
}


\author{ Laurent Gatto }

\keyword{ misc }
