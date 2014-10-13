\name{getRatioProbes}
\alias{getRatioProbes}

\title{ Get the names of degradation control probes on the array }
\description{
  This function returns the probes names used for degradation control
  that are located on the given GeneChip.
}
\usage{
  getRatioProbes(object,onlyFirst)
}

\arguments{
  \item{object}{ An object of class \code{"AffyBatch"} or
    \code{"ExpressionSet"} }
  \item{onlyFirst}{ Boolean defining of only first or all instances
    found should be returned. Default is set to TRUE. Warnings are
    returned if more than one probe is found. The function stops with
    an error if no probe is found.}
}

\value{
  An object of type \code{"character"} with all the Affymetrix
  degradation control probe names. 
}

\examples{
    library(yaqcaffy)
    ## load a dataset
    library(affydata)
    data(Dilution)
    getRatioProbes(Dilution)
}

\seealso{getSpikeProbes,getBioProbes}

\author{ Laurent Gatto }

\keyword{ misc }
