\name{getBioProbes}
\alias{getBioProbes}

\title{ Get the names of the Bio spike probes on the array }
\description{
  This function returns all the AFFX-Bio probes names that
  are located on the given GeneChip.
}
\usage{
  getBioProbes(object,onlyFirst)
}

\arguments{
  \item{object}{ An object of type \code{"AffyBatch"} or
    \code{"ExpressionSet"}}
  \item{onlyFirst}{ Boolean defining of only first or all instances
    found should be returned. Default is set to TRUE. Warnings are
    returned if more than one probe is found. The function stops
    with an error if no probe is found.}
}

\value{
  An object of type \code{"character"} with the Affymetrix hybridation
  (bio) probe names for the given chip type. }

\examples{
    ## load a dataset
    library(affydata)
    data(Dilution)
    getBioProbes(Dilution)
}

\seealso{getSpikeProbes,getRatioProbes}

\author{ Laurent Gatto }

\keyword{ misc }
