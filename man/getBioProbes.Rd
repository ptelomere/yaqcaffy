\name{getBioProbes}
\alias{getBioProbes}

\title{ Get the names of the Bio spike probes on the array }
\description{
    This function returns all the AFFX-Bio probes names that are located on the given GeneChip.
}
\usage{
  getBioProbes(affyBatchObject)
}

\arguments{
  \item{affyBatchObject}{ An object of type \code{"AffyBatch"} }
}

\value{
  An object of type \code{"character"} with all the Affymetrix Bio probe names. 
}

\examples{
    ## load a dataset
    library(affydata)
    data(Dilution)
    getBioProbes(Dilution)
}

\seealso{getSpikeProbes,getRatioProbes}

\author{ Laurent Gatto }

\keyword{ misc }
