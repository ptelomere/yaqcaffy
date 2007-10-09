\name{getSpikeProbes}
\alias{getSpikeProbes}

\title{ Get the names of all spike probes on the array }
\description{
    This function returns all the spike probes (i.e. BioB-3', BioD-5', Lys-3, ...) 
    that are located on the given GeneChip.
}
\usage{
  getSpikeProbes(object)
}

\arguments{
   \item{object}{ An object of type \code{AffyBatch} or \code{ExpressionSet}. }
}

\value{
  An object of class \code{character} containing all Affymetrics spike 
  probe names.
}

\author{ Laurent Gatto }

\examples{
    ## load a dataset
    library(affydata)
    data(Dilution)
    getSpikeProbes(Dilution)
}
\seealso{getBioProbes,getRatiosProbes}

\keyword{ misc }
