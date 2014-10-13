\name{getSpikeProbes}
\alias{getSpikeProbes}

\title{ Get the names of all spike probes on the array }
\description{
  This function returns all the spike probes
  (i.e. BioB-3', BioD-5', Lys-3, ...) that
  are located on the given GeneChip.
}
\usage{
  getSpikeProbes(object,onlyFirst)
}

\arguments{
   \item{object}{ An object of type \code{AffyBatch} or
     \code{ExpressionSet}. }
   \item{onlyFirst}{ Boolean defining of only first or all
     instances found should be returned. Default is set to TRUE.
     Warnings are returned if more than one probe is found. The
     function stops with an error if no probe is found.}
}

\value{
  An object of class \code{character} containing all (hybridization
  and labelling) Affymetrix spike probe names.
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
