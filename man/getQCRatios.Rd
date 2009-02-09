\name{getQCRatios}
\alias{getQCRatios}
\title{Compute qc probe ratios using GCOS intensity values}
\description{
    This function computes the 3'/5' ratios of the GAPDh and \eqn{\beta}-actin
    qc probes using the GCOS intensity values.
}
\usage{getQCRatios(YAQCStatsObject)}
\arguments{
    \item{YAQCStatsObject}{ an object of class YAQCStats}
}

\value{An object of type \code{"matrix"} with two qc ratios per array.}

\seealso{getRatioProbes}

\examples{
    ## load a dataset
    library(affydata)
    data(Dilution)
    ## create yaqc object
    qobj <- yaqc(Dilution)
    getQCRatios(qobj)
}

\author{Laurent Gatto}

\keyword{misc}
