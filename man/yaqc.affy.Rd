\name{yaqc.affy}
\alias{yaqc.affy}
\alias{yaqc}

\title{Generate QC stats from an AffyBatch object }

\description{Generate YAQC metrics for Affymetrix data.}

\details{
Affymetrix recommend a set of quality control metrics to check the quality of 
GeneChips expression arrays. This function applies the guidelines described 
in the Affymetrix Microarray Quality Control Consortium (MAQC) protocols to 
asses the succes of the hybridization. See the package vignette for more 
details.

This function takes a raw (unnormalised) \code{AffyBatch} object or
and \code{ExpressionSet} object. In the first case, it computess GCOS intensity 
values, MAS expression calls(see \code{call.exprs}) and other quality-related 
metrics to generate an YAQCStats. If an \code{ExpressionSet} object is provided,
only the \eqn{\beta}-actin, GAPDH  and internal control values are computed.
}

\usage{
yaqc.affy(object,tgt=100,tau=0.015,logged=FALSE)
}

\arguments{
  \item{object}{a object of type \code{AffyBatch} or \code{ExpressionSet}}
  \item{tgt}{the target intensity to which the chips should be scaled 
            (used to calculated the GCOS intensities)}
  \item{tau}{used by detection p value}
  \item{logged}{to be used with an \code{ExpressionSet} object, defining if the
                expression intensities are logged}
}

\value{
  An YAQCStats object describing the input object
}

\author{ Laurent Gatto }

\examples{
    ## loading data
    library(affydata)
    data(Dilution)
    ## qc analysis
    qobj <- yaqc(Dilution)
    show(qobj)
}

\keyword{ misc }
