\name{yaqc.affy}
\alias{yaqc.affy}
\alias{yaqc}
\alias{yaqc-methods}
\alias{yaqc-method}
\alias{method:yaqc}
\alias{yaqc,eSet-method}
\alias{yaqc,ExpressionSet-method}
\alias{yaqc,AffyBatch-method}

\title{'YAQCStats' constructor}

\description{Generate YAQC metrics for Affymetrix data.}

\details{
  Affymetrix recommends a set of quality control metrics to check
  the quality of GeneChips expression arrays. This function applies
  the guidelines described in the Affymetrix Microarray Quality
  Control Consortium (MAQC) protocols to asses the succes of the
  hybridization. See the package vignette for more details.
  
  This function takes a raw (unnormalised) \code{AffyBatch} object
  or an \code{ExpressionSet} object. In the first case, it computes
  MAS5 intensity values, expression calls(see \code{call.exprs})
  and other quality-related metrics to generate an YAQCStats. If
  an \code{ExpressionSet} object is provided, only the
  \eqn{\beta}-actin, GAPDH  and internal control values are
  computed.
}

\usage{ yaqc.affy(object,
myYaqcControlProbes=NULL,
alphas=NULL,
tgt=100,
tau=0.015,
logged,
verbose) }

\arguments{
  \item{object}{a object of type \code{AffyBatch} or
    \code{ExpressionSet}.}
  \item{myYaqcControlProbes}{ an object of type
    \code{YaqcControlProbes}. If none is supplied (default behaviour), the
    control probes are selected automatically. See the
    \code{YaqcControlProbes} class for more details
    \code{probeSelectionInterface} to generate such an object.}
  \item{alphas}{a numeric of length 2 with the alpha1 and alpha2
    values. Alpha1 and alpha2 are thresholds used to define if a given
    probe shoudl be called present (p<alpha1), marginal
    (alpha1<p<alpha2) or absent (alpha2<p), where p is the p-value
    from the Wilcoxon Signed Rank test used in MAS5. The default is
    to get these values from simpleaffy's qcdef files
    (see simpleaffy's vignette for more details) or to use 0.04 and
    0.06 as default values.}
  \item{tgt}{the target intensity to which the chips should be scaled 
    (used to calculated the MAS5 intensitiy values).}
  \item{tau}{used for detection p-value.}
  \item{logged}{to be used with an \code{ExpressionSet} object, defining
    if the expression intensities are logged.}
  \item{verbose}{logical value. If 'TRUE', it writes out some messages
    indicating progress. If 'FALSE' nothing should be printed.}
}

\value{ An \code{YAQCStats} object describing the input object }

\author{ Laurent Gatto }

\examples{
    ## loading data
    library(affydata)
    data(Dilution)
    ## qc analysis
    qobj <- yaqc(Dilution)
    show(qobj)
}

\keyword{ methods }
