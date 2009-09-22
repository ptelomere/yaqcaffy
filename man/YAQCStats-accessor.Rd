\title{Accessor Methods for YAQCStats class}

\name{moreSpikeInProbes}
\alias{isLog}
\alias{avns}
\alias{moreSpikeInProbes}
\alias{gcosProbes}
\alias{bioCalls}
\alias{arrays}
\alias{plot}

\usage{
isLog(object)
avns(object)
moreSpikeInProbes(object)
gcosProbes(object)
bioCalls(object)
}

\description{
The Microarray Quality Control Consortium guidelines and reference datasets 
are applied to analyse some of the QC metrics recommended by Affymetrix. 
These guidelines, datasets and QC metrics are discussed in detail in the 
package vignette. 

Each of these functions queries a YAQCStats object to retrieve one of
these QC metrics. The YAQCStats class is a subclass of the QCStats object 
as defined in the \pkg{simpleaffy} package. \code{sfs}, 
\code{avbg} and \code{percent.present} methods access to the
\code{scale.factors},  \code{average.background} and
\code{percent.present} attributes respectivley  are defined in the superclass.

Briefly:

isLog:             Returns a logical saying if the expression
intensities are in log2 from.

target:            Returns the target value (for MAS 5.0 normalization).

moreSpikeInProbes: Returns a table of intensities for the internal spike
probes.

gcosProbes:        Returns a table of intensities for GAPDH and
\eqn{\beta}-actin probes.

avbg:              Returns a vector of the average background levels for
each array.

minbg:             Returns a vector of the minimum background levels for
each array.

percent.present:   Returns a vector listing the percentage of probesets
called present on each array.

sfs:               Returns a vector of scale factors for each array (as
produced by the MAS 5.0 algorithm).

bioCalls:          Returns a table of Present/Marginal/Absent calls for
the spike probes.

arrays:            Returns the names of the arrays in the YAQCStats
object.

plot:              Plots the YACQStats object (see
\code{\link[yaqcaffy:yaqc.affy]{yaqc}} for more details). 
}

\seealso{\code{\link[yaqcaffy:yaqc.affy]{yaqc}}, \code{QCStats}}  

\arguments{
  \item{object}{an object of class \code{"YAQCStats"}}
}

\author{Laurnet Gatto}

\keyword{ misc }
