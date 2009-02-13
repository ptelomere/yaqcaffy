\name{YAQCStats-class}
\docType{class}
\alias{YAQCStats-class}


\alias{isLog,YAQCStats-method}
\alias{target,YAQCStats-method}
\alias{avns,YAQCStats-method}
\alias{moreSpikeInProbes,YAQCStats-method}
\alias{gcosProbes,YAQCStats-method}
\alias{bioCalls,YAQCStats-method}

\alias{isLog-methods}
\alias{target-methods}
\alias{avns-methods}
\alias{moreSpikeInProbes-methods}
\alias{gcosProbes-methods}
\alias{bioCalls-methods}
\alias{arrays-methods}

\alias{plot,YAQCStats,ANY-method}
\alias{summary,YAQCStats,ANY-method}
\alias{show,YAQCStats-method}
\alias{merge,YAQCStats-method}
\alias{arrays,YAQCStats-method}

\title{Class "YAQCStats" }
\description{Holds Quality Control data for a set of Affymetrix arrays}
\section{Objects from the Class}{
    Objects can be created by calls of the form \code{yaqc(object)} where object 
    is of class \code{AffyBatch} of \code{ExpressionSet}. 
    YAQCStats is a subclass of \code{QCStats} and uses the \code{scale.factor}, 
    \code{average.background} and \code{percent.present}, \code{arraytype} 
    and \code{target} attributes of it's subclass.
}
\section{Slots}{
  \describe{
    \item{\code{log}:}{Object of class \code{"logical"} that specifies if expression values are in log2 form}
    \item{\code{average.noise}:}{Object of class \code{"numeric"} The average noise for the arrays }
    \item{\code{morespikes}:}{Object of class \code{"matrix"} More spiked in probes (e.g. r2biob5, r2biob3,...) }
    \item{\code{gcos.probes}:}{Object of class \code{"matrix"} GAPDH and \eqn{\beta}-actin qc probes 
    (e.g. gapdh 3,5,M,...) containing the GCOS values}
    \item{\code{bio.calls}:}{Object of class \code{"matrix"} BioB 5',3',M and BioC, BioC 5',3' 
    present/absent/marginal calls }
  }
}

\section{Methods}{
  \describe{
    \item{isLog}{\code{signature(object = "YAQCStats")}: are intensities in log2 form}
    \item{avns}{\code{signature(object = "YAQCStats")}: average noise}
    \item{moreSpikeInProbes}{\code{signature(object = "YAQCStats")}: more spike-in probes }
    \item{gcosProbes}{\code{signature(object = "YAQCStats")}: GCOS values of gapdh and actin QC probes }
    \item{bioCalls}{\code{signature(object = "YAQCStats")}: P/M/A calls for BioB, BioC and BioD spikes  }
    \item{arrays}{\code{signature(object = "YAQCStats")}: shows the array names composing the \code{YAQCStats} object }
    \item{plot}{\code{signature(object = "YAQCStats")}: visual representation of the qc metrics }
    \item{summary}{\code{signature(object = "YAQCStats")}: data frame summary of the outliers }
    \item{show}{\code{signature(object = "YAQCStats")}: displays the content of the object as a data frame }
    \item{merge}{\code{signature(object = "YAQCStats")}: merges two compatible YAQCStats objects, i.e. that have the same values for the \code{log}, \code{target} and \code{arraytype} slots. }
    \item{arrays}{\code{signature(object = "YAQCStats")}: show the array names of an YAQCStats objects. }
  }
}

\author{ Laurent Gatto }

\seealso{QCStats-class}

\keyword{classes}
