\name{YAQCStats-class}
\docType{class}
\alias{YAQCStats-class}
\alias{class:YAQCStats}

\alias{isLog,YAQCStats-method}
\alias{isLog}
\alias{target,YAQCStats-method}
\alias{avns,YAQCStats-method}
\alias{avns}
\alias{avbg,YAQCStats-method}
\alias{minbg,YAQCStats-method}
\alias{percent.present,YAQCStats-method}
\alias{sfs,YAQCStats-method}
\alias{moreSpikeInProbes,YAQCStats-method}
\alias{moreSpikeInProbes}
\alias{gcosProbes,YAQCStats-method}
\alias{gcosProbes}
\alias{bioCalls,YAQCStats-method}
\alias{bioCalls}
\alias{objectVersion,YAQCStats-method}
\alias{objectVersion}
\alias{getYaqcControlProbes,YAQCStats-method}
\alias{getYaqcControlProbes}

\alias{objectVersion,YAQCStats-method}
\alias{getYaqcControlProbes,YAQCStats-method}
\alias{getYaqcControlProbes}

\alias{show,YAQCStats-method}
\alias{yaqc.summary}
\alias{summary,YAQCStats-method}
\alias{merge,YAQCStats,YAQCStats-method}
\alias{merge,YAQCStats-method}
\alias{arrays,YAQCStats-method}
\alias{arrays}

\title{Class "YAQCStats" }

\description{Holds Quality Control Data for a set of Affymetrix Arrays}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{\link{yaqc}}
  providing \code{AffyBatch} or \code{ExpressionSet} instances as
  arguments. 
  YAQCStats is a subclass of \code{QCStats} and uses the
  \code{scale.factor}, \code{average.background}, 
  \code{percent.present}, \code{arraytype} and \code{target}
  attributes of it's super-class. 
}

\section{Slots}{
  Class-specific slots:
  \describe{    
    \item{\code{log}:}{Object of class \code{"logical"} that
      specifies if expression values are in log2 form.}
    \item{\code{average.noise}:}{Object of class \code{"numeric"}.
      The average noise for the arrays.}
    \item{\code{morespikes}:}{Object of class \code{"matrix"}.
      More spiked in probes (e.g. r2biob5, r2biob3,...).}
    \item{\code{gcos.probes}:}{Object of class \code{"matrix"}.
      GAPDH and \eqn{\beta}-actin qc probes 
      (e.g. gapdh 3,5,M,...) containing the GCOS values.}
    \item{\code{bio.calls}:}{Object of class \code{"matrix"}.
      BioB 5',3',M and BioC, BioC 5',3' 
      present/absent/marginal calls .}
    \item{\code{objectVersion}}{Character describing the version
      of the library used to generate the \code{YAQCStats} object.}
    \item{\code{yaqcControlProbes}}{Object of class \code{YaqcControlProbes}
      that defines the different probes used for the quality
      control.}
  }
  See also \code{\linkS4class{QCStats}} for slots inherited from
  super-class. 
}

\section{Methods}{
  Methods inherited from the super-class:
  \describe{
    \item{target}{\code{signature(object = "YAQCStats")}: Returns a
      numeric target value for MAS 5.0 normalization. }
    \item{avbg}{\code{signature(object = "YAQCStats")}: Returns a vector
      of the average background levels for each array. }
    \item{minbg}{\code{signature(object = "YAQCStats")}: Returns a
      vector of the minimum background levels for each array. }
    \item{percent.present}{\code{signature(object = "YAQCStats")}:
      Returns a vector listing the percentage of probesets called
      present on each array. }
    \item{sfs}{\code{signature(object = "YAQCStats")}: Returns a vector
      of scale factors for each array (as produced by the MAS 5.0
      algorithm). } 
  }  
  Class-specific methods:
  \describe{
    \item{isLog}{\code{signature(object = "YAQCStats")}: Returns a
      logical specifying if the expression intensities are in log2
      from. }    
    \item{moreSpikeInProbes}{\code{signature(object = "YAQCStats")}:
      Returns a \code{matrix} of intensities for the internal spike
      probes. }  
    \item{gcosProbes}{\code{signature(object = "YAQCStats")}: Returns a
      \code{matrix} of intensities for GAPDH and \eqn{\beta}-actin
      probes. } 
    \item{avns}{\code{signature(object = "YAQCStats")}: Returns a vector
      listing the average noise levels for each array. }
    \item{bioCalls}{\code{signature(object = "YAQCStats")}: Returns a
      \code{matrix} of Present(P)/Marginal(M)/Absent(A) calls for the spike
      probes. }
    \item{arrays}{\code{signature(object = "YAQCStats")}: Returns the
      names of the arrays in the \code{YAQCStats} instance. }
    \item{plot}{\code{signature(object = "YAQCStats")}: visual
      representation of the qc metrics. (see \code{\link{yaqc.plot}} for
      more details). }     
    \item{summary}{\code{signature(x = "YAQCStats", latex = "logical")}:
      The outliers of the \code{\linkS4class{YAQCStats}} quality
      control metrics are summerized and returned as a data frame. If
      \code{latex} is set to \code{TRUE} (default), the data frame is
      returned as a latex table (requires the \code{xtable} package). }
    \item{show}{\code{signature(object = "YAQCStats")}: displays the content
      of the object as a data frame.} 
    \item{merge}{\code{signature(x = "YAQCStats", y ="YAQCStats")}:
      merges two compatible YAQCStats objects, i.e. that have the same
      values for the \code{log}, \code{target} and \code{arraytype}
      slots.} 
    \item{arrays}{\code{signature(object = "YAQCStats")}: shows the array
      names of an YAQCStats objects. }
    \item{objectVersion}{\code{signature(object = "YAQCStats")}: Returns
      the version of the \code{yaqcaffy} package as a \code{character}
      used to create the given object. }
    \item{getYaqcControlProbes}{\code{signature(object = "YAQCStats")}:
      Returns the \code{\linkS4class{YaqcControlProbes}} object that has
      been used to generate the current \code{YAQCStats} object.} 
  }
}

\author{ Laurent Gatto }

\seealso{
  \code{\linkS4class{QCStats}} from package \code{simpleaffy} and
  \code{\linkS4class{YaqcControlProbes}}.
}

\keyword{classes}
