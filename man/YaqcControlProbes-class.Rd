\name{YaqcControlProbes-class}
\Rdversion{1.1}
\docType{class}
\alias{YaqcControlProbes-class}
\alias{class:YaqcControlProbes}
\alias{YaqcBioProbes-class}
\alias{class:YaqcBioProbes}
\alias{YaqcSpkProbes-class}
\alias{class:YaqcSpkProbes}
\alias{YaqcDegProbes-class}
\alias{class:YaqcDegProbes}
\alias{bio,YaqcControlProbes-method}
\alias{bio,YaqcBioProbes-method}
\alias{bio-methods}
\alias{bio}
\alias{deg,YaqcControlProbes-method}
\alias{deg,YaqcDegProbes-method}
\alias{deg-methods}
\alias{deg}
\alias{info<-,YaqcControlProbes-method}
\alias{info,YaqcControlProbes-method}
\alias{info}
\alias{info<-}
\alias{show,YaqcControlProbes-method}
\alias{show,YaqcBioProbes-method}
\alias{show,YaqcSpkProbes-method}
\alias{show,YaqcDegProbes-method}
\alias{spk,YaqcControlProbes-method}
\alias{spk,YaqcSpkProbes-method}
\alias{spk-methods}
\alias{spk}

\title{Class "YaqcControlProbes" }
\description{
  The \code{\linkS4class{YAQCStats}} class stores the probes used for the
  quality control as a special class, namely the
  \code{YaqcControlProbes} class. This class encapsulated
  the probe names that are used the generate an
  \code{\linkS4class{YAQCStats}} object. Objects of this class are
  created with the \code{\link{probeSelectionInterface}} function.  The
  hybridization, labelling and degradation probes are encapsulated in
  \code{YaqcBioProbes}, \code{YaqcSpkProbes} and \code{YaqcDegProbes}
  objects respectively. These can be retrieved with their respective
  accessors, as described below.
  Furthermore, an \code{info} function allows
  to retrieve or set a free text slot to describe the
  \code{YaqcControlProbes} object.
}

\section{Slots}{
  \describe{
    \item{\code{bio}:}{Object of class \code{"YaqcBioProbes"}
      encapsulating the 'bio' (BioB5, BioB3, BioBM, BioC5,...)
      hybridization probes.}
    \item{\code{spk}:}{Object of class \code{"YaqcSpkProbes"}
      encapsulating the labelling probes (dap5, dap3, dap3, phe5,...).} 
    \item{\code{deg}:}{Object of class \code{"YaqcDegProbes"}
      encapsulating the degradation probes used to assess the 3'/5'
      ratio.}
    \item{\code{info}:}{Object of class \code{"character"} providing
      general information about the \code{YaqcControlProbes}
      object.}
  }
}

\section{Methods}{
  \describe{
    \item{bio}{\code{signature(object = "YaqcControlProbes")}: returns
      the 'bio' (BioB5, BioB3, BioBM, BioC5,...) hybridization probes of
      the current object, as an \code{YaqcBioProbes} instance. }
    \item{bio}{\code{signature(object = "YaqcBioProbes")}: returns
      the 'bio' (BioB5, BioB3, BioBM, BioC5,...) hybridization probes of
      the current object, as characters. }
    \item{spk}{\code{signature(object = "YaqcControlProbes")}: returns
      the labelling probes (dap5, dap3, dap3, phe5,...) of the current
      object, as an \code{YaqcSpkProbes} instance. }
    \item{spk}{\code{signature(object = "YaqcSpkProbes")}: returns
      the labelling probes (dap5, dap3, dap3, phe5,...) of the current
      object, as characters}
    \item{deg}{\code{signature(object = "YaqcControlProbes")}: returns 
      the degradation probes used to assess the 3'/5'ratio, as an
      \code{YaqcDegProbes} instance. }
    \item{deg}{\code{signature(object = "YaqcDegProbes")}: returns 
      the degradation probes used to assess the 3'/5'ratio, as
      characters. }    
    \item{info}{\code{signature(object = "YaqcControlProbes")}: returns
      the inforation slot of the current object.}
    \item{info<-}{\code{signature(object = "YaqcControlProbes")}: sets
      the information slot of the current object.}
    \item{show}{\code{signature(object = "YaqcControlProbes")}: shows
      the current object.}
    }
}
\author{ Laurent Gatto }

\seealso{
  \code{\link{probeSelectionInterface}} and 
  \code{\linkS4class{YAQCStats}}
}

\examples{
showClass("YaqcControlProbes")
showClass("YaqcBioProbes")
showClass("YaqcSpkProbes")
showClass("YaqcDegProbes")
}

\keyword{classes}
