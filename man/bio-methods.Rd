\name{bio-methods}
\docType{methods}
\alias{bio-methods}

\title{ Returns the 'bio' hybridization probes of YaqcControlProbes or
  YaqcBioProbes objects.} 

\description{
  The 'bio' method can be called with ojects of classes
  \code{\linkS4class{YaqcBioProbes}} or
  \code{\linkS4class{YaqcControlProbes}} (the latter containing the
  former) to get the labelling probes (BioB3, BioB5, BioBM, BioC5,...).  
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "YaqcControlProbes")}}{ Returns the labelling
      probes as an instance of the \code{\linkS4class{YaqcBioProbes}}.}
    \item{\code{signature(object = "YaqcBioProbes")}}{ Returns the labelling
      probes as characters.}
  }
}

\keyword{methods}
