\name{deg-methods}
\docType{methods}
\alias{deg-methods}

\title{ Returns the degradation probes of YaqcControlProbes or
  YaqcDegProbes objects.}

\description{
  The 'deg' method can be called with ojects of classes
  \code{\linkS4class{YaqcDegProbes}} or
  \code{\linkS4class{YaqcControlProbes}} (the latter containing the
  former) to get the labelling probes (actin3, actin5, actinM,
  gapdh5,...).
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "YaqcControlProbes")}}{ Returns the
      degratation probes as an instance of the
      \code{\linkS4class{YaqcDegProbes}}.} 
    \item{\code{signature(object = "YaqcDegProbes")}}{ Returns the
      degradation probes as characters.}
  }
}

\keyword{methods}

