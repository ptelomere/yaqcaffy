\name{spk-methods}
\docType{methods}
\alias{spk-methods}


\title{ Returns the labelling probes of YaqcControlProbes or
  YaqcSpkProbes objects.}

\description{
  The 'spk' method can be called with ojects of classes
  \code{\linkS4class{YaqcSpkProbes}} or
  \code{\linkS4class{YaqcControlProbes}} (the latter containing the
  former) to get the hybridization probes (dap3, dap5, dapM, phe5,...).
}


\section{Methods}{
  \describe{
    \item{\code{signature(object = "YaqcControlProbes")}}{ Returns the
      hybridization probes as an instance of the
      \code{\linkS4class{YaqcSpkProbes}}.}
    \item{\code{signature(object = "YaqcSpkProbes")}}{ Returns the
      hybridization probes as characters.}
  }
}

\keyword{methods}
