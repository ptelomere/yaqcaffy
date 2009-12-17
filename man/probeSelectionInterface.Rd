\name{probeSelectionInterface}
\alias{probeSelectionInterface}

\title{Tcltk Interface to Generate an Instance of YaqcControlProbes for a given Chip Set}

\description{
  \code{probeSelectionInterface} starts a tcltk graphical user
  interface (GUI) that allows the user to choose the probes to be used
  for subsequent quality analyses with the \code{yaqcaffy} package. The
  probes are selected on basis of the features of a given set of
  Affymetrix Genechips provided as input.  The list of probes can be
  pre-filtered to display only control probes (i.e starting by AFFX) or
  all probes on the Genechip can be shown.
}

\usage{
probeSelectionInterface(object,
returnVar="yaqcControlProbes",
filter=TRUE)
}

\arguments{
  \item{object}{an object of class \code{linkS4Class{AffyBatch}} or
    \code{linkS4Class{ExpressionSet}}.}  
  \item{returnVar}{a string defining the name of the variable the
    returned object will be saved as in the global environment. The
    default variable name is 'yaqcControlProbes'. If such a variable
    name already exists, a warning will be issued and the user can
    cancel the function.}
  \item{filter}{logical value. If 'TRUE', the feature names of the
    input object are filtered out (see details). If 'FALSE', all
    features are listed for all control probes.}
}
 
\details{
  Three tabs are displayed, one for the hybridization (bio) probes,
  labelling probes (dap, phe, thr and lys) and the degradation probes
  (actin and gapdh) respectively.  If the user uses the 'Close' button,
  no return object is saved in the global environment.  An object is
  saved as \code{returnVar} if the user presses 'Ok'.  If such a variable
  name already exists, a warning will be issued and the user can close
  the interface and cancel the function.

  If filtering is applied, the hybridization menus will list probes
  that match the given probe (BioB, BioC or BioD) and position (5, 3 or
  M).  Similarly, only matching labelling probes (dap, phe, thr and lys)
  and positions will be displayed.  As the pattern for the degradation
  probes are less strict, all the 'AFFX' probes, except those already
  selected as hybridization and labelling probes, will be displayed
  in the drop-down menus.
}

\value{Returns an object of class \code{\linkS4class{YaqcBioProbes}}. }

\examples{
 \dontrun{
   library(affydata)
   data(Dilution)
   probeSelectionInterface(Dilution)
 }
}

\author{Laurent Gatto}

\keyword{misc}
