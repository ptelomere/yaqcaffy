\name{yaqc.plot}
\alias{yaqc.plot}
\alias{plot,YAQCStats}
\alias{plot,YAQCStats-method}
\alias{plot,YAQCStats,missing-method}


\title{Plots a YAQCStats object}
\description{
  Generates a visual summary of the various Affymetrix QC statistics.
}
\usage{yaqc.plot(YAQCStatsObject, which,...)}
\arguments{
    \item{YAQCStatsObject}{an object of class \code{"YAQCStats"}.}
    \item{which}{which quality metrics should be plotted (all by default).}
    \item{...}{other arguments.}
}
\details{
    The quality control metrics of the \code{YAQCStatsObject} are plotted 
    in a serie of graphs with the recommended ranges. 

    The scale factors are represented through a dot chart and the upper and
    lower limits are defined with vertical red lines.
    The other qc metircs are shown using dot plots. For the upper row box 
    plots (average background, average noise, percent present and 
    \eqn{\beta}-actin and GAPDH ratios, the mean is represented by a dashed 
    red line and the mean +/- 2 stdev by red dotted lines. For the lower box 
    plots, featuring the internal controls, grey boxes defines the mean 
    (middle segment) +/- 2 stdev.

    Individual plots can also be generated with the \code{which}
    argument: 'sfs' for the scale factor, 'avbg' and 'avns' for the
    average background and noise, 'pp' for the percentage of present
    calls, 'gapdh' and 'actin' for the GAPDH and \eqn{\beta}-actin
    ratios, 'bio' for the hybridization controls and 'spikes' for the
    retro-transciption spiked controls. 

    If the YAQCStatsObject has been generated with an Expression Set objects, 
    the scale factors, average noise and background and percent present can not
    be computed and the respective plots are removed from the final graph.
}

\examples{
  ## load data
  library(affydata)
  data(Dilution)
  ## create the yaqc object
  ## and plot it
  qobj <- yaqc(Dilution)
  plot(qobj)
}

\author{Laurent Gatto}
 
\keyword{misc}
