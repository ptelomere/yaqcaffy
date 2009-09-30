\name{reprodPlot}
\alias{reprodPlot}
\title{Plot human whole genome GeneChips reproductibility}
\description{
    Compares Affymetrix Human Genome U133 Plus 2.0 Arrays to a subset 
    of the MAQC arrays for a RNA reference.
}
\usage{reprodPlot(userAffyBatchObject,ref,
                  normalize=c("rma","gcrma","mas5","none"),
                  main="MAQC reference reproducibility",
                  cex,...)}
\arguments{
    \item{userAffyBatchObject}{a set of Human Genomue U133 Plus 2.0 arrays 
    provided as an AffBatch object,}
    \item{ref}{a string ("refA", "refB", "refC", or "refD") defining the RNA 
    reference to compare the userAffyBatchObject to,}
    \item{normalize}{a string defining the algorithm used for data
      normalization:  
    \code{rma} (default) for RMA (as imlplemented in the \pkg{affy} library), 
    \code{gcrma} for GCRMA (as imlplemented in the
    \code{\link[gcrma]{gcrma}} library),  
    \code{mas5} for MAS5 (as imlplemented in the \pkg{affy} library) 
    or 'none' for no normalization,}
    \item{main}{an overall title for the plot,}
    \item{cex}{size of text on the plot,}
    \item{...}{other arguments.}
}
\details{
   The plot shows all the pairwise scatterplots (plotted with 
   \pkg{graphics})'s \code{\link[graphics]{smoothScatter}}
   function)  
   with Pearson's correlation factor and MAplots (plotted with 
   \pkg{affy}'s \code{\link[affy:MAplot]{ma.plot}} function).
   The subset of the MAQC arrays are 1 randomly chosen .CEL file out 
   of the 5 replicates for the 6 different test site.
}
\value{Outputs a graph on the available graphical device}

\examples{
 \dontrun{
   ## loading data
   library(MAQCsubsetAFX)
   data(refB)
   d<-refB[,1]
   ## testing the reproductibility against ref A
   reprodPlot(d,"refA",normalize="rma")
 }
}

\author{Laurent Gatto}

\keyword{misc}
