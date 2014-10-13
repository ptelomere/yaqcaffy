## This function checks if the users input array is of type
## hgu133plus2 before merging it using .mergeAffyBatchObjects()
.isHgu133plus2<-function(userAffyBatchObject) {
	if (cdfName(userAffyBatchObject) == "HG-U133_Plus_2") return(TRUE)
	else return(FALSE)
}


## This function merges the users Affymetrix hgu133plus2 with the
## corresponding RNA reference set, saved in the data directory 
## as an R object
.mergeAffyBatchObjects <- function(userAffyBatchObject,ref) {
    if(.isHgu133plus2(userAffyBatchObject)) {
	     if (ref=="refA") refData <- eval(parse(text=data(refA)))
	else if (ref=="refB") refData <- eval(parse(text=data(refB)))
	else if (ref=="refC") refData <- eval(parse(text=data(refC)))
	else if (ref=="refD") refData <- eval(parse(text=data(refD)))
        else if (ref=="test") {
          return(eval(parse(text=data(refA)))[,1:2])
        }
	else {
	    cat("Reference RNA must be refA, refB, refC or refD.\n")
	    cat("See the vignette for more details\n")
	}
    } else {
	cat("Your array must be a Human Genome U133 Plus 2.0\n")
	cat("to be compared with the MAQC reference dataset.\n")
    }
    return(merge(userAffyBatchObject,refData))
}


## This function calculates and plots the repoductibility statistics,
## i.e the coefficient variations and the Pearson correlation factors.
reprodPlot <- function (userAffyBatchObject,ref,
                        normalize=c("rma","gcrma","mas5","none"),
                        main="MAQC reference reproducibility",
                        cex=1,...) {
  ## checking if MAQCsubset is available
  require("MAQCsubsetAFX") || stop("Requires MAQCsubsetAFX to continue")
  
  ## preparing the data
  data <- .mergeAffyBatchObjects(userAffyBatchObject,ref)

  normalize <- match.arg(normalize)
       if (normalize=="gcrma") e<-exprs(gcrma(data))
  else if (normalize=="mas5")  e<-log2(exprs(mas5(data)))
  else if (normalize=="none")  e<-log2(exprs(data))
  else                         e<-exprs(rma(data))

  dim <- ncol(e)   ## plot dimension
  cor <- cor(e)    ## Pearson correlation factors
  labels <- sampleNames(data)
  
  old.par<-par(no.readonly=TRUE) ## save original parameters
  on.exit(par(old.par))
  par(mfrow=c(dim,dim),mgp=c(0,0.2,0),mar=c(0,0,0,0),oma=c(1.8,1.8,1.8,1.8))

  
  ## plotting
  for (i in 1:(dim-1)){
    ## plotting sample names on diagonal
    par(mfg=c(i,i))
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
    text(1,1,labels[i],cex=1.2*cex)
    for (j in (i+1):dim) {
      ## correlation factors above the diagonale
      par(mfg=c(i,j))
      txt <- format(cor[i,j],digits=3)
      plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
      text(1,1,txt,cex=1.3*cex)
      ## scatter plots below the diagonale
      par(mfg=c(j,i))
           if (i==1 && j%%2!=0) smoothScatter(e[,c(i,j)],xlab="",ylab="",xaxt="n")
      else if (j==dim && i%%2==0) smoothScatter(e[,c(i,j)],xlab="",ylab="",yaxt="n")
      else smoothScatter(e[,c(i,j)],xlab="",ylab="",xaxt="n",yaxt="n")
      .plotdiag()
    }    
  }

  ## last sample name  
  par(mfg=c(dim,dim))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  text(1,1,labels[dim],cex=cex)
  ## axis labels
  mtext(main,3,outer=TRUE,cex=1.4*cex)
}



.plotdiag <- function() {
  fc <- c(2,4,8)
  k <- 1
  for ( i in fc) {
    ## upper diagonal and text
    abline(log2(i),1,lwd=0.7,col="darkgray")
    ## lower diagonal and text
    abline(-log2(i),1,lwd=0.7,col="darkgray")
  }
}
