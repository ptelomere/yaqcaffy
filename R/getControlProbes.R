.checkNumberOfProbes <- function(probe,msg="") {
  if (length(probe)==0) stop(c(msg,"No probe found."),call.=FALSE)
  if (length(probe)>1) warning(c(msg,paste(probe,collapse=","),"\nMore than 1 probes found:"),call.=FALSE)
  invisible(length(probe))
}


## =================================================
## Instead of statically define the probes to be
## used for the quality analyses, as in
## simpleaffy (qcdef files) or yaqcaffy
## (morespike.tab), the functions below will
## 'automatically' select the probes based on
## defined patterns.
## 'r2' probes are selected preferentially. If no
## probes are found, a warning is issued and '!r2'
## probes are searched for. If no probes are found, 
## an error is issued.
## If more than 1 probe is found, the first is selected
## unless onlyFirst=FALSE is set and a warning
## is issued.
## - - - - - - - - - - - - - - - - - - - - - - - - -

getBio <- function(x,y,affyobject,onlyFirst=TRUE) {
  x <- casefold(x,upper=TRUE)
  y <- casefold(y,upper=TRUE)
  if (!all(x %in% c("B","C","D"))) stop("Only bio[B|C|D] probes available\n")
  if (!all(y %in% c("3","M","5"))) stop("Only [3|M|5] probes available\n")
  d <- featureNames(affyobject)
  toReturn <- c()
  for (i in x) {
    for (j in y) {
      ## matching r2 probes
      pattern <- paste("AFFX-.*r2.+[B|b]io",i,".+",j,sep="")
      probe <- grep(pattern,d,value=TRUE,ignore.case=TRUE) 
      if (length(probe)==0) {
        ## matching non-r2 probes
        warning(paste("No Bio",x,y," 'r2' probe found. Searching for !r2 probe.\n",sep=""),call.=FALSE)
        pattern <- paste("AFFX.+[B|b]io",i,".+",j,sep="")
        probe <- grep(pattern,d,value=TRUE,ignore.case=TRUE) 
      }
      .checkNumberOfProbes(probe,paste("bio",i,j,": ",sep=""))
      if (onlyFirst) probe <- probe[1]
      toReturn <- c(toReturn,probe)
    }
  }
  return(toReturn)
}


getSpk <- function(x,y,affyobject,onlyFirst=TRUE) {
  x <- casefold(x,upper=FALSE)
  y <- casefold(y,upper=TRUE)
  if (!all(x %in% c("dap","phe","thr","lys")))
    stop("Only [dap|phe|thr|lys] probes available\n")
  if (!all(y %in% c("3","M","5")))
    stop("Only [3|M|5] probes available\n")
  d <- featureNames(affyobject)
  toReturn <- c()
  for (i in x) {
    for (j in y) {
      ## matching r2 probes
      pattern <- paste("AFFX-.*r2.+",i,".+",j,sep="")
      probe <- grep(pattern,d,value=TRUE,ignore.case=TRUE) 
      if (length(probe)==0) {
        ## matching non-r2 probes
        warning(paste("No ",x,y," 'r2' probe found. Searching for !r2 probe.\n",sep=""),call.=FALSE)
        pattern <- paste("AFFX.+",i,".+",j,sep="") 
        probe <- grep(pattern,d,value=TRUE,ignore.case=TRUE) 
      }
      .checkNumberOfProbes(probe,paste(i,j,": ",sep=""))
      if (onlyFirst) probe <- probe[1]
      toReturn <- c(toReturn,probe)      
    }
  }
  return(toReturn)
}

getDeg <- function(x,y,affyobject,onlyFirst=TRUE) {
  x <- casefold(x,upper=FALSE)
  y <- casefold(y,upper=TRUE)
  if (!all((x %in% c("actin","gapdh"))))
    stop("Only [actin|gapdh] probes available\n")
  if (!all((y %in% c("3","M","5"))))
    stop("Only [3|M|5] probes available\n")
  ## special cases: human
  if (substr(cleancdfname(cdfName(affyobject)),1,3)=="hgu") {
    if (any(sel <- (x %in% "actin"))) x[sel] <- "hsac"
  }
  ## special cases: anopheles/plasmodium
  if (cleancdfname(cdfName(affyobject))=="plasmodiumanophelescdf") {
    if (any(sel <- (x %in% "actin"))) x[sel] <- "act"
  }
  ## special cases: Xenopus laevis
  if (cleancdfname(cdfName(affyobject))=="xenopuslaeviscdf") {
    if (any(sel <- (x %in% "actin"))) x[sel] <- "act"
    if (any(sel <- (x %in% "gapdh"))) x[sel] <- "gapd"
  }
  ## special cases: yeast_2
  if (substr(cleancdfname(cdfName(affyobject)),1,5)=="yeast") {
    if (any(sel <- (x %in% "actin"))) x[sel] <- "YFL039"
    if (any(sel <- (x %in% "gapdh"))) x[sel] <- "YER022"
    warning("Yeast: using RNApolII instead of GAPDH as degradation control!",
            call.=FALSE)
  }
  ## special cases: ygs98
  if (cleancdfname(cdfName(affyobject))=="ygs98cdf") {
    if (any(sel <- (x %in% "actin"))) x[sel] <- "YFL039"
    if (any(sel <- (x %in% "gapdh"))) x[sel] <- "YER022"
    warning("Yeast: using RNApolII instead of GAPDH as degradation control!",
            call.=FALSE)
  }
  ## special cases: bovine
  if (cleancdfname(cdfName(affyobject))=="bovinecdf") {
    x <- substr(x,1,4) 
  }
  d <- featureNames(affyobject)
  toReturn <- c()
  for (i in x) {
    for (j in y) {
      ## matching r2 probes
      pattern <- paste("affx-.*r2-.+",i,".*[-|_]",j,"_",sep="")
      probe <- grep(pattern,d,value=FALSE,ignore.case=TRUE) 
      if (length(probe)==0) {
        ## matching non-r2 probes
        pattern <- paste("affx.+",i,".*",j,"_",sep="")
        probe <- grep(pattern,d,value=FALSE,ignore.case=TRUE) 
      }
      probe <- featureNames(affyobject)[probe]
      .checkNumberOfProbes(probe,paste(i,j,": ",sep=""))
      if (onlyFirst) probe <- probe[1]
      toReturn <- c(toReturn,probe)
    }
  }
  return(toReturn)
}
