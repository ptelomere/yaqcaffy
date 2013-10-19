## accessor methods for the YAQCStats class
setMethod("avns",             "YAQCStats",function(object) object@average.noise)
setMethod("moreSpikeInProbes","YAQCStats",function(object) object@morespikes)
setMethod("gcosProbes",       "YAQCStats",function(object) object@gcos.probes)
setMethod("bioCalls",         "YAQCStats",function(object) object@bio.calls)
setMethod("arrays",           "YAQCStats",function(object) names(object@average.noise))
setMethod("isLog",            "YAQCStats",function(object) object@log)
setMethod("getYaqcControlProbes","YAQCStats",function(object) object@yaqcControlProbes)
setMethod("objectVersion",    "YAQCStats",function(object) object@objectVersion)


setMethod("summary","YAQCStats",
          function(object,...) yaqc.summary(object, ...))

yaqc.summary <- function(YAQCStatsObject, latex = FALSE) {
  if (latex) {
    require("xtable") || stop("'xtable' library is required to generate a latex summary table.")
  }
  qcm <- c("sfs","avbg","avns","pp","actin","gapdh")
  l <- vector("list",length(qcm))                   
  bio <- c("biob","bioc","biod")
  biol <- vector("list",3)
  spk <- c("dap","phe","lys","thr")
  spkl <- vector("list",4)
  
  for (i in 1:length(qcm)) {
    if (!is.null(getOutliers(YAQCStatsObject,qcm[i]))) {
      l[[i]] <- paste(names(getOutliers(YAQCStatsObject,qcm[i])),collapse=",")
    }
  }
  for (i in 1:length(bio)) {
    if (!is.null(getOutliers(YAQCStatsObject,bio[i])))
      biol[[i]] <- paste(names(getOutliers(YAQCStatsObject,bio[i])),collapse=",")
  }
  for (i in 1:length(spk)) {
    if (!is.null(getOutliers(YAQCStatsObject,spk[i])))
      spkl[[i]] <- paste(names(getOutliers(YAQCStatsObject,spk[i])),collapse=",")
  }
  l <- sub(".present","",l)
  df <- data.frame(rbind(as.matrix(l),paste(unlist(biol),collapse=" "),
                         paste(unlist(spkl),collapse=" ")
                         ))
  colnames(df) <- "samples"
  rownames(df) <- c(qcm,"bio","spikes")
  if (latex) { return(xtable(df)) }
  else { return(df) }
}


## This accessor method overloads simpleaffy's sfs().
## Here, we add the names to the numeric vector
setMethod("sfs","YAQCStats",
        function(object) {
        foo <- object@scale.factors
        names(foo) <- arrays(object)
        return(foo)
        })


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Changed in versions > 1.7.0 to use getBio, getSpk and getDeg
## functions that dynamically select the control probes based
## on hard-coded patterns
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##
## Functions similar to the ones implemented in the quality control 
## related file (qc.stats.R) of the simpleaffy package, but modifed 
## and renamed for the YAQCStats Object
##getSpikeProbes <- function(object) {
##  cdfn<-cleancdfname(object@annotation)
##  r <- get("morespikes",envir=.yaqcEnv)[cdfn,]
##  return(r[!is.na(r)])
##}
##getBioProbes <- function(affyBatchObject) {
##  cdfn<-cleancdfname(cdfName(affyBatchObject))
##  return(c(as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biob5"]),
##           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biobm"]),
##           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biob3"]),
##           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2bioc5"]),
##           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2bioc3"]),
##           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biod5"]),
##           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biod3"]))
##  )
##}
##getRatioProbes <- function(object) {
##  qlist <- getQcdefData(object)
##  actin <- c(qlist[[9]][3:4],qlist[[8]][4])
##  gapdh <- c(qlist[[11]][3:4],qlist[[10]][4])
##  return(c(actin,gapdh))
##}

## updated in versions > 1.7.0
getSpikeProbes <- function(object,onlyFirst=TRUE) {
  return(c(getBio("b",c("3","5","m"),object,onlyFirst=onlyFirst),
           getBio(c("c","d"),c("3","5"),object,onlyFirst=onlyFirst),
           getSpk(c("lys","phe","thr","dap"),c("3","5","m"),object,onlyFirst=onlyFirst))
         )
}

getSpkProbes <- function(object,onlyFirst=TRUE) {
  return(c(getSpk(c("lys","phe","thr","dap"),c("3","5","m"),object,onlyFirst=onlyFirst)))
}

## updated in versions > 1.7.0
getBioProbes <- function(object,onlyFirst=TRUE) {
  return(c(getBio("b",c("3","5","m"),object,onlyFirst=onlyFirst),
           getBio(c("c","d"),c("3","5"),object,onlyFirst=onlyFirst)))
}

## updated in versions > 1.7.0
getRatioProbes <- function(object,onlyFirst=TRUE) {
  return(c(getDeg(c("actin","gapdh"),c("3","5","m"),object,onlyFirst=onlyFirst)))
}
getDegProbes <- getRatioProbes



## - - - - - - - - - - - - - - - - - - - - - - - -
## Getting alphas from simpleaffy qcdef files
getAlpha1 <- function(object) {
  qlist <- getQcdefData(object)
  return(qlist[[2]][2])
}

getAlpha2 <- function(object) {
  qlist <- getQcdefData(object)
  return(qlist[[3]][2])
}

getQcdefData <- function(object) {
  ## this function retrieves a list of information
  ## about the specific type of affyBatchObject
  cdfn<-cleancdfname(object@annotation)
  fn <- paste(cdfn,"qcdef",sep=".")
  qcfile <- system.file("extdata",fn,package="simpleaffy")
  fl <- file(qcfile,"r")
  if(!isOpen(fl)) {
    warning(paste("Couldn't open file",qcfile,"."))
    return(NULL)
  }
  lines <- readLines(fl)
  lines <- strsplit(lines,"\\s+")
  close(fl)
  return(lines)
}

getQCRatios <- function(YAQCStatsObject) {
   vals <- YAQCStatsObject@gcos.probes
   unique.names <- rownames(vals)
   ##unique.names <- sub("[_-]?5_?.?_at$","",unique.names,perl=T)
   ##unique.names <- sub("[_-]?3_?.?_at$","",unique.names,perl=T)
   ##unique.names <- sub("[_-]?M_?.?_at$","",unique.names,perl=T)
   ##unique.names <- unique(unique.names)
   ##p3 <- .namegrep3(unique.names,rownames(vals))
   ##p5 <- .namegrep5(unique.names,rownames(vals))
   p3 <- grep("3",rownames(vals))
   p5 <- grep("5",rownames(vals))
   if (isLog(YAQCStatsObject)) { 
       res <- (vals[p3,] - vals[p5,])
       res <- 2^res
   }
   else { res <- (vals[p3,]/vals[p5,]) }
   rownames(res) <- paste(c("Actin","GAPDH"),".3'/5'",sep="")
   return(res)
}

## Removed in version > 1.7.0
## these 2 functions were in the older versions of simpleaffy
## but have disappeared in the latest one (2.13.z)
## as getQCRatios calls them, I have copied them from an older
## version of simpleaffy
##.namegrep3 <- function(stems,all) {
##  sapply(stems,function(stem) {
##    ## grep(paste(stem,"[-_wC]3.?_?.?_at$",sep=""),all,value=T)
##    grep(paste(stem,"[-_wC]?3_?.?_at$",sep=""),all,value=T)
##  });
##}
##.namegrep5 <- function(stems,all) {
##  sapply(stems,function(stem) {
##    ## grep(paste(stem,"[-_wC]5.?_?.?_at$",sep=""),all,value=T)
##    grep(paste(stem,"[-_wC]?5_?.?_at$",sep=""),all,value=T)
##  });
##}



getAllInt <- function(YAQCStatsObject,pattern) {
  vals <- YAQCStatsObject@morespikes
  unique.names <- rownames(vals)
  if (pattern=="biob") pattern <- "b[5|3|m]"
  if (pattern=="bioc") pattern <- "c[5|3]"
  if (pattern=="biod") pattern <- "d[5|3]"
  lst <- grep(pattern,unique.names,ignore.case=T,value=T)
  if (length(lst)<1) stop("No probe found with '",pattern,"' pattern.",sep="")
  int<-apply(YAQCStatsObject@morespikes[lst,],2,mean)
  return(int)
}

.bcalls<-function(YAQCStatsObject,pattern) {
    calls <- rownames(bioCalls(YAQCStatsObject))
    lst <- grep(pattern,calls,ignore.case=T,value=T)
    apply(YAQCStatsObject@bio.calls[lst,],2,function(x) {
	x[x!="P"] <- 0
        x[x=="P"] <- 1
        x<-as.numeric(x)
        return(100 * sum(x)/length(x))
    })
}


###############################################################################
#                   creation of an YAQCStats object                           #
###############################################################################

yaqc.affy<-function(object,      ## affybatch or expressionSet object
                    myYaqcControlProbes=NULL, ## YaqcControlProbes object
                    alphas=NULL, ## a numeric of size 2, with alpha1 and alpha2 values
                    tgt=100,     ## target value to scale arrays
                    tau=0.015,
                    logged=FALSE, ## is it in log2
                    verbose=FALSE) {
  ## If object is of class AffyBatch, then an ExpressionSet is calculated
  ## using the MAS5 algorithm (with simpleaffy's justMAS function).
  ## In that case, the scale factors (sfs) and target values (tgt)
  ## are set in the YAQCstats object and the final probe set intensities
  ## are stored unlogged to get the same values than Affymetrix's GCOS software.
  
  if (class(object)=="AffyBatch")
    yaqc.affy.affybatch(object,myYaqcControlProbes,alphas,tgt,tau,verbose=verbose)
  
  ## If there is no AffyBatch object, an ExpressionSet object must be
  ## provided. In this case, sfs, tgt, percentage present,
  ## average background, average noise, and probe calls are not computed
  ## and are not set in the YAQCStats object.
  
  else if (class(object)=="ExpressionSet")
    yaqc.affy.exprSet(object,myYaqcControlProbes,logged)
  
  ## Else, complain.
  
  else stop("Argument must be either an AffyBatch or ExpressionSet object.")
}

yaqc.affy.exprSet<-function(exprSetObject,myYaqcControlProbes,logged=logged) {
  if (is.null(myYaqcControlProbes)) {
    mybio <- c(getBio("b",c("5","3","m"),exprSetObject,onlyFirst=TRUE),
               getBio(c("c","d"),c("5","3"),exprSetObject,onlyFirst=TRUE))
    names(mybio) <- c("b5","b3","bm","c5","c3","d5","d3")
    myspk <- c(getSpk(c("dap","thr","lys","phe"),c("5","3","m"),exprSetObject,onlyFirst=TRUE))
    names(myspk) <-nms <- c("dap5","dap3","dapm",
                            "thr5","thr3","thrm",
                            "lys5","lys3","lysm",
                            "phe5","phe3","phem")
    mydeg <- c(getDeg(c("actin","gapdh"),c("5","3","m"),exprSetObject,onlyFirst=TRUE))
    names(mydeg) <- nms <- c("act5","act3","actm","gap5","gap3","gapm")
    myYaqcControlProbes <- new("YaqcControlProbes",
                               bio=new("YaqcBioProbes",bio=mybio),
                               spk=new("YaqcSpkProbes",spk=myspk),
                               deg=new("YaqcDegProbes",deg=mydeg),
                               info=paste("Generated automatically by getBio/Spk/Deg() in yaqcaffy version",
                                 package.version("yaqcaffy")))
  }
  qc.probenames <- deg(myYaqcControlProbes)
  spike.probenames <- c(bio(myYaqcControlProbes),spk(myYaqcControlProbes))
  qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) exprs(exprSetObject)[y,] )))
  spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) exprs(exprSetObject)[y,] )))
  ## returning a nice YAQCStats object
  return( new("YAQCStats",
              morespikes=t(spike.vals),
              gcos.probes=t(qc.probe.vals),
              log=logged,
              arraytype=exprSetObject@annotation,
              yaqcControlProbe=myYaqcControlProbes,
              objectVersion=package.version("yaqcaffy"))
         )
}


yaqc.affy.affybatch <- function(affyBatchObject,
                                myYaqcControlProbes,
                                alphas,
                                tgt,tau,
                                verbose=verbose) {
  ##cdfname<-cleancdfname(cdfName(affyBatchObject))
  ## creating simpleaffy's .qcEnv
  ##setQCEnvironment(cdfname)
  ## computing MAS5 values
  if (verbose) cat("MAS5 normalisation... ")
  eset.mas5 <- justMAS(affyBatchObject,tgt=tgt,scale=TRUE) 
  if (verbose) cat("done\n")
  ## getting names of the bio/spk/deg probes, either
  ## automatically with the getBio/getSpk/getDeg functions
  if (verbose) cat("Getting probe names... ")
  if (is.null(myYaqcControlProbes)) {
    ## create YaqcControlProbes object
    mybio <- c(getBio("b",c("5","3","m"),affyBatchObject,onlyFirst=TRUE),
               getBio(c("c","d"),c("5","3"),affyBatchObject,onlyFirst=TRUE))
    names(mybio) <- c("b5","b3","bm","c5","c3","d5","d3")
    myspk <- c(getSpk(c("dap","thr","lys","phe"),c("5","3","m"),affyBatchObject,onlyFirst=TRUE))
    names(myspk) <-nms <- c("dap5","dap3","dapm",
                            "thr5","thr3","thrm",
                            "lys5","lys3","lysm",
                            "phe5","phe3","phem")
    mydeg <- c(getDeg(c("actin","gapdh"),c("5","3","m"),affyBatchObject,onlyFirst=TRUE))
    names(mydeg) <- nms <- c("act5","act3","actm","gap5","gap3","gapm")
    myYaqcControlProbes <- new("YaqcControlProbes",
                               bio=new("YaqcBioProbes",bio=mybio),
                               spk=new("YaqcSpkProbes",spk=myspk),
                               deg=new("YaqcDegProbes",deg=mydeg),
                               info=paste("Generated automatically by getBio/Spk/Deg() in yaqcaffy version",
                                 package.version("yaqcaffy")))
  }
  qc.probenames <- deg(deg(myYaqcControlProbes))
  spike.probenames <- c(bio(bio(myYaqcControlProbes)),spk(spk(myYaqcControlProbes)))
  calls <- bio(bio(myYaqcControlProbes))
  if (verbose) cat("done\n")
  if (verbose) cat("Extracting data... ")
  sfs <- experimentData(eset.mas5)@preprocessing$sfs
  tgt <- experimentData(eset.mas5)@preprocessing$tgt 
  exprSetObject<-eset.mas5
  ## getting probe calls to compute percentage present
  ## using detection.p.value from the simpleaffy package
  if (is.null(alphas)) {
    ## getting alpha1 and alpha2 from
    ## simpleaffy's qcdef fille
    alpha1 <- getAlpha1(affyBatchObject)
    alpha2 <- getAlpha2(affyBatchObject)
    ## if no aqdef file found, setting default values
    ## and a warning is returned by getAlpha[1|2]
    if (is.null(alpha1)) alpha1 <- 0.04 
    if (is.null(alpha2)) alpha1 <- 0.06
  } else { 
    alpha1 <- alphas[1] ## user defined value for alpha1
    alpha2 <- alphas[2] ## user defined value for alpha2
  }
  det<-detection.p.val(affyBatchObject,
                       tau=tau,
                       alpha1=alpha1,
                       alpha2=alpha2)
  ppv<-apply(det$call,2,function(x) {
    x[x!="P"] <- 0
    x[x=="P"] <- 1
    x<-as.numeric(x)
    return(100 * sum(x)/length(x))
  })
  ## computing background and noise averages
  bgstats<-simpleaffy:::.bg.stats(affyBatchObject)
  meanbg <- apply(bgstats$zonebg,1,mean)
  meanns <- apply(bgstats$zonesd,1,mean)
  ## getting the spike and qc probes intensities
  qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) 2^(exprs(exprSetObject))[y,] )))
  spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) 2^(exprs(exprSetObject))[y,] )))
  ## getting probe calls for spike probes
  biocalls<-matrix(ncol=length(sampleNames(affyBatchObject)),nrow=length(calls))
  rownames(biocalls)<-sub("at","at_call",calls)
  colnames(biocalls)<-sampleNames(affyBatchObject)
  for (i in 1:length(calls)) {
    if(!is.na(calls[i])) { biocalls[i,] <- det$call[calls[i],] } 
  }
  if (verbose) cat("done\n")
  if (verbose) cat("Generation YAQCStats object...\n")
  ## returning a nice YAQCStats object
  return(new("YAQCStats",
             scale.factors=sfs,                   ## inherited from QCStats class
             percent.present=ppv,                 ## inherited from QCStats class
             average.background=meanbg,           ## inherited from QCStats class
             target=tgt,                          ## inherited from QCStats class
             arraytype=affyBatchObject@annotation,## inherited from QCStats class
             average.noise=meanns,
             morespikes=t(spike.vals),
             gcos.probes=t(qc.probe.vals),
             bio.calls=biocalls,
             log=FALSE,
             yaqcControlProbes=myYaqcControlProbes,
             objectVersion=package.version("yaqcaffy"))
         )
}

setMethod("yaqc","eSet",function(object,...) yaqc.affy(object,...))

###############################################################################
#                       quality control summary                               #
###############################################################################

setMethod("show","YAQCStats",function(object) yaqc.show(object))
setMethod("merge",c("YAQCStats","YAQCStats"),
                    function(x,y,...) merge.yaqc(x,y,...))


yaqc.show <- function(object) {
  samples<-colnames(gcosProbes(object))
  ## show a full YAQCStats object
  if(length(sfs(object))>0) {
    temp.table<-rbind(#object@target,
                      #object@arraytype,
                      object@scale.factors,
                      object@average.background,
                      object@average.noise,
                      object@percent.present,
                      object@morespikes,
                      object@gcos.probes,
                      object@bio.calls)
    rownames(temp.table)[1:4]<-c(#"target",
                                 #"array.type",
                                 "scale.factors",
                                 "average.background",
                                 "average.noise",
                                 "percent.present")
    colnames(temp.table)<-c(samples)
  }
  ## show a small YAQCStats object
  else {
    temp.table<-rbind(object@morespikes,
                      object@gcos.probes)
    colnames(temp.table)<-c(samples)
  }
  print(temp.table)
}


merge.yaqc <- function(o1,o2) {
  if (!identical(o1@yaqcControlProbes,o2@yaqcControlProbes)) {
    stop("YAQC Objects can not be merged because their have been generated using different contol probes\n")
  }
  if ((o1@log==o2@log) && (o1@target==o2@target) && (o1@arraytype==o2@arraytype) ) {
    meanns <- c(o1@average.noise,o2@average.noise)
    spike.vals <- matrix(,dim(o1@morespikes)[1],sum(dim(o1@morespikes)[2],dim(o2@morespikes)[2]))
    rownames(spike.vals) <- rownames(o1@morespikes)
    colnames(spike.vals) <- c(arrays(o1),arrays(o2))
    for (i in 1:dim(o1@morespikes)[1]) spike.vals[i,] <- c(o1@morespikes[i,],o2@morespikes[i,])
    
    gcos.probes.vals <- matrix(,dim(o1@gcos.probes)[1],sum(dim(o1@gcos.probes)[2],dim(o2@gcos.probes)[2]))
    rownames(gcos.probes.vals) <- rownames(o1@gcos.probes)
    colnames(gcos.probes.vals) <- c(arrays(o1),arrays(o2))
    for (i in 1:dim(o1@gcos.probes)[1]) gcos.probes.vals[i,] <- c(o1@gcos.probes[i,],o2@gcos.probes[i,])
    
    biocalls <- matrix(,dim(o1@bio.calls)[1],sum(dim(o1@bio.calls)[2],dim(o2@bio.calls)[2]))
    rownames(biocalls) <- rownames(o1@bio.calls)
    colnames(biocalls) <- c(arrays(o1),arrays(o2))
    for (i in 1:dim(o1@bio.calls)[1]) biocalls[i,] <- c(o1@bio.calls[i,],o2@bio.calls[i,])
    
    sfs <- c(o1@scale.factors,o2@scale.factors)
    ppv <- c(o1@percent.present,o2@percent.present)
    meanbg <- c(o1@average.background,o2@average.background)
    tgt <- o1@target
    log.vals <- o1@log
    atype <- o1@arraytype
    
    return(new("YAQCStats",
               scale.factors=sfs,
               percent.present=ppv,
               average.background=meanbg,
               target=tgt,
               arraytype=atype,
               average.noise=meanns,
               morespikes=as.matrix(spike.vals),
                gcos.probes=as.matrix(gcos.probes.vals),
               bio.calls=as.matrix(biocalls),
               log=log.vals,
               objectVersion=c(o1@objectVersion,o2@objectVersion),
               yaqcControlProbes=o1@yaqcControlProbes
               )
           )
  } else {
    stop("YAQC Objects can not be merged because of their (un)log values and/or target values
            or because they originate from different array types.\n")
  }
}


getOutliers <- function(YAQCStatsObject,slot) {
    if ( length(target(YAQCStatsObject))==0  ) {
         ## this YAQCStatsObject has been generated from an ExpressionSet input
         if (slot %in% c("pp","avns","avbg","sfs") )
             stop("This YAQCStats object has been generated using en ExpressionSet object
                   and does not have the 'pp', 'avns', 'avbg' and 'sfs' slots. Only 
                   'actin', 'gapdh', 'biob', 'bioc', 'biod', 'dap', 'phe' 'the' and 'lys'
                   have been computed.")
    }
         if (slot=="actin") data <- getQCRatios(YAQCStatsObject)[1,]
    else if (slot=="gapdh") data <- getQCRatios(YAQCStatsObject)[2,]
    else if (slot=="biob")  data <- getAllInt(YAQCStatsObject,"b[5|3|m]")
    else if (slot=="bioc")  data <- getAllInt(YAQCStatsObject,"c[5|3]")
    else if (slot=="biod")  data <- getAllInt(YAQCStatsObject,"d[5|3]")
    else if (slot=="dap")   data <- getAllInt(YAQCStatsObject,"dap")
    else if (slot=="thr")   data <- getAllInt(YAQCStatsObject,"thr")
    else if (slot=="phe")   data <- getAllInt(YAQCStatsObject,"phe")
    else if (slot=="lys")   data <- getAllInt(YAQCStatsObject,"lys")
    else if (slot=="pp")    data <- percent.present(YAQCStatsObject)
    else if (slot=="avbg")  data <- avbg(YAQCStatsObject)
    else if (slot=="avns")  data <- avns(YAQCStatsObject)
    else if (slot=="sfs") { data <- sfs(YAQCStatsObject)
                            ll <- mean(data)/2
                            ul <- mean(data)*1.5
                            smpl <- arrays(YAQCStatsObject)
                            nbSmpl <- length(smpl)
                            res <- numeric()
                            for (i in 1:nbSmpl) {
                              if (data[i]<ll || data[i]>ul) res[smpl[i]]<-data[i]
                            }
                            return(res)
                          }
    else { stop("Argument 'slot' not valid") }
    ll <- qnorm(.025,mean(data),sd(data))
    ul <- qnorm(.975,mean(data),sd(data))
    return(.removeOK(data,ll,ul))
}

.removeOK <- function(data,ll,ul) {
  res <- numeric()
  for ( i in 1:length(data) ) {
    if (data[i]<ll || ul<data[i]) res <- append(res,data[i])
  }
  return(res)
}

