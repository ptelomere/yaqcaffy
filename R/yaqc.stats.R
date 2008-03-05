## accessor methods for the YAQCStats class
setMethod("avns",             "YAQCStats",function(object) object@average.noise)
setMethod("moreSpikeInProbes","YAQCStats",function(object) object@morespikes)
setMethod("gcosProbes",       "YAQCStats",function(object) object@gcos.probes)
setMethod("bioCalls",         "YAQCStats",function(object) object@bio.calls)
setMethod("arrays",           "YAQCStats",function(object) names(object@average.noise))
setMethod("isLog",            "YAQCStats",function(object) object@log)

## Functions similar to the ones implemented in the quality control 
## related file (qc.stats.R) of the simpleaffy package, but modifed 
## and renamed for the YAQCStats Object

getSpikeProbes <- function(object) {
  cdfn<-cleancdfname(annotation(object))
  r <- get("morespikes",envir=.yaqcEnv)[cdfn,]
  return(r[!is.na(r)])
}

getBioProbes <- function(affyBatchObject) {
  cdfn<-cleancdfname(cdfName(affyBatchObject))
  return(
         c(as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biob5"]),
           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biobm"]),
           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biob3"]),
           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2bioc5"]),
           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2bioc3"]),
           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biod5"]),
           as.character(get("morespikes",envir=.yaqcEnv)[cdfn,"r2biod3"]))
         )
}

getRatioProbes <- function(object) {
  qlist <- getQcdefData(object)
  actin <- c(qlist[[9]][3:4],qlist[[8]][4])
  gapdh <- c(qlist[[11]][3:4],qlist[[10]][4])
  return(c(actin,gapdh))
}

getAlpha1 <- function(object) {
  qlist <- getQcdefData(object)
  return(qlist[[2]][2])
}

getAlpha2 <- function(object) {
  qlist <- getQcdefData(object)
  return(qlist[[3]][2])
}

getQcdefData <- function(object) {
  ## this funciton retrieves a list of information
  ## about the specific type of affyBatchObject
  cdfn<-cleancdfname(annotation(object))
  fn <- paste(cdfn,"qcdef",sep=".")
  qcfile <- system.file("extdata",fn,package="simpleaffy")
  fl <- file(qcfile,"r")
  if(!isOpen(fl)) { stop(paste("Couldn't open file",qcfile,"\n.")) }
  lines <- readLines(fl)
  lines <- strsplit(lines,"\\s+")
  close(fl)
  return(lines)
}
  
getQCRatios <- function(YAQCStatsObject) {
   vals <- YAQCStatsObject@gcos.probes
   unique.names <- rownames(vals)
   unique.names <- sub("[_-]5.?_?.?_at$","",unique.names,perl=T)
   unique.names <- sub("[_-]3.?_?.?_at$","",unique.names,perl=T)
   unique.names <- sub("[_-]M.?_?.?_at$","",unique.names,perl=T)
   unique.names <- unique(unique.names)
   p3 <- .namegrep3(unique.names,rownames(vals))
   p5 <- .namegrep5(unique.names,rownames(vals))
   if (isLog(YAQCStatsObject)) { 
       res <- rbind(c(),(vals[p3,] - vals[p5,])) 
       res <- 2^res
   }
   else { res <- rbind(c(),(vals[p3,]/vals[p5,])) }
   rownames(res) <- paste(unique.names,".3'/5'",sep="")
   return(res)

}


## these 2 functions were in the older versions of simpleaffy
## but have disappeared in the latest one (2.13.z)
## as getQCRatios calls them, I have copied them from an older
## version of simpleaffy
.namegrep3 <- function(stems,all) {
  sapply(stems,function(stem) {
    grep(paste(stem,"[-_wC]3.?_?.?_at$",sep=""),all,value=T)
  });
}

.namegrep5 <- function(stems,all) {
  sapply(stems,function(stem) {
    grep(paste(stem,"[-_wC]5.?_?.?_at$",sep=""),all,value=T)
  });
}



getAllInt <- function(YAQCStatsObject,pattern) {
    vals <- YAQCStatsObject@morespikes
    unique.names <- rownames(vals)
    lst <- grep(pattern,unique.names,ignore.case=T,value=T)
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
#                       calculation GCOS values                               #
###############################################################################

.gcos <- function(affyBatchObject,tgt=100) {
    mas<-justMAS(affyBatchObject,tgt=tgt,scale=TRUE)
    return(mas)
}

###############################################################################
#                   creation of an YAQCStats object                           #
###############################################################################

yaqc.affy<-function(object,      ## affybatch or expressionSet object
                    tgt=100,     ## target value to scale arrays
                    tau=0.015,
                    logged=FALSE ## is it in log2
                    ) {
    ## If object is of class AffyBatch, then an ExpressionSet is calculated
    ## using the MAS5 algorithm (with simpleaffy's justMAS function).
    ## In that case, the scale factors (sfs) and target values (tgt)
    ## are set in the YAQCstats object and the final probe set intensities
    ## are stored unlogged to get the same values than Affymetrix's GCOS software.
  
    if (class(object)=="AffyBatch") yaqc.affy.affybatch(object,tgt,tau)
  
    ## If there is no AffyBatch object, an expressionSet object must be
    ## provided. In this case, sfs, tgt, percentage present,
    ## average background, average noise, and probe calls are not computed
    ## and are not set in the YAQCStats object.
    
    else if (class(object)=="ExpressionSet") yaqc.affy.exprSet(object,logged)

    ## Else, complain.
    
    else stop("Argument must be either an AffyBatch or ExpressionSet object.")
}
  
yaqc.affy.exprSet<-function(exprSetObject, logged) {
    qc.probenames <- getRatioProbes(exprSetObject)
    qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) exprs(exprSetObject)[y,] )))
    spike.probenames <- getSpikeProbes(exprSetObject)
    spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) exprs(exprSetObject)[y,] )))
    ## returning a nice YAQCStats object
    return( new("YAQCStats",
        morespikes=t(spike.vals),
        gcos.probes=t(qc.probe.vals),
        log=logged,
        arraytype=exprSetObject@annotation
        )
    )
}


yaqc.affy.affybatch<-function(affyBatchObject,tgt,tau) {
    cdfname<-cleancdfname(cdfName(affyBatchObject))
    ## creating simpleaffy's .qcEnv
    setQCEnvironment(cdfname)
    ## computing gcos values as in gcos function (see above)
    eset.mas5<-.gcos(affyBatchObject,tgt=tgt)
    sfs <- experimentData(eset.mas5)@preprocessing$sfs
    tgt <- experimentData(eset.mas5)@preprocessing$tgt 
    exprSetObject<-eset.mas5
    ## getting probe calls to compute percentage present
    ## using detection.p.value from the simpleaffy package
    det<-detection.p.val(affyBatchObject,
                         tau=tau,
                         alpha1=getAlpha1(affyBatchObject),
                         alpha2=getAlpha2(affyBatchObject))
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
    qc.probenames <- getRatioProbes(affyBatchObject)
    qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) 2^(exprs(exprSetObject))[y,] )))
    spike.probenames <- getSpikeProbes(affyBatchObject)
    spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) 2^(exprs(exprSetObject))[y,] )))
    ## getting probe calls for spike probes
    calls <- getBioProbes(affyBatchObject)
    biocalls<-matrix(ncol=length(sampleNames(affyBatchObject)),nrow=length(calls))
    rownames(biocalls)<-sub("at","at_call",calls)
    colnames(biocalls)<-sampleNames(affyBatchObject)
    for (i in 1:length(calls)) {
        if(!is.na(calls[i])) { biocalls[i,] <- det$call[calls[i],] } 
    }
    ## returning a nice YAQCStats object
    return( new("YAQCStats",
        scale.factors=sfs,                   ## inherited from QCStats class
        percent.present=ppv,                 ## inherited from QCStats class
        average.background=meanbg,           ## inherited from QCStats class
        target=tgt,                          ## inherited from QCStats class
        arraytype=affyBatchObject@annotation,## inherited from QCStats class
        average.noise=meanns,
        morespikes=t(spike.vals),
        gcos.probes=t(qc.probe.vals),
        bio.calls=biocalls,
        log=FALSE
        )
    )
}


setMethod("yaqc","eSet",function(object,...) yaqc.affy(object,...))

###############################################################################
#                       quality control summary                               #
###############################################################################

setMethod("show","YAQCStats",function(object) yaqc.show(object))
setMethod("merge","YAQCStats",function(x,y,...) merge.yaqc(x,y,...))


yaqc.show <- function(object) {
  samples<-colnames(gcosProbes(object))
  ## show a full YAQCStats object
  if(length(sfs(object))>0) {
    temp.table<-rbind(
                      #object@target,
                      #object@arraytype,
                      object@scale.factors,
                      object@average.background,
                      object@average.noise,
                      object@percent.present,
                      object@morespikes,
                      object@gcos.probes,
                      object@bio.calls)
    rownames(temp.table)[1:4]<-c(
                                 #"target",
                                 #"array.type",
                                 "scale.factors",
                                 "average.background",
                                 "average.noise",
                                 "percent.present")
    colnames(temp.table)<-c(samples)
  }
  ## show a small YAQCStats object
  else {
    temp.table<-rbind(
                      object@morespikes,
                      object@gcos.probes)
    colnames(temp.table)<-c(samples)
  }
  print(temp.table)
}


merge.yaqc <- function(o1,o2) {
  if ((o1@log == o2@log) && (o1@target==o2@target) && (o1@arraytype==o2@arraytype) ) {
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
    tgt <- c(o1@target,o2@target)
    log.vals <- o1@log
    atype <- c(o1@arraytype,o2@arraytype)

    return( new("YAQCStats",
                scale.factors=sfs,
                percent.present=ppv,
                average.background=meanbg,
                target=tgt,
                arraytype=atype,
                average.noise=meanns,
                morespikes=as.matrix(spike.vals),
                gcos.probes=as.matrix(gcos.probes.vals),
                bio.calls=as.matrix(biocalls),
                log=log.vals
                )
           )
  }
  else cat("YAQC Objects can not be merged because of their (un)log values and/or target values
            or because they originate from different array types.\n")
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
         if (slot=="actin") data<-getQCRatios(YAQCStatsObject)[1,]
    else if (slot=="gapdh") data<-getQCRatios(YAQCStatsObject)[2,]
    else if (slot=="biob")  data<-getAllInt(YAQCStatsObject,"biob")
    else if (slot=="bioc")  data<-getAllInt(YAQCStatsObject,"bioc")
    else if (slot=="biod")  data<-getAllInt(YAQCStatsObject,"biod")
    else if (slot=="dap")   data<-getAllInt(YAQCStatsObject,"dap")
    else if (slot=="thr")   data<-getAllInt(YAQCStatsObject,"thr")
    else if (slot=="phe")   data<-getAllInt(YAQCStatsObject,"phe")
    else if (slot=="lys")   data<-getAllInt(YAQCStatsObject,"lys")
    else if (slot=="pp")    data<-percent.present(YAQCStatsObject)
    else if (slot=="avbg")  data<-avbg(YAQCStatsObject)
    else if (slot=="avns")  data<-avns(YAQCStatsObject)
    else if (slot=="sfs") { data <- sfs(YAQCStatsObject)
                            ll<-mean(data)/2
                            ul<-mean(data)*1.5
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
  res<-c();
  for ( i in 1:length(data) ) {
    if (data[i]<ll || ul<data[i]) res <- append(res,data[i])
  }
  return(res)
}

