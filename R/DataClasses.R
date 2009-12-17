## =======================================
## YaqcControlProbes and related classes
## - - - - - - - - - - - - - - - - - - - -
setClass("YaqcBioProbes",
         representation(bio="character"),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           nms <- c("b5","b3","bm",
                    "c5","c3",
                    "d5","d3")
           l <- length(nms)
           if (length(object@bio)!=l)
             msg <- validMsg(msg, paste("Number of bio probes != ",l,sep=""))
           if (is.null(names(object@bio))) {
             msg <- validMsg(msg, "Bio probes are unnamed.")
             if (!all(names(object@bio) %in% nms))
               msg <- validMsg(msg, "Invalid bio probes names.")
           }
           if (is.null(msg)) TRUE
           else msg
         })

setClass("YaqcSpkProbes",
         representation(spk="character"),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           nms <- c("dap5","dap3","dapm",
                    "thr5","thr3","them",
                    "lys5","lys3","lysm",
                    "phe5","phe3","phem")
           l <- length(nms)
           if (length(object@spk)!=l)
             msg <- validMsg(msg, paste("Number of spk probes != ",l,sep=""))
           if (is.null(names(object@spk))) {
             msg <- validMsg(msg, "Spk probes are unnamed.")
             if (!all(names(object@spk) %in% nms))
               msg <- validMsg(msg, "Invalid spk probes names.")
           }
           if (is.null(msg)) TRUE
           else msg
         })


setClass("YaqcDegProbes",
         representation(deg="character"),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           nms <- c("act5","act3","actm",
                    "gap5","gap3","gapm")
           l <- length(nms)
           if (length(object@deg)!=l)
             msg <- validMsg(msg, paste("Number of deg probes != ",l,sep=""))
           if (is.null(names(object@deg))) {
             msg <- validMsg(msg, "Deg probes are unnamed.")
             if (!all(names(object@deg) %in% nms))
               msg <- validMsg(msg, "Invalid deg probes names.")
           }
           if (is.null(msg)) TRUE
           else msg
         })

setClass("YaqcControlProbes",
         representation(bio="YaqcBioProbes",
                        spk="YaqcSpkProbes",
                        deg="YaqcDegProbes",
                        info="character"),
         validity = function(object) {
           validObject(object@bio) & validObject(object@spk) & validObject(object@deg) 
         })


# ==========================================================================
# YAQCStats class: holds the data necessary for the quality control
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setClass("YAQCStats", 
         contains="QCStats",
         representation(average.noise="numeric",
                        morespikes="matrix",
                        gcos.probes="matrix",
                        bio.calls="matrix",
                        log="logical",
                        objectVersion="character", ## new in versions > 1.7.0 - library version
                        yaqcControlProbes="YaqcControlProbes"), ## new in versions > 1.7.0
         validity=function(object) {
           ## new in versions > 1.7.0
           ## to be updated
           ## things to add: (1) check that probes in yaqcControlProbes
           ## correspond to those in the moresikes, bio.calls and
           ## gcos.probes splots
           ## Note: YaqcContolProbes class has validity method
           return(TRUE)
         })


