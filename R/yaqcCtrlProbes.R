## =================
## accessor methods
## - - - - - - - - -
setMethod("bio","YaqcControlProbes", function(object) object@bio)
setMethod("bio","YaqcBioProbes",     function(object) object@bio)
setMethod("spk","YaqcControlProbes", function(object) object@spk)
setMethod("spk","YaqcSpkProbes",     function(object) object@spk)
setMethod("deg","YaqcControlProbes", function(object) object@deg)
setMethod("deg","YaqcDegProbes",     function(object) object@deg)
setMethod("info","YaqcControlProbes",function(object) object@info)

setReplaceMethod("info", signature(object="YaqcControlProbes"),
                 function(object, value){
                   object@info <- value
                   return(object)
                 })

## ===================
## show methods
## - - - - - - - - - -
setMethod("show","YaqcControlProbes",function(object) {
  cat("Hybridization probes:\n")
  show(bio(object))
  cat("Labeling probes:\n")
  show(spk(object))
  cat("Degradation probes:\n")
  show(deg(object))
  cat(info(object),"\n")
})

setMethod("show","YaqcSpkProbes",function(object) {
  output <- c(paste("Dap[5,3,m]: '",object@spk["dap5"],"', '",object@spk["dap3"],"', '",object@spk["dapm"],"'",sep=""),
              paste("Thr[5,3,m]: '",object@spk["thr5"],"', '",object@spk["thr3"],"', '",object@spk["thrm"],"'",sep=""),
              paste("Phe[5,3,m]: '",object@spk["phe5"],"', '",object@spk["phe3"],"', '",object@spk["phem"],"'",sep=""),
              paste("Lys[5,3,m]: '",object@spk["lys5"],"', '",object@spk["lys3"],"', '",object@spk["lysm"],"'",sep=""))
  print(output)
})

setMethod("show","YaqcDegProbes",function(object) {
  output <- c(paste("Actin[5,3,m]: '",object@deg["act5"],"', '",object@deg["act3"],"', '",object@deg["actm"],"'",sep=""),
              paste("Gapdh[5,3,m]: '",object@deg["gap5"],"', '",object@deg["gap3"],"', '",object@deg["gapm"],"'",sep=""))
  print(output)
})

          
setMethod("show","YaqcBioProbes",function(object) {
  output <- c(paste("BioB[5,3,m]: '",object@bio["b5"],"', '",object@bio["b3"],"', '",object@bio["bm"],"'",sep=""),
              paste("BioC[5,3]:   '",object@bio["c5"],"', '",object@bio["c3"],"'",sep=""),
              paste("BioD[5,3]:   '",object@bio["d5"],"', '",object@bio["d3"],"'",sep=""))
  print(output)
})

## =================
## other
## - - - - - - - - -
setMethod("identical","YaqcControlProbes", function(x,y) 
          return(identical(bio(x),bio(y)) & identical(spk(x),spk(y)) & identical(deg(x),deg(y)))
          )

