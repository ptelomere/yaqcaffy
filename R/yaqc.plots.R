###############################################################################
#                       quality control plot                                  #
###############################################################################

##====================================================##
## helper functions                                   ##
##----------------------------------------------------##

.plotrectangle <- function(data,where) {
    x1<-where-.1
    x2<-where+.1
    m<-mean(data)
    s<-sd(data)
    blq <- qnorm(.025,m,s)
    buq <- qnorm(.975,m,s)
    rect(x1,blq,x2,buq,border="grey95")    
    segments(x1,m,x2,m,col="grey95")
}

.plotdata<-function(data,title,color="white",outliers=TRUE,...) {
    i<-1
    bx<-boxplot(data,main=title,col=color,cex=.5,...)
    if (outliers) {
        for(i in 1:length(bx$out)) { 
            text(1,bx$out[i],names(bx$out[i]),cex=0.7,pos=1) 
        }
    }
    abline(h=mean(data),col="red",lty=5)
    abline(h=qnorm(.025,mean(data),sd(data)),lty=3,col="red")
    abline(h=qnorm(.975,mean(data),sd(data)),lty=3,col="red")
    mtext(paste("CV:",round(sd(data)/mean(data),2)),3,cex=.5)    
}

.plotQCRatios <- function(YAQCStatsObject,
                          which=c("all","gapdh","actin"),
                          ...) {
  which <- match.arg(which)
  if (which=="all") {
    .plotQCRatios(YAQCStatsObject,which="actin",...)
    .plotQCRatios(YAQCStatsObject,which="gapdh",...)
  }
  if (which=="actin") {
    ## plotting beta-actin 3'/5' ratios
    data<-getQCRatios(YAQCStatsObject)[1,]
    .plotdata(data,title="beta-actin 3'/5'",...)
  }
  if (which=="gapdh") {
    ## plotting GAPDH 3'/5' ratios
    data<-getQCRatios(YAQCStatsObject)[2,]
    .plotdata(data,title="GAPDH 3'/5'",...)
  }
}

.plotSfs <- function(YAQCStatsObject,...) {
  ## plotting scale factors results
  data<-sfs(YAQCStatsObject)
  min<-mean(data)/2
  max<-mean(data)*1.5
  dotchart(data,xlim=c(min*0.1,max*1.1),main="scale\nfactors",...)
  abline(v=max,col="red")
  abline(v=min,col="red")
  mtext(paste("CV:",round(sd(data)/mean(data),2)),3,cex=.5)
}

.plotAvbg <- function(YAQCStatsObject,...) {
  ## plotting average background
  data<-YAQCStatsObject@average.background
  .plotdata(data,title="average\nbackground",...)
}

.plotAvns <- function(YAQCStatsObject,...) {
  ## plotting average noise
  data<-YAQCStatsObject@average.noise
  .plotdata(data,title="average\nnoise",...)
}

.plotPp <- function(YAQCStatsObject,...) {
  ## plotting percent present
  data<-YAQCStatsObject@percent.present
  .plotdata(data,title="% present",...)
}
  
.plotBioQC<-function(YAQCStatsObject,calls=TRUE,...) {
    ## plotting BioB spike values
    bb<-getAllInt(YAQCStatsObject,"biob")
    ## plotting BioC spike values
    bc<-getAllInt(YAQCStatsObject,"bioc")
    ## plotting BioD spike values
    bd<-getAllInt(YAQCStatsObject,"biod")
    ## BioB, BioC and BioD boxplots
    low<-min(bb,bc,bd)*0.9
    up <-max(bb,bc,bd)*1.1
    i<-1
    ## get the graph
    bx<-boxplot(bb,col="red",boxwex=0.3,at=0.7,ylim=c(low,up),...)
    ## firest rectangle
    .plotrectangle(bb,0.7)
    ## first boxplot (again)
    bx<-boxplot(bb,col="red",boxwex=0.3,at=0.7,ylim=c(low,up),add=T,...)
    for(i in 1:length(bx)) { text(0.7,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    ## second rectange
    .plotrectangle(bc,1)
    ## second boxplot
    bx<-boxplot(bc,col="blue",boxwex=0.3,at=1,add=T,...)
    for(i in 1:length(bx)) { text(1,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    ## third rectangle
    .plotrectangle(bd,1.3)
    ## third boxplot
    bx<-boxplot(bd,col="green",boxwex=0.3,at=1.3,add=T,...)
    for(i in 1:length(bx)) { text(1.3,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    if (calls) {
        b<-paste("BioB -",round(mean(.bcalls(YAQCStatsObject,"B")),0),"% present")
        c<-paste("BioC -",round(mean(.bcalls(YAQCStatsObject,"C")),0),"% present")
        d<-paste("BioD -",round(mean(.bcalls(YAQCStatsObject,"D")),0),"% present")
    }
    else {
        b<-"BioB"
        c<-"BioC"
        d<-"BioD"
    }
    b <- paste(b,"- CV:",round(sd(bb)/mean(bb),2))
    c <- paste(c,"- CV:",round(sd(bc)/mean(bc),2))
    d <- paste(d,"- CV:",round(sd(bd)/mean(bd),2))    
    legend("topleft",c(b,c,d),fill=c("red","blue","green"),cex=1)
  }

.plotSpikesQC<-function(YAQCStatsObject,...) {
    ## plotting lys values
    lys<-getAllInt(YAQCStatsObject,"lys")
    ## plotting dap values
    dap<-getAllInt(YAQCStatsObject,"dap")
    ## plotting phe values
    phe<-getAllInt(YAQCStatsObject,"phe")
    ## plotting thr values 
    thr<-getAllInt(YAQCStatsObject,"thr")
    ## lys, dap, phe and thr boxplots
    low<-min(lys,dap,phe,thr)*0.9
    up <-max(lys,dap,phe,thr)*1.1
    i<-0
    bx<-boxplot(dap,col="purple",boxwex=0.25,at=0.6,ylim=c(low,up),...)
    .plotrectangle(dap,0.6)
    bx<-boxplot(dap,col="purple",boxwex=0.25,at=0.6,ylim=c(low,up),add=T,...)
    for(i in 1:length(bx)) { text(0.6,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    .plotrectangle(thr,0.85)    
    bx<-boxplot(thr,col="maroon",boxwex=0.25,at=0.85,add=T,...)
    for(i in 1:length(bx)) { text(0.85,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    .plotrectangle(phe,1.15)
    bx<-boxplot(phe,col="yellow",boxwex=0.25,at=1.15,add=T,...)
    for(i in 1:length(bx)) { text(1.15,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    .plotrectangle(lys,1.4)    
    bx<-boxplot(lys,col="orange",boxwex=0.25,at=1.4,add=T,...)
    for(i in 1:length(bx)) { text(1.4,bx$out[i],names(bx$out[i]),cex=0.75,pos=1) }
    dd <- paste("Dap - CV:", round(sd(dap)/mean(dap),2))
    tt <- paste("Thr - CV:", round(sd(thr)/mean(thr),2))
    pp <- paste("Phe - CV:", round(sd(phe)/mean(phe),2))
    ll <- paste("Lys - CV:", round(sd(lys)/mean(lys),2))
    legend("topright",c(dd,tt,pp,ll),fill=c("purple","maroon","yellow","orange"),cex=1)
  }

##====================================================##
## main plotting functions                            ##
##----------------------------------------------------##

yaqc.plot<-function(YAQCStatsObject,
                    which=c("all","sfs","avbg","avns","pp","gapdh","actin","bio","spikes"),
                    ...) {
  which <- match.arg(which)
  if (which=="all") {
    if( length(sfs(YAQCStatsObject))>0 ) yaqc.plot.all(YAQCStatsObject,...)
    else yaqc.plot.part(YAQCStatsObject,...)
  }
  if (which=="sfs")   .plotSfs(YAQCStatsObject,...)
  if (which=="avbg")  .plotAvbg(YAQCStatsObject,...)
  if (which=="avns")  .plotAvns(YAQCStatsObject,...)
  if (which=="pp")    .plotPp(YAQCStatsObject,...)
  if (which=="actin") .plotQCRatios(YAQCStatsObject,"actin",...)
  if (which=="gapdh") .plotQCRatios(YAQCStatsObject,"gapdh",...)
  if (which=="bio")   .plotBioQC(YAQCStatsObject,calls=TRUE,...)
  if (which=="spikes").plotSpikesQC(YAQCStatsObject,...)
}
  
yaqc.plot.part<-function(YAQCStatsObject) {
    layout(matrix(c(1:2,3,4,3,4),2))
    ## plotting beta-actin and GAPDH 3'/5' ratios
    .plotQCRatios(YAQCStatsObject)
    ## plotting Bio QC metrics
    .plotBioQC(YAQCStatsObject,calls=FALSE)
    ## plotting spike QC metrics
    .plotSpikesQC(YAQCStatsObject)
}

yaqc.plot.all<-function(YAQCStatsObject) {
    op<-par(no.readonly=T) ## save original parameters
    par(mar=c(2,4,4,2))
    layout(matrix(c(1:6,rep(7,3),rep(8,3)),2,byrow=T))
    ## plotting AffyBatch object specific QC metrics
    .plotSfs(YAQCStatsObject,labels="")
    .plotAvbg(YAQCStatsObject)
    .plotAvns(YAQCStatsObject)
    .plotPp(YAQCStatsObject)
    ## plotting beta-actin and GAPDH 3'/5' ratios
    .plotQCRatios(YAQCStatsObject)
    par(mar=c(5,4,4,1))
    ## plotting Bio QC metrics
    .plotBioQC(YAQCStatsObject,calls=TRUE)
    ## plotting spike QC metrics
    .plotSpikesQC(YAQCStatsObject)
    par(op) ## restore original parameters
}

## plot.YAQCStats
## <- function(YAQCStatsObject,...) yaqc.plot(YAQCStatsObject,...)

setMethod("plot",c("YAQCStats","missing"),
          function(x,y,...) yaqc.plot(x,...)
          )

