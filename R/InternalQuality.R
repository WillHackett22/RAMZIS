#' InternalQuality determines quality of Internal similarity
#'
#' @param filename abundance data of sample group
#' @param BootSet sampling of sample group
#' @param SimilarityObj Similarity Object from Test-Null Comparison
#' @param PlotTitle Title for Plot
#' @param GroupName Sample Group name
#' @param Int Internal Similarity object if already calculated. Default=NULL
#' @param kmin minimum number of observations needed to be considered real. Default=1
#' @param rel Boolean indicating relative abundance transformation
#' @param MVCorrection Boolean missing value correction, if disabled it turns 0s to NAs. Default=TRUE
#' @param mn Default=FALSE. In default settings it adjusts the distance scaling by average presence. Setting it to a numeric value will use that as a constant instead. (Recommended 1:2)
#' @param logopt Boolean indicating the use of the log transform before relative scaling of abundance. Default=TRUE
#' @param verbose Boolean indicating output. Default=FALSE. If TRUE sends Internal Similarity object to console
#' @param legendpos Legend position
#' @param legendval Legend Value BOOLEAN to produce Original comparison. Default=TRUE
#' @param legendcexval Legend text size
#' @param xbound Xbound vector default=c(0,1)
#' @param ybound ybound vector default has max at max peak height
#' @param logval log transform boolean DEFAULT=TRUE
#'
#' @return Plots Internal Similarity
#' @export
#'
#' @examples #
InternalQuality<-function(filename,BootSet,SimilarityObj,PlotTitle,GroupName,Int=NULL,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,verbose=FALSE,legendpos='topleft',legendval=FALSE,legendcexval=1.1,xbound=c(0,1),ybound=T,logval=T){
  if (!is.null(Int)){
    #Use Int as Int object
  } else {
    Int<-InternalSimilarity(filename,BootSet,kmin,rel,MVCorrection,mn)
  }

  k<-100
  #build density
  TDis<-SimilarityObj$Summary$Tanimoto
  TDis[is.infinite(TDis)]<-NA
  dT<-DENSITY_RD(TDis)
  sdTy<-k*dT$y/sum(dT$y)
  if (logval==T){
    sdTy<-log(sdTy+1)
    if (any(is.infinite(sdTy))){
      sdTy[which(is.infinite(sdTy))]<-0
    }
  }
  NDis<-SimilarityObj$NullOut$NullTani
  NDis[is.infinite(NDis)]<-NA
  dN<-density(NDis,from=-0.1,to=1.1,na.rm=T)
  sdNy<-k*dN$y/sum(dN$y)
  if (logval==T){
    sdNy<-log(sdNy+1)
    if (any(is.infinite(sdNy))){
      sdNy[which(is.infinite(sdNy))]<-0
    }
  }
  IDis<-Int$InternalTanimoto
  IDis[is.infinite(IDis)]<-NA
  dI<-density(IDis,from=-0.1,to=1.1,na.rm=T)
  sdIy<-k*dI$y/sum(dI$y)
  if (logval==T){
    sdIy<-log(sdIy+1)
    if (any(is.infinite(sdIy))){
      sdIy[which(is.infinite(sdIy))]<-0
    }
  }
  TAct<-SimilarityObj$Actual
  OverlapData<-OverlapCalculator(IDis,TDis)
  #plot densities
  if (max(TDis)==0){
    mh<-max(c(k*dI$y/sum(dI$y)))
    if (ybound==T){
      yubound<-c(0,mh)
    } else {
      yubound<-c(0,ybound)
    }
    plot(c(0,0,0.1,0.1),c(0,mh,mh,0),xlim=xbound,ylim=yubound,main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution')
    polygon(c(0,0,0.1,0.1),c(0,mh,mh,0),col=rgb(1,0,0,0.5))
  } else {
    mh<-max(c(k*dT$y/sum(dT$y),k*dI$y/sum(dI$y)))
    if (ybound==T){
      yubound<-c(0,mh)
    } else {
      yubound<-c(0,ybound)
    }
    plot(dT$x,k*dT$y/sum(dT$y),ylim=yubound,main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution',xlim=xbound)
    polygon(c(-0.1001,dT$x,1.1001),c(0,k*dT$y/sum(dT$y),0),col=rgb(1,0,0,0.5))
  }
  #lines(dN$x,k*dN$y/sum(dN$y))
  #polygon(dN$x,k*dN$y/sum(dN$y),col=rgb(0,0,1,0.5))
  lines(dI$x,k*dI$y/sum(dI$y))
  polygon(c(-0.1001,dI$x,1.1001),c(0,k*dI$y/sum(dI$y),0),col=rgb(0,1,0,0.5))
  lines(rep(TAct,2),c(0,mh),col=1,lwd=3)
  #percentile location
  CompPerc<-ecdf(TDis)(TAct)
  JointPerc<-ecdf(NDis)(TAct)

  Overlap<-OverlapData$PercOverlap
  AlphaValue<-round(OverlapData$FP,3)
  BetaValue<-round(OverlapData$FN,3)

  if (AlphaValue>(10^-10)){
    #make alpha area
    for (j in 1:(length(OverlapData$FPi)/2)){
      jdxh<-j*2
      jdxl<-(j-1)*2+1
      fpi<-seq(OverlapData$FPi[jdxl],OverlapData$FPi[jdxh])
      xi<-mean(c(dT$x[fpi[1]],dT$x[fpi[1]-1]),na.rm=T)
      x0<-mean(c(dT$x[fpi[length(fpi)]],dT$x[fpi[length(fpi)]+1]),na.rm=T)
      polygon(c(xi,xi,dT$x[fpi],x0,x0),c(0,sdTy[fpi][1],sdTy[fpi],sdTy[fpi[length(fpi)]],0),col='black')
    }
  }
  if (BetaValue>(10^-10)){
    #make beta area
    for (j in 1:(length(OverlapData$FNi)/2)){
      jdxh<-j*2
      jdxl<-(j-1)*2+1
      fni<-seq(OverlapData$FNi[jdxl],OverlapData$FNi[jdxh])
      xi<-mean(c(dT$x[fni[1]],dT$x[fni[1]-1]),na.rm=T)
      x0<-mean(c(dT$x[fni[length(fni)]],dT$x[fni[length(fni)]+1]),na.rm=T)
      polygon(c(xi,dN$x[fni],x0,x0),c(0,sdIy[fni],sdIy[fni[length(fni)]],0),col='darkgrey')
    }
  }

  deltaIT<-(mean(IDis,na.rm=T)-mean(TDis,na.rm = T))
  Confidence<-round(deltaIT/sd(IDis,na.rm=T)*10^(-AlphaValue-BetaValue),2)
  if (legendval){
    legend(legendpos,legend=c(paste('Observed Similarity of ',round(TAct,2),' in the',round(ecdf(TDis)(TAct),2)*100,'Percentile'),paste('Internal of',GroupName,'with Score=',Confidence),'Test Distribution',paste0(AlphaValue*100,'% False Positive Rate'),paste0(BetaValue*100,'% False Negative Rate')),fill=c(NA,3,2,'black','darkgrey'),lty=c(1,rep(NA,4)),lwd=c(3,NA,NA,NA,NA),density=c(0,NA,NA,NA,NA),border=c(NA,1,1,1,1),cex=legendcexval)
  }
  if (verbose==T){
    XPoint<-NA
    return(list(Int,Overlap,XPoint,CompPerc,JointPerc,Confidence,AlphaValue,BetaValue))
  }


}
