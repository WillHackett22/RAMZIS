#' SimPlot plots the Test and Null distributions
#'
#' @param PlotTitle Plot Title string
#' @param SimilarityObj Test and Null Similarity Object
#' @param legendpos legend position
#' @param verbose Boolean to send overlap data to console
#' @param legendval Boolean to produce legend DEFAULT=TRUE
#' @param xbound xbound vector default c(0,1)
#' @param ybound ybound vector default goes to max peak height
#' @param logval Boolean indicating log transform before relative abundance. DEFAULT=TRUE
#' @param legendcexval legend text size
#'
#' @return Makes plot.
#' @export
#'
#' @examples #
SimPlot<-function(PlotTitle,SimilarityObj,legendpos='topleft',verbose=F,legendval=FALSE,xbound=c(0,1),ybound=T,logval=FALSE,legendcexval=1.1){
  k<-100
  #build density
  TDis<-SimilarityObj$Summary$Tanimoto
  dT<-density(TDis,from=-0.1,to=1.1,na.rm=T)
  sdTy<-k*dT$y/sum(dT$y)
  if (logval==T){
    sdTy<-log(sdTy+1)
    if (any(is.infinite(sdTy))){
      sdTy[which(is.infinite(sdTy))]<-0
    }
  }

  NDis<-SimilarityObj$NullOut$NullTani
  dN<-density(NDis,from=-0.1,to=1.1,na.rm=T)
  sdNy<-k*dN$y/sum(dN$y)
  if (logval==T){
    sdNy<-log(sdNy+1)
    if (any(is.infinite(sdNy))){
      sdNy[which(is.infinite(sdNy))]<-0
    }
  }
  TAct<-SimilarityObj$Actual

  OverlapData<-OverlapCalculator(SimilarityObj$NullOut$NullTani,SimilarityObj$Summary$Tanimoto)
  if (logval==T){
    #OverlapData<-OverlapCalculator(SimilarityObj$NullOut$NullTani,SimilarityObj$Summary$Tanimoto)
  }
  #plot densities
  if (max(TDis)==0){
    mh<-max(sdNy)
    if (ybound==T){
      yubound<-c(0,mh)
    } else {
      yubound<-c(0,ybound)
    }
    plot(c(0,0,0.05,0.05),c(0,mh,mh,0),xlim=xbound,ylim=yubound,main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution')
    polygon(c(0,0,0.05,0.05),c(0,mh,mh,0),col=rgb(1,0,0,0.5))
  } else {
    mh<-max(c(sdNy,sdTy))
    if (ybound==T){
      yubound<-c(0,mh)
    } else {
      yubound<-c(0,ybound)
    }
    plot(dT$x,sdTy,xlim=xbound,ylim=yubound,main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution')
    polygon(c(-0.1001,dT$x,1.1001),c(0,sdTy,0),col=rgb(1,0,0,0.5))
  }
  lines(dN$x,sdNy)
  polygon(c(-0.1001,dN$x,1.1001),c(0,sdNy,0),col=rgb(0,0,1,0.5))
  lines(rep(TAct,2),c(0,mh),col=1,lwd=3)


  #percentile location
  Overlap<-OverlapData$PercOverlap
  AlphaValue<-OverlapData$FP
  BetaValue<-OverlapData$FN
  if (legendval){
    legend(legendpos,cex=legendcexval,legend=c(paste('Observed Similarity of ',round(TAct,2),' in the',round(ecdf(TDis)(TAct),2)*100,'Percentile'),'Test Distribution','Null Distribution',paste0(round(AlphaValue,3)*100,'% False Positive Rate'),paste0(round(BetaValue,3)*100,'% False Negative Rate')),fill=c(NA,rgb(1,0,0,0.5),rgb(0,0,1,0.5),'black','darkgray'),lty=c(1,rep(NA,4)),density=c(0,NA,NA,NA,NA),border=c(NA,1,1,1,1))
  }

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
      polygon(c(xi,dN$x[fni],x0,x0),c(0,sdNy[fni],sdNy[fni[length(fni)]],0),col='darkgrey')
    }
  }

  if (verbose==T){
    return(list(Overlap,TAct,'FalsePositiveRate'=AlphaValue,'FalseNegativeRate'=BetaValue))
  }

}
