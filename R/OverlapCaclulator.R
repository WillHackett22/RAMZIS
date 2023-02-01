#' Overlap Calculator
#'
#' @param ReferenceDis reference distribution, generally higher
#' @param TestDis comparison distribution
#'
#' @return Overlap statistics
#' @export
#'
#' @examples #
OverlapCalculator<-function(ReferenceDis,TestDis){
  #Let the distributions be the needed input
  dR<-density(ReferenceDis,from=-0.1,to=1.1)
  dT<-density(TestDis,from=-.1,to=1.1,na.rm=T)
  if (max(TestDis)==0){
    dT$y<-rep(0,length(dT$y))
    dT$y[which(dT$x>=-0.001 & dT$x<=0.001)]<-1
  }
  sdT<-dT$y/sum(dT$y)
  sdR<-dR$y/sum(dR$y)
  #FP == Ref greater than Test
  FalsePosIDX<-which((sdT<sdR) & (sdT>0.000009))
  FP_T_ecdf<-ecdf(TestDis)
  #FN == Test greater than Ref
  FalseNegIDX<-which((sdT>=sdR) & (sdR>0.000009))
  FN_R_ecdf<-ecdf(ReferenceDis)

  FPcont<-which(diff(FalsePosIDX)!=1)
  FNcont<-which(diff(FalseNegIDX)!=1)
  if (length(FPcont)>0){
    internalFPidx<-c(1)
    tempidx<-2
    for (j in 1:length(FPcont)){
      internalFPidx[tempidx]<-FPcont[j]
      internalFPidx[tempidx+1]<-FPcont[j]+1
      tempidx<-tempidx+2
    }
    internalFPidx[tempidx]<-length(FalsePosIDX)
  } else {
    internalFPidx<-c(1,length(FalsePosIDX))
  }
  if (length(FNcont)>0){
    internalFNidx<-c(1)
    tempidx<-2
    for (j in 1:length(FNcont)){
      internalFNidx[tempidx]<-FNcont[j]
      internalFNidx[tempidx+1]<-FNcont[j]+1
      tempidx<-tempidx+2
    }
    internalFNidx[tempidx]<-length(FalseNegIDX)
  } else {
    internalFNidx<-c(1,length(FalseNegIDX))
  }
  FP<-0
  if (length(FalsePosIDX)>0){
    for (j in 1:(length(internalFPidx)/2)){
      jdxl<-(j-1)*2+1
      jdxh<-j*2
      fptemph<-FP_T_ecdf(dT$x[FalsePosIDX[internalFPidx[jdxh]]])
      fptempl<-FP_T_ecdf(dT$x[FalsePosIDX[internalFPidx[jdxl]]])
      FP<-FP+(fptemph-fptempl)
    }
  }

  FN<-0
  if (length(FalseNegIDX)>0){
    for (j in 1:(length(internalFNidx)/2)){
      jdxl<-(j-1)*2+1
      jdxh<-j*2
      fntemph<-FN_R_ecdf(dR$x[FalseNegIDX[internalFNidx[jdxh]]])
      fntempl<-FN_R_ecdf(dR$x[FalseNegIDX[internalFNidx[jdxl]]])
      FN<-FN+(fntemph-fntempl)
    }
  }
  Overlap<-FN+FP
  OverlapPerc<-Overlap/(2-Overlap)
  FP<-FP/(2-Overlap)
  FN<-FN/(2-Overlap)
  FP_ibounds<-FalsePosIDX[internalFPidx]
  FN_ibounds<-FalseNegIDX[internalFNidx]
  Output<-list('FP'=FP,'FN'=FN,'Overlap'=Overlap,'PercOverlap'=OverlapPerc,'FPi'=FP_ibounds,'FNi'=FN_ibounds)
  return(Output)

}
