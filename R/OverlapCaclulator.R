#' Overlap Calculator
#'
#' @param ReferenceDis reference distribution, generally higher
#' @param TestDis comparison distribution
#'
#' @return Overlap statistics
#' @export
#'
#' @examples #
OverlapCalculator<-function(ReferenceDis,TestDis,zerohandling=F){
  #Let the distributions be the needed input
  dR<-DENSITY_RD(ReferenceDis)
  dT<-DENSITY_RD(TestDis)
  if (max(TestDis)==0){
    dT$y<-rep(0,length(dT$y))
    if (zerohandling){
      dT$y[c(43,44)]<-1/(dT$x[44]-dT$x[43])
    } else {
      dT$y[which(dT$x>=-0.001 & dT$x<=0.001)]<-1
    }
  }
  sdT<-dT$y/sum(dT$y)
  sdR<-dR$y/sum(dR$y)
  #FP == Ref greater than Test
  FalsePosIDX<-which((sdT<sdR) & (sdT>0.0001*max(sdT)))
  FP_T_ecdf<-ecdf(TestDis)
  #FN == Test greater than Ref
  FalseNegIDX<-which((sdT>=sdR) & (sdR>0.0001*max(sdR)))
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

  if ((max(TestDis)==0) | (max(ReferenceDis==0))){
    if (max(ReferenceDis==0) & max(TestDis==0)){
      FP<-0.5
      FN<-0.5
    } else if (max(TestDis==0) & any(ReferenceDis==0)){
      FN<-FN+sum(ReferenceDis==0)/length(ReferenceDis)
    } else if (max(ReferenceDis==0) & any(TestDis==0)){
      FP<-FP+sum(TestDis==0)/length(TestDis)
    } else if (max(TestDis==0) & all(ReferenceDis!=0)){
      FP<-0
      FN<-0
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
