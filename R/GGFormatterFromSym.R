#' GGFormatterFromSym takes similarity distributions and turns them in TIDY compatible density data
#'
#' @param TestDist Tanimoto similarity vector from the Test Comparison
#' @param RefDist Tanimoto similarity vector from the Reference Comparison
#'
#' @return TIDY formatted data
#' @export
#'
#' @examples #
GGFormatterFromSym<-function(TestDist,RefDist,RefName){
  TDis<-RelDensity(TestDist)
  RDis<-RelDensity(RefDist)
  OverlapData<-OverlapCalculator(RefDist,TestDist)
  dfTemp<-data.frame(matrix(0,nrow=length(TDis$x),ncol=5))
  colnames(dfTemp)<-c("Sim","Test","Ref","Beta","Alpha")
  dfTemp$Sim<-TDis$x
  dfTemp$Test<-TDis$y
  dfTemp$Ref<-RDis$y
  FPi<-AlphaIdxs(OverlapData)
  if (length(FPi)>1){
    dfTemp$Alpha[FPi]<-TDis$y[FPi]
  } else {
    if (!is.na(FPi)){
      dfTemp$Alpha[FPi]<-TDis$y[FPi]
    }
  }
  FNi<-BetaIdxs(OverlapData)
  if (length(FNi)>1){
    dfTemp$Beta[FNi]<-RDis$y[FNi]
  } else {
    if (!is.na(FNi)){
      dfTemp$Beta[FNi]<-RDis$y[FNi]
    }
  }
  #make into tidy format
  tl<-length(TDis$x)
  dfTidy<-data.frame(matrix(0,nrow=tl*4,ncol=3))
  colnames(dfTidy)<-c('Sim','Dens','Dist')
  dfTidy$Sim<-rep(dfTemp$Sim,4)
  colnames(dfTemp)<-c("Sim","Test",RefName,"Beta","Alpha")
  for (j in 1:4){
    dfTidy$Dens[((j-1)*tl+1):(j*tl)]<-dfTemp[,1+j]
    dfTidy$Dist[((j-1)*tl+1):(j*tl)]<-rep(colnames(dfTemp)[1+j],tl)
  }
  dfTidy$Dist<-factor(dfTidy$Dist,unique(dfTidy[order(dfTidy$Dens,decreasing=T),"Dist"]))
  return(dfTidy)
}

#' GGFormatterGeneral takes similarity distributions and turns them in TIDY compatible density data
#'
#' @param TestDist Tanimoto similarity vector from the Test Comparison
#' @param RefDist Tanimoto similarity vector from the Null Comparison
#'
#' @return TIDY formatted data
#' @export
#'
#' @examples #
GGFormatterGeneral<-function(TestDist,RefDist){
  dfTidy<-GGFormatterFromSym(TestDist,RefDist,"Null")
  return(dfTidy)
}

#' GGFormatterInternal takes similarity distributions and turns them in TIDY compatible density data
#'
#' @param TestDist Tanimoto similarity vector from the Test Comparison
#' @param RefDist Tanimoto similarity vector from the Internal Comparison
#' @param IntGrpNum A numeric (1 or 2) indicating which file the internal comparison is in relation to
#'
#' @return TIDY formatted data
#' @export
#'
#' @examples #
GGFormatterInternal<-function(TestDist,RefDist,IntGrpNum){
  dfTidy<-GGFormatterFromSym(TestDist,RefDist,paste0('Internal',IntGrpNum))
  return(dfTidy)
}

