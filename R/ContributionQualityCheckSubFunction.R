#' ContributionQualityCheckSubFunction used by CQCFunction
#'
#' @param InternalCont Internal Ranking info
#' @param TestCont Test Ranking info
#'
#' @return QC results
#'
#' @examples #
ContributionQualityCheckSubFunction<-function(InternalCont,TestCont){
  rl<-dim(TestCont)[2]
  il<-dim(InternalCont)[2]
  tempcont<-merge(TestCont,InternalCont,by=0,all=T)
  tempcont<-tempcont[,-1]
  row.names(tempcont)<-row.names(TestCont)
  iN<-row.names(InternalCont)
  tempcont[is.na(tempcont)]<-0
  tempout<-data.frame(matrix(NA,nrow=dim(TestCont)[1],ncol=1),row.names=row.names(TestCont))
  for (j in 1:dim(InternalCont)[1]){
    tempout[iN[j],]<-OverlapCalculator(unlist(tempcont[iN[j],(rl+1):(rl+il)]),unlist(tempcont[iN[j],1:rl]))$PercOverlap
  }
  return(tempout)
}
