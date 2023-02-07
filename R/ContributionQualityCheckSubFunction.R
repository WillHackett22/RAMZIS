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
  tempcont<-MatrixMerge_Helper(TestCont,InternalCont,na.rm=F)
  iN<-row.names(tempcont)
  tempout<-data.frame(matrix(NA,nrow=dim(tempcont)[1],ncol=1),row.names=row.names(tempcont))
  tempout2<-data.frame(matrix(NA,nrow=dim(tempcont)[1],ncol=1),row.names=row.names(tempcont))
  tempout3<-data.frame(matrix(NA,nrow=dim(tempcont)[1],ncol=1),row.names=row.names(tempcont))
  for (j in 1:dim(tempout)[1]){
    intval<-mean(unlist(tempcont[iN[j],(rl+1):(rl+il)]),na.rm=T)
    testval<-mean(unlist(tempcont[iN[j],1:rl]),na.rm=T)
    if ( (!is.na(testval))){
      if (!is.na(intval)){
        tempout[iN[j],]<-OverlapCalculator(unlist(tempcont[iN[j],(rl+1):(rl+il)]),unlist(tempcont[iN[j],1:rl]))$PercOverlap
        relval<-intval-testval
        tempout2[iN[j],]<-relval
        tempout3[iN[j],]<-testval
      } else {
        tempout[iN[j],]<-1
        relval<-0-testval
        tempout2[iN[j],]<-relval
        tempout3[iN[j],]<-testval
      }
    } else if ( (is.na(testval)) & (!is.na(intval))){
      tempout[iN[j],]<-NA
      tempout2[iN[j],]<-intval
      tempout3[iN[j],]<-NA
    } else {
      tempout[iN[j],]<-NA
      tempout2[iN[j],]<-NA
      tempout3[iN[j],]<-NA
    }
  }
  return(list("Perc"=tempout,"Rel"=tempout2,"Test"=tempout3))
}
