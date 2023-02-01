#' ZScoreContributionFunction determines ZScore ranking between test and null
#'
#' @param ReferenceCont Null distribution
#' @param TestCont Test similarity distribution
#'
#' @return z scores for determining ranks
#' @export
#'
#' @examples #
ZScoreContributionFunction<-function(ReferenceCont,TestCont){
  ZOut<-data.frame(matrix(NA,nrow=dim(ReferenceCont)[1],ncol=dim(TestCont)[2]),row.names = row.names(ReferenceCont))
  ZN<-row.names(ZOut)
  for (j in 1:(dim(ZOut)[1])){
    ZOut[ZN[j],]<-ZScoreDistributionFunction(ReferenceCont[ZN[j],],TestCont[ZN[j],])
  }
  return(ZOut)
}
