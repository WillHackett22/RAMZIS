#' ContributionQualityCheckFunction
#'
#' @param SimilarityObj TestnullSimilarity Object
#' @param InternalSimObj1 Internal Similarity object of first group
#' @param InternalSimObj2 Internal similarity object of second group
#'
#' @return Contribution QC
#' @export
#'
#' @examples #
ContributionQualityCheckFunction<-function(SimilarityObj,InternalSimObj1,InternalSimObj2){
  TestCont<-SimilarityObj$RankInfoFinal
  PercOut1<-ContributionQualityCheckSubFunction(InternalSimObj1$InternalRankingInfo,TestCont)
  PercOut2<-ContributionQualityCheckSubFunction(InternalSimObj2$InternalRankingInfo,TestCont)
  PercOut<-merge(PercOut1,PercOut2,by=0,all=T)
  row.names(PercOut)<-PercOut[,1]
  PercOut<-PercOut[,-1]
  colnames(PercOut)<-c("SampleGroup1","SampleGroup2")
  PassOut<-PercOut<=0.25
  return(list("PassIndicator"=PassOut,"PercentOverlaps"=PercOut))
}
