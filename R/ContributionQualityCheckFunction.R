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
  ConObj1<-ContributionQualityCheckSubFunction(InternalSimObj1$InternalRankingInfo,TestCont)
  PercOut1<-ConObj1$Perc
  RelOut1<-ConObj1$Rel
  ConObj2<-ContributionQualityCheckSubFunction(InternalSimObj2$InternalRankingInfo,TestCont)
  PercOut2<-ConObj2$Perc
  RelOut2<-ConObj2$Rel
  PercOut<-MatrixMerge_Helper(PercOut1,PercOut2,na.rm=F)
  colnames(PercOut)<-c("SampleGroup1","SampleGroup2")
  RelOut<-MatrixMerge_Helper(RelOut1,RelOut2,na.rm=F)
  colnames(RelOut)<-c("SampleGroup1","SampleGroup2")
  TestTemp<-MatrixMerge_Helper(ConObj1$Test,ConObj2$Test)
  TestOut<-rowMeans(TestTemp)
  PassOut<-PercOut<=0.25
  return(list("PassIndicator"=PassOut,"PercentOverlaps"=PercOut,"RelativeDiff"=RelOut,"TestConMean"=TestOut))
}
