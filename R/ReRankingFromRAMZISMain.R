#' ReRankingFromRAMZISMain reformats the data object from RAMZISMain and sends it into RankingQualityFunctionV3. It is meant to allow easier reanalysis with changed ranking info
#'
#' @param SimilarityObject Output from RAMZISMain
#' @param RankingInfo Default="Numerator" is the standard ranking info, "WeightedContributions" changes to the weighted contribution as ranking info
#' @param QualityInfo Default="WeightedContributions" is the standard quality info, "Numerator" changes to the weighted contribution as quality info
#'
#' @return #RankingSummary
#' @export
#'
#' @examples #
ReRankingFromRAMZISMain<-function(SimilarityObject,RankingInfo="Numerator",QualityInfo="WeightedContributions"){
  TestSimObj<-list("Numerator"=SimilarityObject$RankInfo$RankData$Test,"WeightedContributions"=SimilarityObject$WeightedContributions$Test)
  NullSimObj<-list("Numerator"=SimilarityObject$RankInfo$RankData$Null,"WeightedContributions"=SimilarityObject$WeightedContributions$Null)
  InternalSimObj1<-list("Numerator"=SimilarityObject$RankInfo$RankData$Internal1,"WeightedContributions"=SimilarityObject$WeightedContributions$Internal1)
  InternalSimObj2<-list("Numerator"=SimilarityObject$RankInfo$RankData$Internal2,"WeightedContributions"=SimilarityObject$WeightedContributions$Internal2)
  ActualSimObj<-list("Contribution"=SimilarityObject$RankInfo$RankData$Actual)
  RankingDataSummary<-RankingQualityFunctionV3(TestSimObj = TestSimObj,NullSimObj = NullSimObj,ActualSimObj = ActualSimObj,InternalSimObj1 = InternalSimObj1,InternalSimObj2 = InternalSimObj2,QualityInfo=QualityInfo,RankingInfo=RankingInfo)
  return(RankingDataSummary)
}
