#' RankingQualityFunctionV3 combines ranking and QC data and produces the final output
#'
#' @param TestSimObj Test Similarity Object
#' @param NullSimObj Null Similarity Object
#' @param ActualSimObj Observed/Actual Similarity Object
#' @param InternalSimObj1 Internal Similarity Object 1
#' @param InternalSimObj2 Internal Similarity Object 2
#' @param QualityInfo Default="WeightedContributions"
#' @param RankingInfo Default="Numerator"
#'
#' @return Quality and ranking object
#' @export
#'
#' @examples #
RankingQualityFunctionV3<-function(TestSimObj,NullSimObj,ActualSimObj,InternalSimObj1,InternalSimObj2,QualityInfo="WeightedContributions",RankingInfo="Numerator"){
  CQData<-ContributionQualityCheckFunctionV3(TestSimObj,InternalSimObj1,InternalSimObj2,QualityInfo)
  RData<-RankingByZScoreV3(TestSimObj,NullSimObj,RankingInfo)
  RankingDataSummary<-MatrixMerge_Helper(CQData,RData,na.rm=F)
  colnames(RankingDataSummary)<-c("QualityGroup1","QualityGroup2","Overlap%_1","Overlap%_2","RelDiff_1","RelDiff_2","TestConMean","ZScore")
  tempval1<-RankingDataSummary$QualityGroup1
  tempval1[is.na(tempval1)]<-FALSE
  tempval2<-RankingDataSummary$QualityGroup2
  tempval2[is.na(tempval2)]<-FALSE
  RankingDataSummary$PassOverall<-(tempval1&tempval2)
  gplist<-colnames(ActualSimObj$Contribution)
  MinObs<-row.names(RankingDataSummary) %in% gplist
  RankingDataSummary$MinimumObs<-MinObs
  RankingDataSummary<-RankingDataSummary[order(RankingDataSummary$ZScore),]
  return(RankingDataSummary)
}
