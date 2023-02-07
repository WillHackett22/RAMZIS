#' RankingQualityFunction combines ranking and QC data and produces the final output
#'
#' @param SimilarityObj Similarity Object
#' @param InternalSimObj1 Internal Similarity Object 1
#' @param InternalSimObj2 Internal Similarity Object 2
#'
#' @return Quality and ranking object
#' @export
#'
#' @examples #
RankingQualityFunction<-function(SimilarityObj,InternalSimObj1,InternalSimObj2){
  CQData<-ContributionQualityCheckFunction(SimilarityObj,InternalSimObj1,InternalSimObj2)
  RData<-RankingByZScore(SimilarityObj)
  RankingDataSummary<-MatrixMerge_Helper(CQData,RData,na.rm=F)
  colnames(RankingDataSummary)<-c("QualityGroup1","QualityGroup2","Overlap%_1","Overlap%_2","RelDiff_1","RelDiff_2","TestConMean","ZScore")
  tempval1<-RankingDataSummary$QualityGroup1
  tempval1[is.na(tempval1)]<-FALSE
  tempval2<-RankingDataSummary$QualityGroup2
  tempval2[is.na(tempval2)]<-FALSE
  RankingDataSummary$PassOverall<-(tempval1&tempval2)
  gplist<-colnames(SimilarityObj$RankInfoActual)
  MinObs<-row.names(RankingDataSummary) %in% gplist
  RankingDataSummary$MinimumObs<-MinObs
  RankingDataSummary<-RankingDataSummary[order(RankingDataSummary$ZScore),]
  return(RankingDataSummary)
}
