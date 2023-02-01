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
  RankingDataSummary<-merge(CQData,RData,by=0,all=T)
  row.names(RankingDataSummary)<-RankingDataSummary[,1]
  RankingDataSummary<-RankingDataSummary[,-1]
  colnames(RankingDataSummary)<-c("QualityGroup1","QualityGroup2","Overlap%_1","Overlap%_2","ZScore")
  tempval1<-RankingDataSummary$QualityGroup1
  tempval1[is.na(tempval1)]<-TRUE
  tempval2<-RankingDataSummary$QualityGroup2
  tempval2[is.na(tempval2)]<-TRUE
  RankingDataSummary$PassOverall<-(tempval1&tempval2)
  RankingDataSummary<-RankingDataSummary[order(RankingDataSummary$ZScore),]
  return(RankingDataSummary)
}
