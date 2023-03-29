#' RankingByZScoreV2 determines rank
#'
#' @param SimilarityObj Similarity object
#'
#' @return rankings
#' @export
#'
#' @examples #
RankingByZScoreV2<-function(SimilarityObj){
  ZData<-ZScoreContributionFunction(SimilarityObj$RankInfo$Null,SimilarityObj$RankInfo$Test)
  ZMeans<-rowMeans(ZData,na.rm=T)
  return(sort(ZMeans))
}
