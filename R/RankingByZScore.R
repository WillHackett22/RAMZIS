#' RankingByZScore determines rank
#'
#' @param SimilarityObj Similarity object
#'
#' @return rankings
#' @export
#'
#' @examples #
RankingByZScore<-function(SimilarityObj){
  ZData<-ZScoreContributionFunction(SimilarityObj$NullRankInfoFinal,SimilarityObj$RankInfoFinal)
  ZMeans<-rowMeans(ZData,na.rm=T)
  ZMeans[is.na(ZMeans)]<-0
  return(sort(ZMeans))
}
