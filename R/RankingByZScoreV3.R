#' RankingByZScoreV3 determines rank
#'
#' @param TestSimObj Test Similarity object
#' @param NullSimObj Null Similarity object
#'
#' @return rankings
#' @export
#'
#' @examples #
RankingByZScoreV3<-function(TestSimObj,NullSimObj,RankingInfo="Numerator"){
  ZData<-ZScoreContributionFunction(NullSimObj[[RankingInfo]],TestSimObj[[RankingInfo]])
  ZMeans<-rowMeans(ZData,na.rm=T)
  return(sort(ZMeans))
}
