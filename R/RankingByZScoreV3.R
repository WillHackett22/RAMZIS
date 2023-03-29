#' RankingByZScoreV3 determines rank
#'
#' @param TestSimObj Test Similarity object
#' @param NullSimObj Null Similarity object
#'
#' @return rankings
#' @export
#'
#' @examples #
RankingByZScoreV3<-function(TestSimObj,NullSimObj){
  ZData<-ZScoreContributionFunction(NullSimObj$Numerator,TestSimObj$Numerator)
  ZMeans<-rowMeans(ZData,na.rm=T)
  return(sort(ZMeans))
}
