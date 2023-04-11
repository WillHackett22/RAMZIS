#' RankingByZScore determines rank
#'
#' @param SimilarityObj Similarity object from RAMZISMain()
#' @param RankInfo Which information to use for ranking. Default="RankInfo", Alternative use is "WeightedContributions"
#'
#' @return rankings
#' @export
#'
#' @examples #
RankingByZScore<-function(SimilarityObj,RankInfo="RankInfo"){
  if (RankInfo=="RankInfo"){
    RankInfoData<-SimilarityObj[[RankInfo]]$RankData
  } else {
    RankInfoData<-SimilarityObj[[RankInfo]]
  }

  ZData<-ZScoreContributionFunction(RankInfoData$Null,RankInfoData$Test)
  ZMeans<-rowMeans(ZData,na.rm=T)
  return(sort(ZMeans))
}
