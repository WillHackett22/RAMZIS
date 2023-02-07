#' ObservedSimilarityStats gets the pertinent summary statistics about the Observed Similarity
#'
#' @param TestSimilarity All Similarity Values in the Test Distribution
#' @param ActualSim The Observed Similarity
#'
#' @return A list of the Observed similarity, its percentile location, and its zscore
#' @export
#' @import stats
#'
#' @examples #
ObservedSimilarityStats<-function(TestSimilarity,ActualSim){
  Perc<-round(ecdf(TestSimilarity)(ActualSim),2)*100
  ActualZ<-ZScoreFunction(mean(TestSimilarity,na.rm=T),stats::sd(TestSimilarity,na.rm=T),ActualSim)
  out<-list("Actual"=ActualSim,"Percentile"=Perc,"ZScore"=ActualZ)
  return(out)
}
