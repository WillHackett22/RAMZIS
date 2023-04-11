#' SimPlotInternal produces a similarity plot of one of the two internal distributions
#'
#' @param PlotTitle The Plot Title
#' @param SimilarityObject Similarity Object Produced by RAMZISMain()
#' @param Dataset Selects which dataset to produce the internal plot of. Dataset=1 produces internal of the first dataset given. Dataset=2 produces the internal of the second dataset.
#' @param legend.bool Boolean to produce legend or not
#' @param themepick Default=theme_RAMZIS. Theme object from ggplot2 for formatting.
#' @param zopt Indicates use of Z-score for simulation validity over percentile. Default=T
#'
#' @return Internal Similarity Plot
#' @export
#'
#' @examples #
SimPlotInternal<-function(PlotTitle,SimilarityObject,Dataset,legend.bool=T,themepick=theme_RAMZIS,zopt=T,...){
  if (Dataset==1){
    Refname<-'Internal1'
  } else if (Dataset==2){
    Refname<-'Internal2'
  } else {
    print('Error: Dataset must be 1 or 2')
    quit()
  }
  dfSim<-SimilarityObject$PlotData[[Refname]]
  OverlapData<-OverlapCalculator(SimilarityObject$Similarity[[Refname]],SimilarityObject$Similarity$Test)
  ConfidenceData<-SimilarityObject$Summary$QualityChecks[[Refname]]$Confidence
  ObservedData<-ObservedSimilarityStats(SimilarityObject$Similarity$Test,SimilarityObject$Similarity$Actual)
  p<-SimPlotGGBase(PlotTitle,dfSim,Refname,OverlapData,ObservedData,Confidence=ConfidenceData,legend.bool=T,...)
  return(p)
}
