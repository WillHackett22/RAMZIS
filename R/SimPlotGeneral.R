#' SimPlotGeneral produces a plot of the densities of the similarity of the Test and Null Similarities with the SimilarityObject produced by RAMZISMain()
#'
#' @param PlotTitle The Plot Title
#' @param SimilarityObject Similarity Object Produced by RAMZISMain()
#' @param legend.bool Boolean to produce legend or not
#' @param themepick Default=theme_RAMZIS. Theme object from ggplot2 for formatting.
#' @param zopt Indicates use of Z-score for simulation validity over percentile. Default=T
#'
#' @return PlaceHolder
#' @export
#'
#' @examples #
SimPlotGeneral<-function(PlotTitle,SimilarityObject,legend.bool=T,themepick=theme_RAMZIS,zopt=T,...){
  Refname<-'Null'
  dfSim<-SimilarityObject$PlotData$General
  OverlapData<-OverlapCalculator(SimilarityObject$Similarity$Null,SimilarityObject$Similarity$Test)
  ObservedSimObj<-ObservedSimilarityStats(SimilarityObject$Similarity$Test,SimilarityObject$Similarity$Actual)
  p<-SimPlotGGBase(PlotTitle,dfSim,Refname,OverlapData,ObservedSimObj,Confidence=NULL,legend.bool,themepick=theme_RAMZIS,zopt=T,...)
  return(p)
}
