#' SimPlotGeneral produces a plot of the densities of the similarity of the Test and Null Similarities
#'
#' @param PlotTitle
#' @param dfSim
#' @param OverlapData
#' @param legend.bool
#'
#' @return
#' @export
#'
#' @examples
SimPlotGeneral<-function(PlotTitle,dfSim,OverlapData,legend.bool=T){
  RefName<-'Null'
  p<-SimPlotGGBase(PlotTitle,dfSim,RefName,OverlapData,Confidence=NULL,legend.bool)
  return(p)
}
