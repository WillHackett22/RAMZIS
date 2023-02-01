#' SimPlotGeneral produces a plot of the densities of the similarity of the Test and Null Similarities
#'
#' @param PlotTitle PlaceHolder
#' @param dfSim PlaceHolder
#' @param OverlapData PlaceHolder
#' @param legend.bool PlaceHolder
#'
#' @return PlaceHolder
#' @export
#'
#' @examples #
SimPlotGeneral<-function(PlotTitle,dfSim,OverlapData,legend.bool=T){
  RefName<-'Null'
  p<-SimPlotGGBase(PlotTitle,dfSim,RefName,OverlapData,Confidence=NULL,legend.bool)
  return(p)
}
