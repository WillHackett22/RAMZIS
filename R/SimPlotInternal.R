#' SimPlotInternal produces a similarity plot of one of the two internal distributions
#'
#' @param PlotTitle The Plot Title
#' @param dfSim Similarity Object Produced by RAMZISMain()
#' @param Dataset Selects which dataset to produce the internal plot of. Dataset=1 produces internal of the first dataset given. Dataset=2 produces the internal of the second dataset.
#' @param OverlapData OverlapData from RAMZISMain() OverlapCalculator()
#' @param Confidence From InternalConfidenceScore()
#' @param legend.bool Boolean to produce legend or not
#'
#' @return Internal Similarity Plot
#' @export
#'
#' @examples #
SimPlotInternal<-function(PlotTitle,dfSim,Dataset,OverlapData,Confidence=NULL,legend.bool=T){
  if (Dataset==1){
    Refname<-'Internal1'
  } else if (Dataset==2){
    Refname<-'Internal2'
  } else {
    print('Error: Dataset must be 1 or 2')
    quit()
  }
  p<-SimPlotGGBase(PlotTitle,dfSim,RefName,OverlapData,Confidence=NULL,legend.bool=T)
  return(p)
}
