#' InternalConfidenceScore
#' Finds the Confidence Score of a Test Similarity Distribution relative to one of its Internal Similarity Distribution
#'
#' @param ReferenceDis One of the Internal Similarity Distributions
#' @param TestDis The Test Similarity Distribution
#' @param OverlapData The known overlap of these distributions from OverlapCalculator()
#'
#' @return PlaceHolder
#' @export
#'
#' @examples #
InternalConfidenceScore<-function(ReferenceDis,TestDis,OverlapData=NULL){
  if (is.null(OverlapData)){
    OverlapData=OverlapCalculator(ReferenceDis,TestDis)
  }
  Alpha<-OverlapData$FP
  Beta<-OverlapData$FN
  deltaR<-abs(mean(ReferenceDis,na.rm=T)-mean(TestDis,na.rm = T))
  Confidence<-round(deltaR/sd(ReferenceDis,na.rm=T)*10^(-Alpha-Beta),2)
  return(Confidence)
}
