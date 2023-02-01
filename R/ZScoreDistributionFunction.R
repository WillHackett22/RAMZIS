#' ZScoreDistributionFunction determines the Z-score of all values in the second distribution compared to the mean and SD of the first
#'
#' @param ReferenceDis First distribution for reference
#' @param TestDis Second distribution for testing
#'
#' @return distribution of z scores
#' @export
#'
#' @examples #
ZScoreDistributionFunction<-function(ReferenceDis,TestDis){
  ReferenceMean<-mean(ReferenceDis,na.rm=T)
  ReferenceSD<-sd(ReferenceDis,na.rm=T)
  Zdis<-sapply(TestDis,ZScoreFunction,ReferenceMean=ReferenceMean,ReferenceSD=ReferenceSD)
  return(Zdis)
}
