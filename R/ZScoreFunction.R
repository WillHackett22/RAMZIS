#' ZScoreFunction determines the z score.
#'
#' @param ReferenceMean A mean
#' @param ReferenceSD the Standard Deviation
#' @param TestNum The Number you want to know the z-score of
#'
#' @return the z-score
#' @export
#'
#' @examples #
ZScoreFunction<-function(ReferenceMean,ReferenceSD,TestNum){
  Zscore<-(TestNum-ReferenceMean)/ReferenceSD
  return(Zscore)
}
