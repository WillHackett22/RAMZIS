#' VectorMatrixMean determines the mean of a dataframe against a single vector
#'
#' @param vector1 numeric vector
#' @param dataframe1 numeric dataframe/matrix
#'
#' @return matrix of means
#' @export
#'
#' @examples #
VectorMatrixMean<-function(vector1,dataframe1){
  Out<-(vector1+dataframe1)/2
  return(Out)
}
