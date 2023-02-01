#' VectorMatrixAddition determines the multiplication of a dataframe against a single vector
#'
#' @param vector1 numeric vector
#' @param dataframe1 numeric dataframe/matrix
#'
#' @return matrix of additions
#' @export
#'
#' @examples #
VectorMatrixAddition<-function(vector1,dataframe1){
  Out<-(vector1+dataframe1)
  return(Out)
}
