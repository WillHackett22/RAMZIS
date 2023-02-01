#' VectorMatrixMultiplication determines the multiplication of a dataframe against a single vector
#'
#' @param vector1 numeric vector
#' @param dataframe1 numeric dataframe/matrix
#'
#' @return matrix of multiplications
#' @export
#'
#' @examples #
VectorMatrixMultiplication<-function(vector1,dataframe1){
  Out<-(vector1*dataframe1)
  return(Out)
}
