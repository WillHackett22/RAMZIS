#' VectorMatrixAddition determines the multiplication of a dataframe against a single vector
#'
#' @param vector1 numeric vector
#' @param vector2 numeric dataframe/matrix
#'
#' @return matrix of additions
#' @export
#'
#' @examples #
VectorMatrixAddition<-function(vector1,vector2){
  Out<-(vector1+vector2)
  return(Out)
}
