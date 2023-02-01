#' ManhattanVectorDistance determines the manhattan distance of a vector compared to a dataframe
#'
#' @param vector1 is a vector of numbers
#' @param dataframe1 is a matrix/dataframe of numbers
#'
#' @return the distance
#' @export
#'
#' @examples #
ManhattanVectorDistance<-function(vector1,dataframe1){
  Out<-abs(vector1-dataframe1)
  return(Out)
}
