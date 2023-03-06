#' MVCount counts the number of missing values
#'
#' @param vec Vector to count zeros in
#' @param val Value for something to be considered missing equal to. Default=0
#'
#' @return Number of zeros in a vector
#'
#' @examples #
MVCount<-function(vec,val=0){
  return(sum(vec==val))
}
