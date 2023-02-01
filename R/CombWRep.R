#' CombWRep determines the number of combinations possible with replacement
#'
#' @param n number of possible members
#' @param r number of sampled points
#'
#' @return number possible
#' @export
#'
#' @examples
#' #if there are 8 samples and we want to sample 4 with replacement, how many possible combinations are there
#' CombWRep(8,4)
CombWRep<-function(n,r){
  return(factorial(n+r-1)/(factorial(n-1)*factorial(r)))
}
