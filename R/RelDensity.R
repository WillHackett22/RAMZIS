#' Converts a similarity to relative density
#'
#' @param SimDist Tanimoto Similarity Vector
#'
#' @return Relative Density object
#' @export
#'
#' @examples #
RelDensity<-function(SimDist){
  SimDens<-DENSITY_RD(SimDist)
  SDensy<-100*SimDens$y/sum(SimDens$y)
  return(list(y=SDensy,x=SimDens$x))
}
