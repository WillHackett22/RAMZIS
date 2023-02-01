#' PresenceCalc_Helper
#'
#' @param gl glycopeptide list
#' @param df dataset
#'
#' @return The presence of each glycopeptide in a dataset
#' @export
#'
#' @examples
PresenceCalc_Helper<-function(gl,df){
  out<-data.frame(apply(df,1,PresenceCalcBase))
  return(out)
}

#' PresenceCalcBase finds the presence of a single vector
#'
#' @param vec
#'
#' @return
#'
#' @examples
PresenceCalcBase<-function(vec){
  if (any(is.na(vec))){
    vec[is.na(vec)]<-0
  }
  pres<-1-sum(vec==0)/length(vec)
  return(pres)
}
