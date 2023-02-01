#' MatrixMerge_Helper
#'
#' @param m1 Matrix 1
#' @param m2 Matrix 2
#'
#' @return A matrix combined by row names
#'
#' @examples
MatrixMerge_Helper<-function(m1,m2){
  out<-merge(m1,m2,by=0,all=TRUE)
  row.names(out)<-out[,1]
  out<-out[,-1]
  out[is.na(out)]<-0
  return(out)
}
