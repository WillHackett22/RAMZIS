#' MatrixMerge_Helper
#'
#' @param m1 Matrix 1
#' @param m2 Matrix 2
#' @param na.rm Default=T turns NA into 0
#'
#' @return A matrix combined by row names
#'
#' @examples #
MatrixMerge_Helper<-function(m1,m2,na.rm=T){
  out<-merge(m1,m2,by=0,all=TRUE)
  row.names(out)<-out[,1]
  out<-out[,-1]
  if (na.rm){
    out[is.na(out)]<-0
  }
  return(out)
}
