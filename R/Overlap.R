#' Overlap: DEPRICATED: finds the overlap statistics of two density distributions
#'
#' @param Dis1 Distribution 1
#' @param Dis2 Distribution 2
#'
#' @return Overlap statistics: Total Percent, False Positive (Alpha), and False Negative (Beta)
#' @export
#'
#' @examples #
Overlap<-function(Dis1,Dis2){
  density1<-density(Dis1,from=-0.1,to=1.1,na.rm=T)
  density2<-density(Dis2,from=-0.1,to=1.1,na.rm=T)
  df <- merge(
    as.data.frame(density1[c("x", "y")]),
    as.data.frame(density2[c("x", "y")]),
    by = "x", suffixes = c(".A", ".B")
  )
  df$comp <- as.numeric(df$y.A > df$y.B)
  df$cross <- c(NA, diff(df$comp))
  tempidx<-which(df$cross!=0)
  CPoint<-tempidx[which(max(df$y.A[tempidx])==df$y.A[tempidx])]
  XPoint<-df$x[CPoint]
  if (length(CPoint)==0){
    Overlap<-0
    Alpha<-0
    Beta<-0
  } else if (df$y.A[CPoint]<(10^-10)){
    Overlap<-0
    Alpha<-0
    Beta<-0
  } else {
    Area1<-(1-ecdf(Dis1)(XPoint))
    Area2<-ecdf(Dis2)(XPoint)
    Overlap1<-round(Area1+Area2,2)
    Overlap<-round(100*Overlap1/(2-Overlap1),2)
    Alpha<-round(Area1/(2-Overlap1),2)
    Beta<-round(Area2/(2-Overlap1),2)
  }
  return(list('PercOverlap'=Overlap,'Alpha'=Alpha,'Beta'=Beta))
}
