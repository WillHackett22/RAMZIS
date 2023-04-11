#' DENSITY_RD
#' This is the default density function of RAMZIS and is really not meant for user usage, but is available for validation.
#'
#' @param data A vector of data to be turned into a density function
#' @param lb The lower bound of the density function. Default=-0.1
#' @param ub The upper bound of the density function. Default=1.1
#' @param k The scaling factor by which densitys are adjusted. The deafult converts from proportion to percentage. Default=100
#'
#' @return Density Object akin to base::density()
#' @export
#' @import stats
#'
#' @examples
#' dist<-rbeta(1000,2,2)
#' #dens<-DENSITY_RD(dist)
DENSITY_RD<-function(data,lb=-0.1,ub=1.1,k=100){
  dens_out<-stats::density(data,from=lb,to=ub,na.rm=T)
  if (max(data)==0){
    dens_out$y<-rep(0,length(dens_out$y))
    dens_out$y[which(dens_out$x>=-0.001 & dens_out$x<=0.001)]<-1
  } else{
    dens_out$y<-k*dens_out$y/sum(dens_out$y)
  }
  return(dens_out)
}
