#' ModalityZ for determining z scores in modality data
#'
#' @param ModeObject modality data from Modality()
#' @param PropObject proportionality data from a MembershipProportion function
#' @param npercent determines allowable variance. (Default is 10%) Default=0.1
#'
#' @return printed statements about use
#' @export
#'
#' @examples #
ModalityZ<-function(ModeObject,PropObject,npercent=.1){
  samples<-dim(PropObject)[2]
  PIdx<-which(ModeObject[[1]]$N_Members>=npercent*sum(ModeObject[[1]]$N_Members))
  reject<-list()
  if (dim(PropObject)[1]>1){
    if (length(PIdx)>1){
      for (j in 1:samples){
        xMN<-apply(PropObject[PIdx,-j]/rowSums(PropObject[PIdx,-j]),1,mean)
        xSD<-apply(PropObject[PIdx,-j]/rowSums(PropObject[PIdx,-j]),1,sd)
        xZ<-(PropObject[PIdx,j]/rowSums(PropObject[PIdx,-j])-xMN)/xSD
        if (any(abs(xZ)>=3)){
          zidx<-which(abs(xZ)>=3)
          ztemp<-matrix(nrow=2,ncol=length(zidx))
          ztemp[2,]<-xZ[zidx]
          ztemp[1,]<-ModeObject[[1]][PIdx,"Peak_x"][zidx]
          reject[[j]]<-ztemp
        }
      }
    } else {
      for (j in 1:samples){
        xMN<-mean(PropObject[PIdx,-j]/sum(PropObject[PIdx,-j]))
        xSD<-sd(PropObject[PIdx,-j]/sum(PropObject[PIdx,-j]))
        xZ<-(PropObject[PIdx,j]/sum(PropObject[PIdx,-j])-xMN)/xSD
        if (any(abs(xZ)>=3)){
          zidx<-which(abs(xZ)>=3)
          ztemp<-matrix(nrow=2,ncol=length(zidx))
          ztemp[2,]<-xZ[zidx]
          ztemp[1,]<-ModeObject[[1]][PIdx,'Peak_x'][zidx]
          reject[[j]]<-ztemp
        }
      }
    }
  }
  if (length(reject)>0){
    return(reject)
  } else {
    return(NULL)
  }
}
