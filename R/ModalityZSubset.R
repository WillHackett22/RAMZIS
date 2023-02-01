#' ModalityZSubset is deprecated and doesn't matter
#'
#' @param ModeObject modality object
#' @param PropObject proportionality object
#' @param npercent threshold
#'
#' @return ignore
#'
#' @examples #
ModalityZSubset<-function(ModeObject,PropObject,npercent=.1){
  samples<-dim(PropObject)[2]
  PIdx<-which(ModeObject[[1]]$N_Members>=npercent*sum(ModeObject[[1]]$N_Members))
  reject<-list()
  i<-1
  if ((samples-2)>=2){
    mMax<-samples-2
    mseq<-mMax:2
  }

  if (dim(PropObject)[1]>1){
    if (length(PIdx)>1){
      for (l in 1:length(mseq)){
        subsets<-combn(samples,mseq[l])
        slen<-dim(subsets)[2]
        for (j in 1:slen){
          sub<-subsets[,j]
          xMN<-apply(PropObject[PIdx,-sub]/rowSums(PropObject[PIdx,-sub]),1,mean)
          xSD<-apply(PropObject[PIdx,-sub]/rowSums(PropObject[PIdx,-sub]),1,sd)
          xZ<-(PropObject[PIdx,sub]/rowSums(PropObject[PIdx,-sub])-xMN)/xSD
          if (any(abs(xZ)>=3)){
            if (any(rowSums(abs(xZ)>=3)==length(sub))){
              zidx<-which(rowSums(abs(xZ)>=3)==length(sub))
              ztemp<-matrix(nrow=2,ncol=length(zidx))
              ztemp[2,]<-xZ[zidx]
              ztemp[1,]<-ModeObject[[1]][PIdx,"Peak_x"][zidx]
              reject[[i]]<-list()
              reject[[i]][[1]]<-as.data.frame(ztemp)
              reject[[i]][[2]]<-sub
              i<-i+1
            }

          }
        }
      }

    } else {
      for (l in 1:length(mseq)){
        subsets<-combn(samples,mseq[l])
        slen<-dim(subsets)[2]
        for (j in 1:slen){
          sub<-subsets[,j]
          xMN<-mean(PropObject[PIdx,-sub]/sum(PropObject[PIdx,-sub]),na.rm=T)
          xSD<-sd(PropObject[PIdx,-sub]/sum(PropObject[PIdx,-sub]),na.rm=T)
          xZ<-(PropObject[PIdx,sub]/sum(PropObject[PIdx,-sub])-xMN)/xSD
          if (any(abs(xZ)>=3)){
            if (any(rowSums(abs(xZ)>=3)==length(sub))){
              zidx<-which(rowSums(abs(xZ)>=3)==length(sub))
              ztemp<-matrix(nrow=2,ncol=length(zidx))
              ztemp[2,]<-xZ[zidx]
              ztemp[1,]<-ModeObject[[1]][PIdx,"Peak_x"][zidx]
              reject[[i]]<-list()
              reject[[i]][[1]]<-as.data.frame(ztemp)
              reject[[i]][[2]]<-sub
              i<-i+1
            }

          }
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
