#' Modality for determining multi-modality
#'
#' @param SimDist The similarity distribution you would like to find the modality of
#' @param threshold The ratio needed to identify a different peak. Default =0.01
#'
#' @return Returns a list of data on the peaks in the similarity distribution
#' @export
#'
#' @examples #
Modality<-function(SimDist,threshold=.01){
  SimDens<-density(SimDist)
  #build sim density distribution
  LocalMinimaIdx<-which(diff(diff(SimDens$y)<=0)<0)+1
  LocalMaximaIdx<-which(diff(diff(SimDens$y)>=0)<0)+1
  #find probable minima and maxima
  MinMat<-data.frame(matrix(NA,nrow=length(LocalMaximaIdx),ncol=9))
  colnames(MinMat)<-c('Peak_x','dPeakLeft','dPeakRight','LeftCheck','RightCheck','Overall','N_Members','InterPeakDistance','Peak_y')
  MinMat[,1]<-SimDens$x[LocalMaximaIdx]
  MinMat[,9]<-SimDens$y[LocalMaximaIdx]
  #Peak Locations
  MinMat[,2]<-(SimDens$y[LocalMaximaIdx]-c(0,SimDens$y[c(LocalMinimaIdx)]))
  #Height Differentials before peak
  MinMat[,3]<-(SimDens$y[LocalMaximaIdx]-c(SimDens$y[c(LocalMinimaIdx)],0))
  #height differentials after peak
  MinMat[,4]<-abs(MinMat[,2])<=threshold*SimDens$y[LocalMaximaIdx]
  MinMat[,5]<-abs(MinMat[,3])<=threshold*SimDens$y[LocalMaximaIdx]
  #chance that the height differentials are by random chance

  TempMin<-c(0,SimDens$x[LocalMinimaIdx],1)
  for (j in 1:length(LocalMaximaIdx)){
    temppeak<-which(TempMin[j]<=SimDist & SimDist<=TempMin[j+1])
    MinMat[j,7]<-length(temppeak)
    #XSD<-sd(SimDist[temppeak])
    #MinMat[j,8]<-abs(TempMin[j+1]-TempMin[j])<=3*XSD
    MinMat[j,8]<-FALSE
  }
  MinMat[,6]<-(sum(MinMat[,4])>0|sum(MinMat[,5])>0)

  while (sum(MinMat[,6])>0){
    rowpick<-which(MinMat[,4]|MinMat[,5])
    targeti<-which(rowSums(abs(MinMat[rowpick,c(2,3)])==min(abs(MinMat[rowpick,c(2,3)])))>0)[1]
    targeti<-rowpick[targeti]
    #find which peaks have the smallest differential, ties go to the first
    targetl<-which(abs(MinMat[targeti,c(2,3)])==min(abs(MinMat[targeti,c(2,3)])))[1]-2
    #find whether it is the preceding or secondary trough, ties go to the preceding
    #if first peak and trough, remove peak and next trough
    #if last peak and last trough, remove peak and prior trough
    #else remove peak and lowest trough
    if (targetl==-1 & targeti==1){
      LocalMaximaIdx<-LocalMaximaIdx[-targeti]
      LocalMinimaIdx<-LocalMinimaIdx[-1]
      #remove peak
    } else if (targetl==0 & targeti==length(LocalMaximaIdx)){
      LocalMaximaIdx<-LocalMaximaIdx[-targeti]
      LocalMinimaIdx<-LocalMinimaIdx[-length(LocalMinimaIdx)]
      #remove peak
    } else {
      LocalMaximaIdx<-LocalMaximaIdx[-targeti]
      #remove peak
      tempidx<-targeti+targetl
      LocalMinimaIdx<-LocalMinimaIdx[-tempidx]
      #remove lower trough
    }
    if (length(LocalMinimaIdx)>0){
      MinMat<-data.frame(matrix(NA,nrow=length(LocalMaximaIdx),ncol=9))
      colnames(MinMat)<-c('Peak_x','dPeakLeft','dPeakRight','LeftCheck','RightCheck','Overall','N_Members','InterPeakDistance','Peak_y')
      MinMat[,1]<-SimDens$x[LocalMaximaIdx]
      MinMat[,9]<-SimDens$y[LocalMaximaIdx]
      #Peak Locations
      MinMat[,2]<-(SimDens$y[LocalMaximaIdx]-c(0,SimDens$y[c(LocalMinimaIdx)]))
      #Height Differentials before peak
      MinMat[,3]<-(SimDens$y[LocalMaximaIdx]-c(SimDens$y[c(LocalMinimaIdx)],0))
      #height differentials after peak
      MinMat[,4]<-abs(MinMat[,2])<=threshold*SimDens$y[LocalMaximaIdx]
      MinMat[,5]<-abs(MinMat[,3])<=threshold*SimDens$y[LocalMaximaIdx]
      #chance that the height differentials are by random chance

      TempMin<-c(0,SimDens$x[LocalMinimaIdx],1)
      for (j in 1:length(LocalMaximaIdx)){
        temppeak<-which(TempMin[j]<=SimDist & SimDist<=TempMin[j+1])
        MinMat[j,7]<-length(temppeak)
        #XSD<-sd(SimDist[temppeak])
        #MinMat[j,8]<-abs(TempMin[j+1]-TempMin[j])<=3*XSD
        MinMat[j,8]<-FALSE
      }
      MinMat[,6]<-(sum(MinMat[,4])>0|sum(MinMat[,5])>0)
    } else {
      MinMat<-data.frame(matrix(NA,nrow=length(LocalMaximaIdx),ncol=9))
      colnames(MinMat)<-c('Peak_x','dPeakLeft','dPeakRight','LeftCheck','RightCheck','Overall','N_Members','InterPeakDistance','Peak_y')
      MinMat[,1]<-SimDens$x[LocalMaximaIdx]
      MinMat[,9]<-SimDens$y[LocalMaximaIdx]
      #Peak Locations
      MinMat[,2]<-(SimDens$y[LocalMaximaIdx]-c(0,SimDens$y[c(LocalMinimaIdx)]))
      #Height Differentials before peak
      MinMat[,3]<-(SimDens$y[LocalMaximaIdx]-c(SimDens$y[c(LocalMinimaIdx)],0))
      #height differentials after peak
      MinMat[,4]<-FALSE
      MinMat[,5]<-FALSE
      #chance that the height differentials are by random chance

      TempMin<-c(0,SimDens$x[LocalMinimaIdx],1)
      for (j in 1:length(LocalMaximaIdx)){
        temppeak<-which(TempMin[j]<=SimDist & SimDist<=TempMin[j+1])
        MinMat[j,7]<-length(temppeak)
        #XSD<-sd(SimDist[temppeak])
        #MinMat[j,8]<-abs(TempMin[j+1]-TempMin[j])<=3*XSD
        MinMat[j,8]<-FALSE
      }
      MinMat[,6]<-(sum(MinMat[,4])>0|sum(MinMat[,5])>0)
    }


  }
  return(list(MinMat,'Maxima'=LocalMaximaIdx,'Minima'=LocalMinimaIdx))
}
