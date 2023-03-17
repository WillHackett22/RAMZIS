#' GATHERER_HelperTIDY
#' Get all the information for similarity comparisons
#'
#' @param gl Glycopeptide List Deprecated
#' @param combo Sampling
#' @param df Dataset
#' @param MVCorrection If False, it will ignore missing values rather than count towards the overall average. Default=TRUE
#'
#' @return Gets data for use in matrix similarity comparison
#'
#' @examples #
GATHERER_HelperTIDY<-function(gl,combo,df,MVCorrection){
  nc<-nrow(combo)
  gl<-row.names(df)
  cl<-row.names(combo)
  Partial<-data.frame(matrix(0,nrow=length(gl)*nc,ncol=5))
  colnames(Partial)<-c("Identification","SampleID","Value","Presence","Type")
  Partial$Identification<-rep(gl,nc)
  Partial$Type<-"Partial"
  Partial$SampleID<-rep(cl,each=length(gl))
  Self<-data.frame(matrix(0,nrow=nc,ncol=3))
  colnames(Self)<-c("SampleID","Value","Type")
  Self$Type<-"Self"
  Self$SampleID<-cl
  for (j in 1:nc){
    #separate datasets for combinations
    temp<-data.frame(df[,unlist(combo[j,])]) # subset of first data
    if (MVCorrection!=TRUE){
      temp[temp==0]<-NA
    }
    pdx<-which(Partial$SampleID==cl[j])
    Partial$Value[pdx]<-rowMeans(temp,na.rm =TRUE)
    Partial$Presence[pdx]<-1-apply(temp,1,MVCount)/ncol(temp)
    sdx<-which(Self$SampleID==cl[j])
    Self$Value[sdx]<-sum((rowMeans(temp,na.rm =TRUE))^2)
  }
  if (any(is.na(Partial$Presence))){
    Partial$Presence[which(is.na(Partial$Presence))]<-0
  }
  out<-list("Self"=Self,"Partial"=Partial)
  return(out)
}
