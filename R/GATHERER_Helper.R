#' GATHERER_Helper
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
GATHERER_Helper<-function(gl,combo,df,MVCorrection){
  out<-c()
  nc<-nrow(combo)
  gl<-row.names(df)
  out$Partial<-data.frame(matrix(0,nrow=length(gl),ncol=nc))
  row.names(out$Partial)<-gl
  out$Self<-data.frame(matrix(0,nrow=1,ncol=nc))
  out$Presence<-data.frame(matrix(0,nrow=length(gl),ncol=nc))
  row.names(out$Presence)<-gl
  for (j in 1:nc){
    #separate datasets for combinations
    temp<-data.frame(df[,unlist(combo[j,])]) # subset of first data
    if (MVCorrection!=TRUE){
      temp[temp==0]<-NA
    }
    out$Partial[,j]<-rowMeans(temp,na.rm =TRUE)
    out$Self[,j]<-sum((rowMeans(temp,na.rm =TRUE))^2)
    for (m in 1:length(gl)){
      nacheck<-1-sum(temp[gl[m],]==0)/ncol(temp)
      if (is.na(nacheck)){
        nacheck<-0
      }
      out$Presence[gl[m],j]<-nacheck
    }
  }
  return(out)
}
