#' SimilarityCalculation_Singular
#' For calculating two datasets similarity without sampling. It is used in RAMZIS_Main to calculate the observed similarity
#'
#' @param df1 The first dataframe to be compared
#' @param df2 The second dataframe to be compared
#' @param mn The scaling factor. Default will use 1+mean(presence) of a glycopeptide/identification. Setting this to a number will override that process.
#'
#' @return
#' @export
#'
#' @examples
SimilarityCalculation_Singular<-function(df1,df2,mn=FALSE){
  gl1<-row.names(df1)
  gl2<-row.names(df2)
  T1_a<-rowMeans(df1,na.rm =TRUE)
  T10a<-sum((rowMeans(df1,na.rm =TRUE))^2)
  T_1a<-rowMeans(df2,na.rm =TRUE)
  T01a<-sum((rowMeans(df2,na.rm =TRUE))^2)
  presence1<-PresenceCalc_Helper(gl1,df1)
  presence2<-PresenceCalc_Helper(gl2,df2)
  presenceA<-MatrixMerge_Helper(presence1,presence2)
  THA<-MatrixMerge_Helper(T1_a,T_1a)
  TRef<-row.names(THA)
  lT<-length(TRef)
  if (mn==FALSE){
    KTerm<-data.frame(matrix((1+rowMeans(presenceA,na.rm=T)),nrow=lT,ncol=1))
    row.names(KTerm)<-row.names(presenceA)
  } else {
    KTerm<-rep(mn,lT)
    row.names(KTerm)<-row.names(presenceA)
  }
  dTA<-data.frame(matrix(abs(THA[,1]-THA[,2]),nrow=lT,ncol=1))
  row.names(dTA)<-TRef
  tanmatA<-data.frame(matrix(0,ncol=lT,nrow=1))
  tanmatAFinal<-data.frame(matrix(0,ncol=lT,nrow=1))
  colnames(tanmatA)<-TRef
  colnames(tanmatAFinal)<-TRef
  T11A<-data.frame(matrix(0,nrow=lT,ncol=1))
  row.names(T11A)<-TRef
  for (m in 1:lT){
    T11A[TRef[m],1]<-THA[TRef[m],1]*THA[TRef[m],2]*(KTerm[TRef[m],1]^(-dTA[TRef[m],1]))
    tanmatA[1,TRef[m]]<-unlist(T11A[TRef[m],])
  }
  for (m in 1:length(TRef)){
    tanmatAFinal[1,TRef[m]]<-unlist(T11A[TRef[m],]/(T10a+T01a-sum(T11A,na.rm=T)))
  }
  taniActualF<-sum(tanmatAFinal)
  out<-list("Similarity"=taniActualF,"Contribution"=tanmatAFinal)
  return(out)
}
