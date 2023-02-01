#' NullSimilarityFunction
#'
#' @param df1 Original Dataset 1
#' @param df2 Original Dataset 2
#' @param ncombos Null Dataset Samplings
#' @param MVCorrection If False, it will ignore missing values rather than count towards the overall average. Default=TRUE
#' @param mn  The scaling factor. Default will use 1+mean(presence) of a glycopeptide/identification. Setting this to a number will override that process.
#'
#' @return Distribution of Null similarities and their contributions
#' @export
#'
#' @examples
NullSimilarityFunction<-function(df1,df2,ncombos,MVCorrection,mn){
  gl1<-c(row.names(df1))
  gl2<-c(row.names(df2))
  glycojoint<-unique(c(gl1,gl2))
  ncn1dim<-dim(ncombos$NDis1)
  ncn2dim<-dim(ncombos$NDis2)
  ncrr<-ncn1dim[1]*ncn2dim[1]

  #first set of null distributions
  T1_NObj<-GATHERER_Helper(glycojoint,ncombos$NDis1,mergedf,MVCorrection)
  T1_N<-T1_NObj$Partial
  PHoldN1<-T1_NObj$Presence

  #second set of null distributions
  T_1NObj<-GATHERER_Helper(glycojoint,ncombos$NDis1,mergedf,MVCorrection)
  T_1N<-T_1NObj$Partial
  PHoldN2<-T_1NObj$Presence

  #Bring T11N related terms together
  T__HoldN<-MatrixMerge_Helper(T1_NObj$Partial,T_1NObj$Partial)
  PHoldN<-MatrixMerge_Helper(T1_NObj$Presence,T_1NObj$Presence)

  ncrtot<-ncn1dim[1]+ncn2dim[1]
  TempSimObject<-SimilarityCalculation(T__HoldN,T1_NObj$Self,T_1NObj$Self,PHoldN,ncn1dim[1],ncn2dim[1],mn)
}
