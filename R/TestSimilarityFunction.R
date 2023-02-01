#' TestSimilarityFunction
#'
#' @param df1 Dataset 1
#' @param combo1 Sampling of dataset 1
#' @param gl1 Glycopeptide list ofdataset 1
#' @param df2 Dataset 2
#' @param combo2 Sampling of dataset 2
#' @param gl2 glycopeptide list of dataset 2
#' @param MVCorrection If False, it will ignore missing values rather than count towards the overall average. Default=TRUE
#' @param mn The scaling factor. Default will use 1+mean(presence) of a glycopeptide/identification. Setting this to a number will override that process.
#'
#' @return
#' @export
#'
#' @examples
TestSimilarityFunction<-function(df1,combo1,gl1,df2,combo2,gl2,MVCorrection,mn=FALSE){
  #iterate through all combinations used in within of df 1
  T1_Obj<-GATHERER_Helper(gl1,combo1,df1,MVCorrection)
  T1_<-T1_Obj$Partial
  T10<-T1_Obj$Self
  PHold1<-T1_Obj$Presence
  #iterate through all combinations used in within of df 2
  T_1Obj<-GATHERER_Helper(gl2,combo2,df2,MVCorrection)
  T_1<-T_1Obj$Partial
  T01<-T_1Obj$Self
  PHold1<-T_1Obj$Presence
  #Create A*B Expression
  T__Hold<-MatrixMerge_Helper(T1_Obj$Partial,T_1Obj$Partial)
  PHold<-MatrixMerge_Helper(T1_Obj$Presence,T_1Obj$Presence)
  TRef<-row.names(T__Hold)

  #Bring T11 related terms together
  TempSimObject<-SimilarityCalculation(T__Hold,T10,T01,PHold,nrow(combo1),nrow(combo2),mn)
  return(TempSimObject)
}
