#' TestSimilarityFunctionTIDY
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
#' @return Similarity between two matrices
#' @export
#'
#' @examples #
TestSimilarityFunctionTIDY<-function(df1,combo1,gl1,df2,combo2,gl2,MVCorrection,mn=FALSE){
  #iterate through all combinations used in within of df 1
  T1_Obj<-GATHERER_HelperTIDY(gl1,combo1,df1,MVCorrection)
  #iterate through all combinations used in within of df 2
  T_1Obj<-GATHERER_HelperTIDY(gl2,combo2,df2,MVCorrection)
  gl<-unique(c(gl1,gl2))
  #Bring T11 related terms together
  TempSimObject<-SimilarityCalculationTIDY(T1_Obj,T_1Obj,nrow(combo1),nrow(combo2),gl,mn)
  return(TempSimObject)
}
