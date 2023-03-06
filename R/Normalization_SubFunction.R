#' Normalization_SubFunction produces the normalization of a dataframe
#'
#' @param df Dataframe to be normalized
#' @param normvector Vector of NUMERIC values that is equal to sample size/number of df columns. Default='None'
#' @param logoption Boolean indicating use of log transformation. Default=TRUE
#'
#' @return Normalized dataframe
#' @export
#'
#' @examples #
Normalization_SubFunction<-function(df,normvector='None',logoption=TRUE){
  if (logoption){
    df<-log(df+1)
  }
  if (length(normvector)==ncol(df)){
    for (i in 1:ncol(df)){
      df[,i]<-as.numeric(df[,i])*normvector[i]
    }
  } else {
    if (normvector=='None'){
      dfout<-df
    } else if (normvector=='Relative'){
      for (i in 1:ncol(df)){
        df[,i]<-as.numeric(df[,i])/sum(as.numeric(df[,i]),na.rm = T)
      }
    } else {
      print("NormVector Invalid")
    }
  }
  dfout<-df
  return(dfout)
}
