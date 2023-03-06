#' Standardization_SubFunction scales a dataframe's columns to all be between 0 and 1
#'
#' @param df Dataframe to be standardized
#' @param rel Choice of Standardization method. 'Within' scales compared to largest total signal size in df. 'AsIs' does no scaling. NUMERIC scales to a specified number. Default='Within'
#'
#' @return Standardized dataframe
#' @export
#'
#' @examples #
Standardization_SubFunction<-function(df,rel='Within'){
  if (rel=='Within'){
    scale<-max(colSums(df,na.rm = T))
  } else if (is.numeric(rel)){
    scale<-rel
  } else if (rel=='AsIs'){
    scale<-1
  } else {
    print(paste0('rel invalid'))
  }
  for (i in 1:ncol(df)){
    df[,i]<- as.numeric(df[,i])/scale
  }
  return(df)
}
