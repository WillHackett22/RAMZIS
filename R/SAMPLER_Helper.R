#' SAMPLER_Helper
#' Generates samples of a single data set and removes duplicates
#'
#' @param coln Number of columns (samples) in a dataset
#' @param cols Number of rows (glycopeptides) in a dataset
#' @param sampn Number of samplings desired
#' @param filename Source of dataset
#'
#' @return
#' @export
#'
#' @examples
SAMPLER_Helper<-function(coln,cols,sampn=100,filename=NULL){
  if (coln<6 & 2<coln){
    Nsmpl<-CombWRep(coln,coln-1)
    combo<-data.frame(matrix(1,nrow=Nsmpl,ncol=(coln)))
    combo<-data.frame(gtools::combinations(coln,(coln),repeats.allowed = TRUE))
  } else if (coln<=2) {
    if (!is.null(filename)){
      print(paste(filename,' has too few samples. This will not work.'))
    }
    stop()
  } else {
    combo<-data.frame(matrix(1,nrow=sampn,ncol=(coln)))
    for (j in 1:100){
      combo[j,]<-sort(sample(coln1,(coln),replace = TRUE))
    }
    if (sum(duplicated(combo))>0){
      combo<-combo[-which(duplicated(combo)),]
    }
  }
  return(combo)
}
