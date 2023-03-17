#' SAMPLER_HelperTIDY
#' Generates samples of a single data set and removes duplicates
#'
#' @param coln Number of columns (samples) in a dataset
#' @param sampn Number of samplings desired
#' @param filename Source of dataset
#' @param GroupID String identifying which group this applies to. ie TestA, TestB, NullA, NullB. Defaults to "Group"
#'
#' @return Matrix of indices for samples
#' @export
#'
#' @examples #
SAMPLER_HelperTIDY<-function(coln,sampn=100,filename=NULL,GroupID="Group"){
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
      combo[j,]<-sort(sample(coln,(coln),replace = TRUE))
    }
    if (sum(duplicated(combo))>0){
      combo<-combo[-which(duplicated(combo)),]
    }
  }
  row.names(combo)<-paste0(1:nrow(combo),"_",GroupID)
  return(combo)
}
