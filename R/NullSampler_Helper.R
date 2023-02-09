#' NullSampler_Helper
#' Generates samples from two datasets to create null datasets and then removes duplicates.
#'
#' @param coln1 number of columns (samples) in first dataset
#' @param coln2 number of columns (samples) in second dataset
#' @param sampn initial number of samplings to begin with. Default=200
#'
#' @return List of dataframes of samplings
#' @export
#'
#' @examples #
NULLSAMPLER_Helper<-function(coln1,coln2,sampn=200){
  #get total number of replicates
  colnt<-coln1+coln2
  #check size of replicates
  if (colnt<=6){
    combot1<-data.frame(gtools::combinations(colnt,coln1,repeats.allowed = TRUE))
    combot2<-data.frame(gtools::combinations(colnt,coln2,repeats.allowed = TRUE))
  } else {
    combot1<-data.frame(matrix(1,nrow=sampn,ncol=coln1))
    combot2<-data.frame(matrix(1,nrow=sampn,ncol=coln2))
    for (j in 1:sampn){
      combot1[j,]<-sort(sample(colnt,coln1,replace = TRUE))
      combot2[j,]<-sort(sample(colnt,coln2,replace = TRUE))
    }
    #remove duplicates
    if (sum(duplicated(combot1))>0){
      combot1<-combot1[-which(duplicated(combot1)),]
    }
    if (sum(duplicated(combot2))>0){
      combot2<-combot2[-which(duplicated(combot2)),]
    }
  }
  row.names(combot2)<-NULL
  row.names(combot1)<-NULL

  #remove all with less than 2 of each side
  if (coln1>3 & coln2>3){
    minsample<-2
  } else {
    minsample<-1
  }
  #check no overlap between null distributions
  if (coln1==coln2){
    combot<-rbind(combot1,combot2,all=T,fill=T)
    if (sum(duplicated(combot))>0){
      combot<-combot[-which(duplicated(combot)),]
    }
    combot<-SAMPLER_Helper_Parity(combot,coln1,coln2,minsample)
    si<-sample(dim(combot)[1],floor(dim(combot)[1]/2))
    combot1<-combot[si,]
    combot2<-combot[-si,]
    row.names(combot2)<-NULL
    row.names(combot1)<-NULL
  } else {
    combot1<-SAMPLER_Helper_Parity(combot1,coln1,coln2,minsample)
    combot2<-SAMPLER_Helper_Parity(combot2,coln1,coln2,minsample)
  }
  null_samples<-c()
  null_samples$NDis1<-combot1
  null_samples$NDis2<-combot2
  return(null_samples)
}

#' Sampler_Helper_Parity is used by Null Sampler Helper function to ensure that the datasets are represented equally
#'
#' @param combo samplings
#' @param coln1 samples in 1
#' @param coln2 samples in 2
#' @param minsample minimum number of samples needed from each
#'
#' @return
#' @export
#'
#' @examples #
SAMPLER_Helper_Parity<-function(combo,coln1,coln2,minsample=2){
  # any(rowSums(combo<=coln1)<minsample) identifies
  ## i. samples from file1: <=
  ## ii. how many of those samples are in a row: rowSums(EXPi)
  ## iii. if a row has less than 2 (or 1 for low sample): any(EXPii<minsample)
  # any(rowSums(combo>coln1)<minsample)
  ## same as above but for samples from file2
  if (any(rowSums(combo<=coln1)<minsample) | any(rowSums(combo>coln1)<minsample)){
    if (any(rowSums(combo<=coln1)<minsample)){
      combo<-combo[-which(rowSums(combo<=coln1)<minsample),]
      row.names(combo)<-NULL
    }
    if (sum(rowSums(combo>coln1)<minsample)>0){
      combo<-combo[-which(rowSums(combo>coln1)<minsample),]
      row.names(combo)<-NULL
    }
  }
  return(combo)
}

