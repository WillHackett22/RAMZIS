#' SimDataClean Standardize and normalize an abundance matrix
#'
#' @param filename The abundance matrix in question. This can also be the matrix/df already in system
#' @param kmin The minimum number of observations needed for a glycopeptide to be real. Default=2
#' @param rel Pick if values should be normalized relative to the max total abundance of any column('Within'), to a given value (NUMERIC), or read as is('AsIs'). Recommended numeric is the max of both columns for a given value. Default = 'Within'. Only runs if no normvector given, or if rel_force=T
#' @param normvector A vector for normalizing signal between samples, usually based on the XIC. Default: 'None' reads values as is; 'Relative' scales samples to their overall signal (assuming equal signal distribution). Providing a vector multiplies sample glycopeptides by that vector
#' @param logoption Boolean indicating use of log transformation on values. This occurs after normalization and before relativity. Default= TRUE
#' @param rel_force Default=FALSE. If true, forces standardization by total signal even when given normvector.
#'
#' @return The normalized and standardized relative log abundance dataframe.
#' @export
#'
#' @examples
#' #From the outputs of GlycReRead example
#' #Sample1DataMatrix<-SimDataClean('OutputSaveFile.csv')
#' ## OR
#' #Sample1DataMatric<-SimDataClean(AbundanceDF)
SimDataClean<-function(filename,kmin=2,rel='Within',normvector='None',logoption=TRUE,rel_force=FALSE){
  #Read in Data and measure dataframe size
  if (typeof(filename)=='character'){
    file1<-read.csv(filename,header=TRUE, row.names=1,stringsAsFactors = FALSE)
  } else if (typeof(filename)=='list'){
    file1<-filename
  }
  file1[is.na(file1)]<-0
  #normalize based off of normalization vector
  data1<-Normalization_SubFunction(file1,normvector,logoption)
  #Clean Dataframes: remove rows where there are one or fewer signals
  #number of replicates
  coln<-dim(data1)[2] #samples
  mvhold<-apply(data1,1,MVCount)
  remhold<-coln-mvhold
  #remove data where less than kmin (default=2) are seen
  if (length(which(remhold<kmin))>0){
    data1<-data1[-which(remhold<kmin),]
  }
  #standardize data by max abundance
  if (normvector!="None" & rel_force==FALSE){
    rel<-"AsIs"
  }
  data1<-Standardization_SubFunction(data1,rel)
  data1[is.na(data1)]<-0
  #remove data where mean is lower than logx=-10
  if (sum(-10>log(apply(data1,1,mean)))>0){
    data1<-data1[-which(-10>log(apply(data1,1,mean))),]
  }
  return(data1)
}
