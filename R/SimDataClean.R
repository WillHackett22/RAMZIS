#' SimDataClean Standardize and normalize an abundance matrix
#'
#' @param filename The abundance matrix in question. This can also be the matrix/df already in system
#' @param kmin The minimum number of observations needed for a glycopeptide to be real. Default=2
#' @param rel Pick if values should be normalized relative to the max total abundance of any column(TRUE) or read as is(FALSE). Default = TRUE
#' @param normvector A vector for normalizing signal between samples, usually based on the XIC. Default: 'None' reads values as is; 'Relative' scales samples to their overall signal (assuming equal signal distribution). Providing a vector multiplies sample glycopeptides by that vector
#' @param logoption Boolean indicating use of log transformation on values. This occurs after normalization and before relativity. Default= TRUE
#'
#' @return The normalized and standardized relative log abundance dataframe.
#' @export
#'
#' @examples
#' #From the outputs of GlycReRead example
#' #Sample1DataMatrix<-SimDataClean('OutputSaveFile.csv')
#' ## OR
#' #Sample1DataMatric<-SimDataClean(AbundanceDF)
SimDataClean<-function(filename,kmin=2,rel=TRUE,normvector='None',logoption=TRUE){
  #Read in Data and measure dataframe size
  if (typeof(filename)=='character'){
    file1<-read.csv(filename,header=TRUE, row.names=1,stringsAsFactors = FALSE)
    #normalize based off of normalization vector
    if (logoption){
      file1<-log(file1+1)
    }
    if (length(normvector)==ncol(file1) ){
      for (i in 1:ncol(file1)){
        file1[,i]<-as.numeric(file1[,i])*normvector[i]
      }
    } else {
      if (normvector=='None'){
        data1<-file1
      } else if (normvector=='Relative'){
        for (i in 1:ncol(file1)){
          file1[,i]<-as.numeric(file1[,i])/sum(as.numeric(file1[,i]),na.rm = T)
        }
      } else {
        print("NormVector Invalid")
      }
    }
    data1<-file1
    #standardize data by max abundance
    if (rel){
      for (i in 1:ncol(file1)){
        data1[,i]<- as.numeric(file1[,i])/max(colSums(file1,na.rm = T))
      }
    }
  } else if(typeof(filename)=='list'){
    file1<-filename
    if (logoption){
      file1<-log(file1+1)
    }
    if (length(normvector)==ncol(file1) ){
      for (i in 1:ncol(file1)){
        file1[,i]<-as.numeric(file1[,i])*normvector[i]
      }
    } else {
      if (normvector=='None'){
        data1<-file1
      } else if (normvector=='Relative'){
        for (i in 1:ncol(file1)){
          file1[,i]<-as.numeric(file1[,i])/sum(as.numeric(file1[,i]),na.rm = T)
        }
      }
    }
    data1<-file1
    if (rel){
      for (i in 1:ncol(file1)){
        data1[,i]<- as.numeric(file1[,i])/max(colSums(file1,na.rm = T))
      }
    }
  }
  file1[is.na(file1)]<-0
  data1[is.na(data1)]<-0


  rown<-dim(file1)[1] #GPs
  coln<-dim(file1)[2] #samples
  #Clean Dataframes: remove rows where there are one or fewer signals
  #number of replicates
  reps<- coln
  num<-0
  #count the missing values in rows
  mvhold<-rep(0,rown)
  for (row in 1:rown){
    rowi<-row-num
    mvcount<-0
    #check how many missing values in a given row
    for (col in 1:coln){
      for (x in file1[row,col]){
        if (x==0){
          mvcount<- mvcount+1
        }
      }
    }
    mvhold[row]<-mvcount
  }
  #turn missing values into observed values
  remhold<-coln-mvhold
  #remove data where less than kmin (default=2) are seen
  if (length(which(remhold<kmin))>0){
    data1<-data1[-which(remhold<kmin),]
  }
  #remove data where mean is lower than logx=-10
  if (sum(-10>log(apply(data1,1,mean)))>0){
    data1<-data1[-which(-10>log(apply(data1,1,mean))),]
  }
  return(data1)
}
