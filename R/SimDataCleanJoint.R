#' SimDataCleanJoint normalizes and standardizes both data sets. By default it standardizes to the joint max signal.
#'
#' @param filename1 First file or dataframe: sample 1
#' @param filename2 Second file or dataframe: sample 2
#' @param kmin Minimum number of identifications needed to be considered real. Default=2
#' @param rel Choice of Standardization method. 'Within' scales compared to largest signal size in each sample. 'Joint' scales to largest signal value between both files. 'AsIs' does no scaling. NUMERIC scales to a specified number. Default='Joint'
#' @param normvector List of Normalization vectors to be multiplied against. Vector length should equal sample size. Default=list('None','None')
#' @param logoption Boolean indicating use of log transformation. Default=TRUE
#'
#' @return Normalized and standardized cleaned data
#' @export
#'
#' @examples #
SimDataCleanJoint<-function(filename1,filename2,kmin=2,rel='Joint',normvector=list('None','None'),logoption=TRUE){
  if (typeof(filename1)=='character'){
    file1<-read.csv(filename1,header=TRUE, row.names=1,stringsAsFactors = FALSE)
  } else if (typeof(filename1)=='list'){
    file1<-filename1
  }
  d1temp<-Normalization_SubFunction(file1,normvector[[1]],logoption)
  if (typeof(filename2)=='character'){
    file2<-read.csv(filename2,header=TRUE, row.names=1,stringsAsFactors = FALSE)
  } else if (typeof(filename2)=='list'){
    file2<-filename2
  }
  d2temp<-Normalization_SubFunction(file2,normvector[[2]],logoption)
  dj<-MatrixMerge_Helper(d1temp,d2temp)
  if (rel=='Joint'){
    rel<-max(colSums(dj,na.rm=T))
  }
  data1<-SimDataClean(filename1,kmin,rel,normvector[[1]],logoption)
  data2<-SimDataClean(filename2,kmin,rel,normvector[[2]],logoption)
  dataout<-list("DF1"=data1,"DF2"=data2)
  return(dataout)
}
