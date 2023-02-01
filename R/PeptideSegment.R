#' PeptideSegment
#'
#' @param filename1 dataset 1
#' @param filename2 dataset 2
#' @param kmin minimum observations needed DEFAULT=2
#' @param simtype similarity type DEFAULT="Tanimoto"
#' @param rel Boolean DEFAULT=TRUE
#' @param OverallName Analysis name string
#' @param SampleName1 Sample name string
#' @param SampleName2 Sample name string
#' @param mn Default=FALSE. In default settings it adjusts the distance scaling by average presence. Setting it to a numeric value will use that as a constant instead. (Recommended 1:2)
#' @param logopt Boolean indicating the use of the log transform before relative scaling of abundance. Default=TRUE
#'
#' @return List of glycopeptides by comparison
#' @export
#'
#' @examples #
PeptideSegment<-function(filename1,filename2,kmin=2,simtype="Tanimoto",rel=TRUE,OverallName,SampleName1,SampleName2,mn=FALSE,logopt=TRUE){
  datfile1<-SimDataCleanLogA(filename1,kmin=kmin,rel = rel)
  datfile2<-SimDataCleanLogA(filename2,kmin=kmin,rel = rel)
  Rnames1<-row.names(datfile1)
  Cleaned1<-gsub("\\s*\\{[^\\}]+\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames1))
  UniCle1<-unique(Cleaned1)
  NameVec1<-rep(0,nrow(datfile1))
  for (j in seq(length(UniCle1))){
    NameVec1[which(Cleaned1 %in% UniCle1[j])]<-j
  }
  Rnames2<-row.names(datfile2)
  Cleaned2<-gsub("\\s*\\{[^\\}]+\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames2))
  UniCle2<-unique(Cleaned2)
  NameVec2<-rep(0,nrow(datfile2))
  for (j in seq(length(UniCle2))){
    NameVec2[which(Cleaned2 %in% UniCle2[j])]<-j
  }
  Output<-list()
  failvec1<-c()
  failvec2<-c()
  if (length(UniCle2)>=length(UniCle1)){
    for (j in seq(length(UniCle2))){
      if (UniCle2[j] %in% UniCle1){
        k<-which(UniCle1 %in% UniCle2[j])
        temp1<-datfile1[which(NameVec1==k),]
        temp2<-datfile2[which(NameVec2==j),]

        Out<-SymmetricalSimBootstrap(temp1,temp2,kmin=1)
        SimPlot(paste0(OverallName,': ',UniCle2[j]),Out)
        append(Output,list(paste0(UniCle2[j]),Out))
        failvec1<-c(failvec1,k)
      } else {
        failvec2<-c(failvec2,j)
      }
    }
    if (length(failvec1)>0){
      print(paste(filename1,'had the following unmatched peptide backbones',paste(UniCle1[-failvec1],collapse=', ')),sep=' ')
    } else {
      print(paste(filename1,'had the following unmatched peptide backbones',paste(UniCle1,collapse=', ')),sep=' ')
    }
    if (length(failvec2)>0){
      print(paste(filename2,'had the following unmatched peptide backbones',paste(UniCle2[failvec2],collapse=', ')),sep=' ')
    } else {
      print(paste(filename2,'had the following unmatched peptide backbones',sep=' '))
    }


  } else {
    for (j in seq(length(UniCle1))){
      if (UniCle1[j] %in% UniCle2){
        k<-which(UniCle2 %in% UniCle1[j])
        temp1<-datfile1[which(NameVec1==j),]
        temp2<-datfile2[which(NameVec2==k),]
        Out<-SymmetricalSimBootstrap(temp1,temp2,kmin=1)
        SimPlot(paste0(OverallName,': ',UniCle2[j]),Out)
        append(Output,list(paste0(UniCle1[j]),Out))
        failvec2<-c(failvec2,k)
      } else {
        failvec1<-c(failvec1,j)
      }
    }
    if (length(failvec2)>0){
      print(paste(filename2,'had the following unmatched peptide backbones',paste(UniCle2[-failvec2],collapse=', ')),sep=' ')
    } else {
      print(paste(filename2,'had the following unmatched peptide backbones',paste(UniCle2,collapse=', ')),sep=' ')
    }
    if (length(failvec1)>0){
      print(paste(filename1,'had the following unmatched peptide backbones',paste(UniCle1[failvec1],collapse=', ')),sep=' ')
    } else {
      print(paste(filename1,'had the following unmatched peptide backbones',sep=' '))
    }
  }
  return(Output)


}
