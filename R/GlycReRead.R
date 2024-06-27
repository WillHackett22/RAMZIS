#' GlycReRead: Read a set of GlycReSoft files and collate them into a joint dataset
#'
#' @param Filelist A vector of exported GlycReSoft results files (.csv) of at least length 1.
#' @param Outfilename The save location for the collective dataset (don't include a filetype). If unused results will not be saved.
#' @param verbose Default: TRUE  Returns dataset to environment
#' @param JoinDuplicates Default: TRUE. Joins duplicate glycopeptides by summation on total_signals
#' @param CleanProtName Default: TRUE. Turns all punctuation in protein names into '_'
#' @param SaveResults Default: TRUE. If FALSE, disables csv saving.
#'
#' @return Matrix of glycopeptide abundances with glycopeptides as rows and samples as columns
#' @export
#' @import utils
#'
#' @examples
#' #files<-c('gsoft1.csv','gsoft2.csv','gsoft3.csv')
#' #AbundanceDF<-GlycReRead(files,'OutputSaveFile')
GlycReRead<-function(Filelist,Outfilename=NULL,verbose=T,JoinDuplicates=T,CleanProtName=T,SaveResults=T){
  if (length(Filelist)==1){
    GFiles<-utils::read.table(Filelist,stringsAsFactors = FALSE)
    GFiles<-GFiles[,1]
  } else {
    GFiles<-Filelist
  }
  GLen<-length(GFiles)
  NameHold<-rep(0,GLen)
  HoldTable<-data.frame(matrix(0,nrow = 1,ncol=(GLen+1)))
  FileG<-GFiles[1]
  FNameA<-strsplit(FileG,'\\/')[[1]][length(strsplit(FileG,'\\/')[[1]])]
  FNameB<-strsplit(FNameA,'\\.csv')[[1]][length(strsplit(FNameA,'\\.csv')[[1]])]
  NameHold[1]<-FNameB
  GFile<-utils::read.csv(FileG,header=TRUE,stringsAsFactors = FALSE)
  if (sum(GFile$total_signal==0)>0){
    GFile<-GFile[-which(GFile$total_signal==0),]
  }
  key<-paste(GFile$protein_name,'::',GFile$glycopeptide)
  if (any(duplicated(key))){
    if (JoinDuplicates==T){
      while(any(duplicated(key))){
        dx<-which(duplicated(key))[1]
        kx<-which(key==key[1])[1]
        GFile$total_signal[kx]<-GFile$total_signal[kx]+GFile$total_signal[dx]
        GFile<-GFile[-dx,]
        key<-paste(GFile$protein_name,'::',GFile$glycopeptide)
      }
    } else {
      print("Duplicate Identifications detected. Please merge or rename.")
    }
  }
  abun<-GFile$total_signal
  tempdata<-data.frame(matrix(nrow=length(key),ncol=1))
  row.names(tempdata)<-key
  colnames(tempdata)<-FNameB
  tempdata[,FNameB]<-abun
  for (j in 2:GLen){
    FileG<-GFiles[j]
    FNameA<-strsplit(FileG,'\\/')[[1]][length(strsplit(FileG,'\\/')[[1]])]
    FNameB<-strsplit(FNameA,'\\.csv')[[1]][length(strsplit(FNameA,'\\.csv')[[1]])]
    NameHold[j]<-FNameB
    GFile<-read.csv(FileG,header=TRUE,stringsAsFactors = FALSE)
    if (sum(GFile$total_signal==0)>0){
      GFile<-GFile[-which(GFile$total_signal==0),]
    }
    key<-paste(GFile$protein_name,'::',GFile$glycopeptide)
    abun<-GFile$total_signal
    while (any(duplicated(key))){
      kix<-which(duplicated(key))[1]
      gix<-which(key[kix]==key)[1]
      abun[gix]<-abun[gix]+abun[kix]
      key<-key[-kix]
      abun<-abun[-kix]
    }
    tempg<-data.frame(matrix(nrow=length(key),ncol=1))
    row.names(tempg)<-key
    colnames(tempg)<-FNameB
    tempg[,FNameB]<-abun
    tempdata<-MatrixMerge_Helper(tempdata,tempg)
  }
  glyprot<-row.names(tempdata)
  prot<-vapply(strsplit(glyprot,'::'), `[`, 1, FUN.VALUE=character(1))
  glypep<-vapply(strsplit(glyprot,'::'), `[`, 2, FUN.VALUE=character(1))
  uniprot<-unique(prot)
  outlist<-list('Full'=tempdata)
  for (j in 1:length(uniprot)){
    idx<-which(prot==uniprot[j])
    protdata<-tempdata[idx,]
    row.names(protdata)<-vapply(strsplit(row.names(protdata),'::'), `[`, 2, FUN.VALUE=character(1))
    if (!is.null(Outfilename)){
      name<-strsplit(uniprot[j],"\\|")[[1]][length(strsplit(uniprot[j],"\\|")[[1]])]
      name<-strsplit(name,' ')[[1]][1]
      if (CleanProtName){
        name<-gsub("[[:punct:]]","_",name)
      }
      utils::write.csv(protdata,paste0(Outfilename,'_',name,'.csv'))
    }
    if (verbose){
      name<-strsplit(uniprot[j],"\\|")[[1]][length(strsplit(uniprot[j],"\\|")[[1]])]
      name<-strsplit(name,' ')[[1]][1]
      if (CleanProtName){
        name<-gsub("[[:punct:]]","_",name)
      }
      outlist[[name]]<-protdata
    }
  }
  unigp<-duplicated(glypep)
  if (sum(unigp)>0){
    tempdata<-tempdata[-which(unigp),]
  }
  gpeprow<-vapply(strsplit(row.names(tempdata),'::'), `[`, 2, FUN.VALUE=character(1))
  row.names(tempdata)<-gpeprow
  if (SaveResults){
    if (!is.null(Outfilename)){
      utils::write.csv(tempdata,paste0(Outfilename,'.csv'))
    }
  }
  if (verbose){
    return(outlist)
  }
}
