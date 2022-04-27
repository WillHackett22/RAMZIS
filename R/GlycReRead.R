GlycReRead<-function(Filelist,Outfilename=NULL,verbose=T){
  if (length(Filelist)==1){
    GFiles<-read.table(Filelist,stringsAsFactors = FALSE)
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
  GFile<-read.csv(FileG,header=TRUE,stringsAsFactors = FALSE)
  if (sum(GFile$total_signal==0)>0){
    GFile<-GFile[-which(GFile$total_signal==0),]
  }
  key<-paste(GFile$protein_name,'::',GFile$glycopeptide)
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
    tempg<-data.frame(matrix(nrow=length(key),ncol=1))
    row.names(tempg)<-key
    colnames(tempg)<-FNameB
    tempg[,FNameB]<-abun
    tempdata<-merge(tempdata,tempg,by=0,all = TRUE)
    row.names(tempdata)<-tempdata[,1]
    tempdata<-tempdata[,-1]
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
      write.csv(protdata,paste0(Outfilename,'_',name,'.csv'))
    }
    if (verbose){
      name<-strsplit(uniprot[j],"\\|")[[1]][length(strsplit(uniprot[j],"\\|")[[1]])]
      name<-strsplit(name,' ')[[1]][1]
      outlist[[name]]<-protdata
    }
  }
  unigp<-duplicated(glypep)
  if (sum(unigp)>0){
    tempdata<-tempdata[-which(unigp),]
  }
  gpeprow<-vapply(strsplit(row.names(tempdata),'::'), `[`, 2, FUN.VALUE=character(1))
  row.names(tempdata)<-gpeprow

  if (!is.null(Outfilename)){
    write.csv(tempdata,paste0(Outfilename,'.csv'))
  }
  if (verbose){
    return(outlist)
  }
}
