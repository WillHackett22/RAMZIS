#RAMZIS Draft Script2
#Will Hackett
#1-14-19
#Functions for Ranking Glycopeptides based off of comparative similarity metrics
#Final Function:
#RAMZIS(FileList1,FileList2,Norm1,Norm2,kmin)
#If your data is pre-normalized: Use a lower level function
#FileList1/2: List of GlycReSoft output files
#Norm1/2: Normalization Factor pulled from PEAKS analysis, defaults to all being equal
#kmin= minimum number of observations a glycopeptide  must have in a sample group to be observed
##
##Previous Expected Input
##Expected input: CSV files that have columns of total signal values from glycresoft output
##Expected input Cont: Files should be glycopeptides for a single protein
##Expected input cont: eg File1 and File2 should only have values from AGP1 to compare AGP1
##Expected input cont: The first column should be the glycopeptide, all others are signal columns from different samples
#
#RAMZIS Function Structure
#
#
#Functions and Dependencies:
#WithinSim
##Takes a single glycopeptide-signal filename (with location) as input. 
##Outputs within similarity distribution
#BetweenSim
##Takes two glycopeptide-signal filenames as input.
##Outputs similarity value comparing the two files
#SimPlot
##Takes two within distributions and a between similarity value
##Outputs a comparative plot of the two distributions
#colMax
##Finds maxima by column
#remrow
##removes rows from dataframe
#SimDataClean
##Loads and Cleans a dataframe
library(MASS)
library(gtools)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
remrow <- function(x, rows) x[-rows,, drop = FALSE]
remcol <- function(x, cols) x[,-cols , drop = FALSE]

#GlycReRead
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
      name<-strsplit(uniprot[j],"\\|")[[1]][2]
      name<-strsplit(name,' ')[[1]][1]
      write.csv(protdata,paste0(Outfilename,'_',name,'.csv'))
    }
    if (verbose){
      name<-strsplit(uniprot[j],"\\|")[[1]][length(strsplit(uniprot[j],"\\|")[[1]])]
      name<-strsplit(name,' ')[[1]][1]
      
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
    return(tempdata)
  }
}

#not yet functional
ByonicRead<-function(Filelist,Outfilename=NULL,verbose=T){
  print('This function is not yet functional')
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
  GFile<-GFile[-which(GFile$total_signal==0),]
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
    GFile<-GFile[-which(GFile$total_signal==0),]
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
  if (!is.null(Outfilename)){
    write.csv(tempdata,Outfilename)
  }
  if (verbose){
    return(tempdata)
  }
}


GlyLineRead<-function(Filelist){}

#SimDataClean
SimDataClean<-function(filename,kmin=2,rel=TRUE,normvector='Default',logopt=FALSE){
  #Read in Data and measure dataframe size
  if (typeof(filename)=='character'){
    file1<-read.csv(filename,header=TRUE, row.names=1,stringsAsFactors = FALSE)
    #normalize based off of normalization vector
    if (length(normvector)==ncol(file1) ){
      for (i in 1:ncol(file1)){
        file1[,i]<-as.numeric(file1[,i])*normvector[i]
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
    if (length(normvector)==ncol(file1) ){
      for (i in 1:ncol(file1)){
        file1[,i]<-as.numeric(file1[,i])*normvector[i]
      }
    }
    data1<-file1
    if (rel){
      for (i in 1:ncol(file1)){
        data1[,i]<- as.numeric(file1[,i])/max(colSums(file1,na.rm = T))
      }
    }
    
    #print('Program Detected List for filename. Attempting use as dataframe')
  }
  file1[is.na(file1)]<-0
  data1[is.na(data1)]<-0
  
  if (logopt){
    data1<-log(data1)
  }
  
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
    data1<-data1[-which(remhold<2),]
  }
  return(data1)
}

SimDataCleanLog<-function(filename,kmin=2,rel=TRUE,normvector='Default'){
  #Read in Data and measure dataframe size
  if (typeof(filename)=='character'){
    file1<-read.csv(filename,header=TRUE, row.names=1,stringsAsFactors = FALSE)
    #normalize based off of normalization vector
    if (length(normvector)==ncol(file1) ){
      for (i in 1:ncol(file1)){
        file1[,i]<-as.numeric(file1[,i])*normvector[i]
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
    if (length(normvector)==ncol(file1) ){
      for (i in 1:ncol(file1)){
        file1[,i]<-as.numeric(file1[,i])*normvector[i]
      }
    }
    data1<-file1
    if (rel){
      for (i in 1:ncol(file1)){
        data1[,i]<- as.numeric(file1[,i])/max(colSums(file1,na.rm = T))
      }
    }
    
    #print('Program Detected List for filename. Attempting use as dataframe')
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

#SimDataCleanLog that does individual relative abundance
SimDataCleanLogA<-function(filename,kmin=2,rel=TRUE,normvector='None'){
  #Read in Data and measure dataframe size
  if (typeof(filename)=='character'){
    file1<-read.csv(filename,header=TRUE, row.names=1,stringsAsFactors = FALSE)
    #normalize based off of normalization vector
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
    #standardize data by max abundance 
    if (rel){
      for (i in 1:ncol(file1)){
        data1[,i]<- as.numeric(file1[,i])/max(colSums(file1,na.rm = T))
      }
    }
  } else if(typeof(filename)=='list'){
    file1<-filename
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
    
    #print('Program Detected List for filename. Attempting use as dataframe')
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


PeptideCollapser<-function(filename,kmin=2,rel=TRUE,sequonvector){
  datfile1<-SimDataClean(filename,kmin,rel)
  #find sequon sections
  #find glycans associated with a sequon
  #merge same glycans for a given sequon
  #return datafile with sequon plus glycan
  Rnames<-row.names(datfile1)
  Cleaned<-gsub("\\s*\\{[^\\}]+\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames))
  Outdat<-data.frame(matrix(NA,nrow=1,ncol=dim(datfile1)[2]))
  idxterm<-0
  NewRnameList<-c()
  for (j in 1:length(sequonvector)){
    Selection<-grep(sequonvector[j],Cleaned)
    Glycans<-gsub(".*\\{|\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames[Selection]))
    UniGly<-unique(Glycans)
    for (l in 1:length(UniGly)){
      RowGlyPep<-paste0(sequonvector[j],'{',UniGly[l],'}')
      GlySel<-grep(UniGly[l],Glycans)
      Outdat[idxterm+l,]<-colSums(datfile1[Selection,][GlySel,],na.rm=T)
      NewRnameList<-c(NewRnameList,RowGlyPep)
    }
    idxterm<-idxterm+length(UniGly)
    
  }
  row.names(Outdat)<-NewRnameList
  return(Outdat)
}

#Presence Ratio Function
GPPresence<-function(dataf,kmin=1,rel=TRUE){
  file1<-dataf
  jlen<-dim(file1)[1] # glycopeptides
  jtot<-dim(file1)[2] # samples
  if (jlen>0){
    phold<-rep(0,jlen) # glycopeptides presence holder
    for (j in seq(jlen)){ #iterate through GPs
      phold[j]<-1-sum(file1[j,]==0)/jtot #1 -number missing/number samples
    }
  } else {
    phold<-0
  }
  
  return(phold)
}

#Modality Test
Modality<-function(SimDist,threshold=.01){
  SimDens<-density(SimDist)
  #build sim density distribution
  LocalMinimaIdx<-which(diff(diff(SimDens$y)<=0)<0)+1
  LocalMaximaIdx<-which(diff(diff(SimDens$y)>=0)<0)+1
  #find probable minima and maxima
  MinMat<-data.frame(matrix(NA,nrow=length(LocalMaximaIdx),ncol=9))
  colnames(MinMat)<-c('Peak_x','dPeakLeft','dPeakRight','LeftCheck','RightCheck','Overall','N_Members','InterPeakDistance','Peak_y')
  MinMat[,1]<-SimDens$x[LocalMaximaIdx]
  MinMat[,9]<-SimDens$y[LocalMaximaIdx]
  #Peak Locations
  MinMat[,2]<-(SimDens$y[LocalMaximaIdx]-c(0,SimDens$y[c(LocalMinimaIdx)]))
  #Height Differentials before peak
  MinMat[,3]<-(SimDens$y[LocalMaximaIdx]-c(SimDens$y[c(LocalMinimaIdx)],0))
  #height differentials after peak
  MinMat[,4]<-abs(MinMat[,2])<=threshold*SimDens$y[LocalMaximaIdx]
  MinMat[,5]<-abs(MinMat[,3])<=threshold*SimDens$y[LocalMaximaIdx]
  #chance that the height differentials are by random chance
  
  TempMin<-c(0,SimDens$x[LocalMinimaIdx],1)
  for (j in 1:length(LocalMaximaIdx)){
    temppeak<-which(TempMin[j]<=SimDist & SimDist<=TempMin[j+1])
    MinMat[j,7]<-length(temppeak)
    #XSD<-sd(SimDist[temppeak])
    #MinMat[j,8]<-abs(TempMin[j+1]-TempMin[j])<=3*XSD
    MinMat[j,8]<-FALSE
  }
  MinMat[,6]<-(sum(MinMat[,4])>0|sum(MinMat[,5])>0)
  
  while (sum(MinMat[,6])>0){
    rowpick<-which(MinMat[,4]|MinMat[,5])
    targeti<-which(rowSums(abs(MinMat[rowpick,c(2,3)])==min(abs(MinMat[rowpick,c(2,3)])))>0)[1]
    targeti<-rowpick[targeti]
    #find which peaks have the smallest differential, ties go to the first
    targetl<-which(abs(MinMat[targeti,c(2,3)])==min(abs(MinMat[targeti,c(2,3)])))[1]-2
    #find whether it is the preceding or secondary trough, ties go to the preceding
    #if first peak and trough, remove peak and next trough
    #if last peak and last trough, remove peak and prior trough
    #else remove peak and lowest trough
    if (targetl==-1 & targeti==1){
      LocalMaximaIdx<-LocalMaximaIdx[-targeti]
      LocalMinimaIdx<-LocalMinimaIdx[-1]
      #remove peak
    } else if (targetl==0 & targeti==length(LocalMaximaIdx)){
      LocalMaximaIdx<-LocalMaximaIdx[-targeti]
      LocalMinimaIdx<-LocalMinimaIdx[-length(LocalMinimaIdx)]
      #remove peak
    } else {
      LocalMaximaIdx<-LocalMaximaIdx[-targeti]
      #remove peak
      tempidx<-targeti+targetl
      LocalMinimaIdx<-LocalMinimaIdx[-tempidx]
      #remove lower trough
    }
    if (length(LocalMinimaIdx)>0){
      MinMat<-data.frame(matrix(NA,nrow=length(LocalMaximaIdx),ncol=9))
      colnames(MinMat)<-c('Peak_x','dPeakLeft','dPeakRight','LeftCheck','RightCheck','Overall','N_Members','InterPeakDistance','Peak_y')
      MinMat[,1]<-SimDens$x[LocalMaximaIdx]
      MinMat[,9]<-SimDens$y[LocalMaximaIdx]
      #Peak Locations
      MinMat[,2]<-(SimDens$y[LocalMaximaIdx]-c(0,SimDens$y[c(LocalMinimaIdx)]))
      #Height Differentials before peak
      MinMat[,3]<-(SimDens$y[LocalMaximaIdx]-c(SimDens$y[c(LocalMinimaIdx)],0))
      #height differentials after peak
      MinMat[,4]<-abs(MinMat[,2])<=threshold*SimDens$y[LocalMaximaIdx]
      MinMat[,5]<-abs(MinMat[,3])<=threshold*SimDens$y[LocalMaximaIdx]
      #chance that the height differentials are by random chance
      
      TempMin<-c(0,SimDens$x[LocalMinimaIdx],1)
      for (j in 1:length(LocalMaximaIdx)){
        temppeak<-which(TempMin[j]<=SimDist & SimDist<=TempMin[j+1])
        MinMat[j,7]<-length(temppeak)
        #XSD<-sd(SimDist[temppeak])
        #MinMat[j,8]<-abs(TempMin[j+1]-TempMin[j])<=3*XSD
        MinMat[j,8]<-FALSE
      }
      MinMat[,6]<-(sum(MinMat[,4])>0|sum(MinMat[,5])>0)
    } else {
      MinMat<-data.frame(matrix(NA,nrow=length(LocalMaximaIdx),ncol=9))
      colnames(MinMat)<-c('Peak_x','dPeakLeft','dPeakRight','LeftCheck','RightCheck','Overall','N_Members','InterPeakDistance','Peak_y')
      MinMat[,1]<-SimDens$x[LocalMaximaIdx]
      MinMat[,9]<-SimDens$y[LocalMaximaIdx]
      #Peak Locations
      MinMat[,2]<-(SimDens$y[LocalMaximaIdx]-c(0,SimDens$y[c(LocalMinimaIdx)]))
      #Height Differentials before peak
      MinMat[,3]<-(SimDens$y[LocalMaximaIdx]-c(SimDens$y[c(LocalMinimaIdx)],0))
      #height differentials after peak
      MinMat[,4]<-FALSE
      MinMat[,5]<-FALSE
      #chance that the height differentials are by random chance
      
      TempMin<-c(0,SimDens$x[LocalMinimaIdx],1)
      for (j in 1:length(LocalMaximaIdx)){
        temppeak<-which(TempMin[j]<=SimDist & SimDist<=TempMin[j+1])
        MinMat[j,7]<-length(temppeak)
        #XSD<-sd(SimDist[temppeak])
        #MinMat[j,8]<-abs(TempMin[j+1]-TempMin[j])<=3*XSD
        MinMat[j,8]<-FALSE
      }
      MinMat[,6]<-(sum(MinMat[,4])>0|sum(MinMat[,5])>0)
    }
    
    
  }
  return(list(MinMat,'Maxima'=LocalMaximaIdx,'Minima'=LocalMinimaIdx))
}

#membership of Null Distribution
NullMembershipProportion<-function(SimDist,BootDis1,BootDis2,ModalityList){
  BootCombinations<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=2))
  BootCombinations[,1]<-rep(seq(1,dim(BootDis1)[1]),each=dim(BootDis2)[1])
  BootCombinations[,2]<-rep(seq(1,dim(BootDis2)[1]),dim(BootDis1)[1])
  #make density
  SimDens<-density(SimDist)
  #find minima
  LocalMinimaIdx<-ModalityList$Minima
  idxL<-length(LocalMinimaIdx)+1
  LBound<-c(0,SimDens$x[LocalMinimaIdx])
  UBound<-c(SimDens$x[LocalMinimaIdx],1)
  #determine the membership of each comparison
  BootMembers<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=dim(BootDis1)[2]+dim(BootDis2)[2]))
  MemberProp<-data.frame(matrix(NA,nrow=dim(BootCombinations)[1],ncol=3))
  MemberProp[,1]<-rowSums(BootDis1[BootCombinations[,1],]<=dim(BootDis1)[2])/dim(BootDis1)[2]
  MemberProp[,2]<-rowSums(BootDis2[BootCombinations[,2],]<=dim(BootDis2)[2])/dim(BootDis2)[2]
  MemberProp[,3]<-MemberProp[,1]*(1-MemberProp[,2])
  for (j in 1:dim(BootMembers)[2]){
    BootMembers[,j]<-rowSums(BootDis1[BootCombinations[,1],]==j)+rowSums(BootDis2[BootCombinations[,2],]==j)
  }
  
  #find members by 
  Members<-data.frame(matrix(NA,nrow=idxL,ncol=dim(BootDis1)[2]+dim(BootDis2)[2]))
  MemberPropCount<-data.frame(matrix(NA,nrow=idxL,ncol=length(unique(unlist(MemberProp)))))
  PropList<-unique(unlist(MemberProp))
  colnames(MemberPropCount)<-PropList
  MemberDegree<-data.frame(matrix(NA,nrow=idxL,ncol=ceiling(dim(BootDis1)[2]+dim(BootDis2)[2])/2))
  colnames(MemberDegree)<-paste0(seq(1,dim(MemberDegree)[2])/dim(MemberDegree)[2])
  for (j in 1:idxL){
    Members[j,]<-colSums(BootMembers[which((SimDist>=LBound[j]&SimDist<=UBound[j])),])
    Origin<-rowSums(BootMembers[which((SimDist>=LBound[j]&SimDist<=UBound[j])),1:dim(MemberDegree)[2]])/(dim(BootDis1)[2]+dim(BootDis2)[2])
    MemberDegree[j,1]<-sum(Origin<=1/dim(MemberDegree)[2])
    
    for (l in 1:dim(MemberPropCount)[2]){
      MemberPropCount[j,l]<-sum(MemberProp[which((SimDist>=LBound[j]&SimDist<=UBound[j])),]==PropList[l])
    }
    for (l in 2:dim(MemberDegree)[2]){
      MemberDegree[j,l]<-sum(Origin<=(l/dim(MemberDegree)[2])&Origin>((l-1)/dim(MemberDegree)[2]))/(length(Origin)*idxL)
    }
  }
  MemberProp<-t(t(Members)/colSums(Members))
  return(MemberProp)
}

#internal Membership
InternalMembershipProportion<-function(IntDist,BootDis1,ModalityList){
  BootCombinations<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis1)[1],ncol=2))
  BootCombinations[,1]<-rep(seq(1,dim(BootDis1)[1]),each=dim(BootDis1)[1])
  BootCombinations[,2]<-rep(seq(1,dim(BootDis1)[1]),dim(BootDis1)[1])
  BootCombinations<-BootCombinations[-which(BootCombinations[,1]==BootCombinations[,2]),]
  
  #make density
  SimDens<-density(IntDist)
  #find minima
  LocalMinimaIdx<-ModalityList$Minima
  idxL<-length(LocalMinimaIdx)+1
  LBound<-c(0,SimDens$x[LocalMinimaIdx])
  UBound<-c(SimDens$x[LocalMinimaIdx],1)
  #determine the membership of each comparison
  BootMembers<-data.frame(matrix(NA,nrow=dim(BootCombinations)[1],ncol=dim(BootDis1)[2]))
  
  for (j in 1:dim(BootMembers)[2]){
    BootMembers[,j]<-rowSums(BootDis1[BootCombinations[,1],]==j)+rowSums(BootDis1[BootCombinations[,2],]==j)
  }
  MemberProp<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))
  
  #find members by 
  
  Members<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))
  PropList<-paste(seq(1,max(BootDis1)))
  colnames(Members)<-PropList
  rownames(Members)<-paste(UBound)
  for (j in 1:idxL){
    Members[j,]<-colSums(BootMembers[which((IntDist>=LBound[j]&IntDist<=UBound[j])),])
  }
  MemberProp<-t(t(Members)/colSums(Members))
  return(MemberProp)
}

#Test Membership
TestMembershipProportion<-function(SimDist,BootDis1,BootDis2,ModalityList){
  BootCombinations<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=2))
  BootCombinations[,1]<-rep(seq(1,dim(BootDis1)[1]),each=dim(BootDis2)[1])
  BootCombinations[,2]<-rep(seq(1,dim(BootDis2)[1]),dim(BootDis1)[1])
  #make density
  SimDens<-density(SimDist)
  #find minima
  LocalMinimaIdx<-ModalityList$Minima
  idxL<-length(LocalMinimaIdx)+1
  LBound<-c(0,SimDens$x[LocalMinimaIdx])
  UBound<-c(SimDens$x[LocalMinimaIdx],1)
  #determine the membership of each comparison
  BootMembers<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=dim(BootDis1)[2]+dim(BootDis2)[2]))
  for (j in 1:dim(BootMembers)[2]){
    if (j<=dim(BootDis1)[2]){
      BootMembers[,j]<-rowSums(BootDis1[BootCombinations[,1],]==j)
    } else{
      BootMembers[,j]<-rowSums(BootDis2[BootCombinations[,1],]==(j-dim(BootDis1)[2]))
    }
  }
  MemberProp<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))
  
  
  Members<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))
  PropList<-paste(seq(1,max(BootDis1)))
  colnames(Members)<-PropList
  rownames(Members)<-paste(UBound)
  for (j in 1:idxL){
    Members[j,]<-colSums(BootMembers[which((SimDist>=LBound[j]&SimDist<=UBound[j])),])
  }
  MemberProp<-t(t(Members)/colSums(Members))
  return(MemberProp)
}

#ModalityTest
ModalityTest<-function(SimObject,IntObject1,IntObject2){
  #SimMode<-Modality(SimObject$Summary$Tanimoto)
  #gather modality data
  #NullMode<-Modality(SimObject$NullOut$NullTani)
  Int1Mode<-Modality(IntObject1$InternalTanimoto)
  Int2Mode<-Modality(IntObject2$InternalTanimoto)
  #gather proportion data
  #TMP<-TestMembershipProportion(SimObject$Summary$Tanimoto,SimObject$Boot[[1]],SimObject$Boot[[2]],SimMode)
  #NMP<-NullMembershipProportion(SimObject$NullOut$NullTani,SimObject$Boot[[3]],SimObject$Boot[[4]],NullMode)
  IMP1<-InternalMembershipProportion(IntObject1$InternalTanimoto,SimObject$Boot[[1]],Int1Mode)
  IMP2<-InternalMembershipProportion(IntObject2$InternalTanimoto,SimObject$Boot[[2]],Int2Mode)
  #TZ<-ModalityZ(SimMode,TMP)
  #NZ<-ModalityZ(NullMode,NMP)
  TZ<-NULL
  NZ<-NULL
  IZ1<-ModalityZ(Int1Mode,IMP1)
  IZ2<-ModalityZ(Int2Mode,IMP2)
  RejectList<-list('Test'=TZ,'Null'=NZ,'Int1'=IZ1,'Int2'=IZ2)
  return(RejectList)
}

ModalityInterpreter<-function(RejectList,IntObject1,Mode1,IntObject2,Mode2,SimObject=NULL){
  Int1Out<-NULL
  Int2Out<-NULL
  if (!is.null(RejectList$Int1)){
    PrimaryPeak1<-Mode1[[1]]$Peak_x[which(Mode1[[1]]$N_Members==max(Mode1[[1]]$N_Members))]
    for (j in 1:length(RejectList$Int1)){
      if (!is.null(RejectList$Int1[[j]])){
        for (l in 1:dim(RejectList$Int1[[j]])[2]){
          if ((RejectList$Int1[[j]][1,l]==PrimaryPeak1 & RejectList$Int1[[j]][2,l]<0)|(RejectList$Int1[[j]][1,l]!=PrimaryPeak1 & RejectList$Int1[[j]][2,l]>0)){
            Int1Out<-c(Int1Out,j)
          }
        }
      }
    }
  }
  if (!is.null(RejectList$Int2)){
    PrimaryPeak2<-Mode2[[1]]$Peak_x[which(Mode2[[1]]$N_Members==max(Mode2[[1]]$N_Members))]
    for (j in 1:length(RejectList$Int2)){
      if (!is.null(RejectList$Int2[[j]])){
        for (l in 1:dim(RejectList$Int2[[j]])[2]){
          if ((RejectList$Int2[[j]][1,l]==PrimaryPeak2 & RejectList$Int2[[j]][2,l]<0)|(RejectList$Int2[[j]][1,l]!=PrimaryPeak2 & RejectList$Int2[[j]][2,l]>0)){
            Int2Out<-c(Int2Out,j)
          }
        }
      }
    }
  }
  if (!is.null(Int1Out)){
    Int1Out<-unique(Int1Out)
  }
  if (!is.null(Int2Out)){
    Int2Out<-unique(Int2Out)
  }
  if (verbose){
    print('Due to over-representation in non-primary peaks, dataset 1 should see the removal of replicates:')
    print(Int1Out)
    print('Due to over-representation in non-primary peaks, dataset 2 should see the removal of replicates:')
    print(Int2Out)
  }
  return(list(Int1Out,Int2Out))
  
}

#Modality ZScore
ModalityZ<-function(ModeObject,PropObject,npercent=.1){
  samples<-dim(PropObject)[2]
  PIdx<-which(ModeObject[[1]]$N_Members>=npercent*sum(ModeObject[[1]]$N_Members))
  reject<-list()
  if (dim(PropObject)[1]>1){
    if (length(PIdx)>1){
      for (j in 1:samples){
        xMN<-apply(PropObject[PIdx,-j]/rowSums(PropObject[PIdx,-j]),1,mean)
        xSD<-apply(PropObject[PIdx,-j]/rowSums(PropObject[PIdx,-j]),1,sd)
        xZ<-(PropObject[PIdx,j]/rowSums(PropObject[PIdx,-j])-xMN)/xSD
        if (any(abs(xZ)>=3)){
          zidx<-which(abs(xZ)>=3)
          ztemp<-matrix(nrow=2,ncol=length(zidx))
          ztemp[2,]<-xZ[zidx]
          ztemp[1,]<-ModeObject[[1]][PIdx,"Peak_x"][zidx]
          reject[[j]]<-ztemp
        }
      }
    } else {
      for (j in 1:samples){
        xMN<-mean(PropObject[PIdx,-j]/sum(PropObject[PIdx,-j]))
        xSD<-sd(PropObject[PIdx,-j]/sum(PropObject[PIdx,-j]))
        xZ<-(PropObject[PIdx,j]/sum(PropObject[PIdx,-j])-xMN)/xSD
        if (any(abs(xZ)>=3)){
          zidx<-which(abs(xZ)>=3)
          ztemp<-matrix(nrow=2,ncol=length(zidx))
          ztemp[2,]<-xZ[zidx]
          ztemp[1,]<-ModeObject[[1]][PIdx,'Peak_x'][zidx]
          reject[[j]]<-ztemp
        }
      }
    } 
  }
  if (length(reject)>0){
    return(reject)
  } else {
    return(NULL)
  }
}

#Modality ZScore Subset
ModalityZSubset<-function(ModeObject,PropObject,npercent=.1){
  samples<-dim(PropObject)[2]
  PIdx<-which(ModeObject[[1]]$N_Members>=npercent*sum(ModeObject[[1]]$N_Members))
  reject<-list()
  i<-1
  if ((samples-2)>=2){
    mMax<-samples-2
    mseq<-mMax:2
  }
  
  if (dim(PropObject)[1]>1){
    if (length(PIdx)>1){
      for (l in 1:length(mseq)){
        subsets<-combn(samples,mseq[l])
        slen<-dim(subsets)[2]
        for (j in 1:slen){
          sub<-subsets[,j]
          xMN<-apply(PropObject[PIdx,-sub]/rowSums(PropObject[PIdx,-sub]),1,mean)
          xSD<-apply(PropObject[PIdx,-sub]/rowSums(PropObject[PIdx,-sub]),1,sd)
          xZ<-(PropObject[PIdx,sub]/rowSums(PropObject[PIdx,-sub])-xMN)/xSD
          if (any(abs(xZ)>=3)){
            if (any(rowSums(abs(xZ)>=3)==length(sub))){
              zidx<-which(rowSums(abs(xZ)>=3)==length(sub))
              ztemp<-matrix(nrow=2,ncol=length(zidx))
              ztemp[2,]<-xZ[zidx]
              ztemp[1,]<-ModeObject[[1]][PIdx,"Peak_x"][zidx]
              reject[[i]]<-list()
              reject[[i]][[1]]<-as.data.frame(ztemp)
              reject[[i]][[2]]<-sub
              i<-i+1
            }
            
          }
        }
      }
      
    } else {
      for (l in 1:length(mseq)){
        subsets<-combn(samples,mseq[l])
        slen<-dim(subsets)[2]
        for (j in 1:slen){
          sub<-subsets[,j]
          xMN<-mean(PropObject[PIdx,-sub]/sum(PropObject[PIdx,-sub]),na.rm=T)
          xSD<-sd(PropObject[PIdx,-sub]/sum(PropObject[PIdx,-sub]),na.rm=T)
          xZ<-(PropObject[PIdx,sub]/sum(PropObject[PIdx,-sub])-xMN)/xSD
          if (any(abs(xZ)>=3)){
            if (any(rowSums(abs(xZ)>=3)==length(sub))){
              zidx<-which(rowSums(abs(xZ)>=3)==length(sub))
              ztemp<-matrix(nrow=2,ncol=length(zidx))
              ztemp[2,]<-xZ[zidx]
              ztemp[1,]<-ModeObject[[1]][PIdx,"Peak_x"][zidx]
              reject[[i]]<-list()
              reject[[i]][[1]]<-as.data.frame(ztemp)
              reject[[i]][[2]]<-sub
              i<-i+1
            }
            
          }
        }
      }
    } 
  }
  if (length(reject)>0){
    return(reject)
  } else {
    return(NULL)
  }
}

TheoreticalDataGenerator<-function(n,g,p,w,alp=2,bet=2,maxim=TRUE,w0=2,w1=10){
  #n samples
  #g identifications
  #p average presence
  ## p determines the presence distribution that appears across the sample
  ## binomial distribution rate of observation is proportional to signal intensity
  ## eg P(1,n,p_adj)=#observations. p_adj=p*exp(x-mean(x))
  #w term controlling the variance in the data
  #alp and bet are terms for general data distribution
  datagen<-data.frame(matrix(0,nrow=g,ncol=n))
  rowsampname<-paste0('Identification',1:g)
  row.names(datagen)<-rowsampname
  colnames(datagen)<-paste0('Sample',1:n)
  
  #initialize datavalues
  datavals<-rep(0,g)
  #scalar or vector for initialization
  if (length(alp)==1 & length(bet)==1){
    datavals<-rbeta(g,alp,bet)
  } else if (length(alp)==g & length(bet)==g){
    for (j in 1:g){
      datavals[j]<-rbeta(1,alp[j],bet[j])
    }
  }
  
  #USE MASS:::.rat(mean)$rat to get alphas and betas
  #use rounding to generate betas for each value so they all are within 0-1
  datmean<-mean(datavals)
  for (j in 1:g){
    if (length(w)==n & length(w)!=g){
      for (l in 1:n){
        #make sample variant data
        FractionObj<-MASS:::.rat(round(datavals[j],w0))$rat
        FractionObj<-FractionObj/(w1^(w[l]-6))
        bet0<-FractionObj[2]-FractionObj[1]
        alp0<-FractionObj[1]
        datagen[j,l]<-rbeta(1,alp0,bet0)
      }
    } else if (length(w)==g) {
      #make glycopeptide variant data
      FractionObj<-MASS:::.rat(round(datavals[j],w0))$rat
      FractionObj<-FractionObj/(w1^(w[j]-6))
      bet0<-FractionObj[2]-FractionObj[1]
      alp0<-FractionObj[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    } else {
      #make uniform variant data
      FractionObj<-MASS:::.rat(round(datavals[j],w0))$rat
      FractionObj<-FractionObj/(w1^(w[1]-6))
      bet0<-FractionObj[2]-FractionObj[1]
      alp0<-FractionObj[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    }
    padj<-p*exp(mean(unlist(datagen[j,])-datmean))
    padj[padj>=1]<-0.99999999
    padj[padj<=0]<-0.00000001
    keepnum<-rbinom(1,(n-1),padj)+1
    datagen[j,-match(sample(unlist(datagen[j,]),keepnum),datagen[j,])]<-0
  }
  if (maxim==TRUE){
    datagen[which(rowMeans(datagen)==max(rowMeans(datagen))),]<-rep(1,n)
  }
  FinalOut<-list("Dataset"=datagen,"parameters"=c(n,g,p),"w"=w,"alp"=alp,"bet"=bet,"DistributionCenters"=datavals)
  return(FinalOut)
}



#WithinSim
CombWRep<-function(n,r){
  return(factorial(n+r-1)/(factorial(n-1)*factorial(r)))
}


###Preserved Symmetry Comparison
SymmetricalSimBootstrap<-function(filename1,filename2,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,bootie=TRUE,logopt=FALSE,normvec=c('None','None')){
  #load data and acquire glycopeptides
  file1<-SimDataCleanLogA(filename1,kmin,rel,normvector = normvec[1])
  file2<-SimDataCleanLogA(filename2,kmin,rel,normvector = normvec[2])
  glycopep1<-c(row.names(file1))
  glycopep2<-c(row.names(file2))
  glycojoint<-unique(c(glycopep1,glycopep2))
  mergedf<-merge(file1,file2,by=0,all=TRUE)
  row.names(mergedf)<-mergedf[,1]
  mergedf<-mergedf[,-1]
  coln1<-dim(file1)[2]
  coln2<-dim(file2)[2]
  cols1<-seq(coln1)
  cols2<-seq(coln2)
  
  Nsmpl2<-CombWRep(coln2,coln2)
  Nsmpl1<-CombWRep(coln1,coln1)

  if (((coln1-2)>1) & ((coln1-2)<20)){
    cmseq<-(coln1-2):2
    cm1fact<-sum(factorial(coln1)/(factorial(cmseq)*factorial(coln1-cmseq)))
    MandatoryCombinations1<-matrix(NA,nrow = cm1fact,ncol=coln1)
    jdx<-1
    kdx<-0
    for (j in 1:length(cmseq)){
      kdx<-kdx+(factorial(coln1)/(factorial(cmseq[j])*factorial(coln1-cmseq[j])))
      MandatoryCombinations1[jdx:kdx,1:cmseq[j]]<-t(combn(coln1,cmseq[j]))
      MandatoryCombinations1[jdx:kdx,(cmseq[j]+1):coln1]<-t(apply(t(combn(coln1,cmseq[j])),1,sample,size=coln1-cmseq[j],replace=T))
      jdx<-kdx+1
    }
  }
  if (((coln2-2)>1) & ((coln2-2)<20)){
    cmseq<-(coln2-2):2
    cm2fact<-sum(factorial(coln2)/(factorial(cmseq)*factorial(coln2-cmseq)))
    MandatoryCombinations2<-matrix(NA,nrow = cm2fact,ncol=coln2)
    jdx<-1
    kdx<-0
    for (j in 1:length(cmseq)){
      kdx<-kdx+(factorial(coln2)/(factorial(cmseq[j])*factorial(coln2-cmseq[j])))
      MandatoryCombinations2[jdx:kdx,1:cmseq[j]]<-t(combn(coln2,cmseq[j]))
      MandatoryCombinations2[jdx:kdx,(cmseq[j]+1):coln2]<-t(apply(t(combn(coln2,cmseq[j])),1,sample,size=coln2-cmseq[j],replace=T))
      jdx<-kdx+1
    }
  }
  
  if (bootie){
    if (coln1<6 & 2<coln1){
      Nsmpl1<-CombWRep(coln1,coln1-1)
      combo1<-data.frame(matrix(1,nrow=Nsmpl1,ncol=(coln1)))
      combo1<-data.frame(gtools::combinations(coln1,(coln1),repeats.allowed = TRUE))
    } else if (coln1<=2) {
      print(paste(filename1,' has too few samples. This will not work.'))
      stop()
    } else {
      print('There be a number of samples here. The bootstrap generation will be slower.')
      combo1<-data.frame(matrix(1,nrow=100,ncol=(coln1)))
      for (j in 1:100){
        combo1[j,]<-sort(sample(coln1,(coln1),replace = TRUE))
      }
      if (sum(duplicated(combo1))>0){
        combo1<-combo1[-which(duplicated(combo1)),]
      }
      
    }
    if (coln2<6 & 2<coln2){
      Nsmpl2<-CombWRep(coln2,coln2)
      combo2<-data.frame(matrix(1,nrow=Nsmpl2,ncol=(coln2)))
      combo2<-data.frame(gtools::combinations(coln2,(coln2),repeats.allowed = TRUE))
    } else if (coln2<=2) {
      print(paste(filename2,' has too few samples. This will not work.'))
      stop()
    } else {
      print('There be a number of samples here. The bootstrap generation will be slower.')
      combo2<-data.frame(matrix(1,nrow=100,ncol=(coln2)))
      for (j in 1:100){
        combo2[j,]<-sort(sample(coln2,(coln2),replace = TRUE))
      }
      if (sum(duplicated(combo1))>0){
        combo2<-combo2[-which(duplicated(combo2)),]
      }
      
    }
  } else {
    combo1<-combn(cols1,(coln1))
    combo2<-combn(cols2,(coln2))
  }
  
  
  jac1<-rep(0,nrow(combo1)*nrow(combo2)) # jaccard holder depricated
  tan1<-rep(0,nrow(combo1)*nrow(combo2)) # tanimoto holder
  tan1Final<-rep(0,nrow(combo1)*nrow(combo2))
  tanmathold<-data.frame(matrix(data=0,nrow=(nrow(combo1)*nrow(combo2)) ,ncol=length(glycojoint)))
  colnames(tanmathold)<-glycojoint
  tanmatholdW<-t(data.frame(matrix(data=0,nrow=(nrow(combo1)*nrow(combo2)) ,ncol=length(glycojoint))))
  row.names(tanmatholdW)<-glycojoint
  #acquire glycopeptides
  gl1<-glycopep1
  gl2<-glycopep2
  
  #iterate through all combinations used in within of file 1
  T1_<-data.frame(matrix(0,nrow=length(gl1),ncol=nrow(combo1)))
  row.names(T1_)<-gl1
  T10<-data.frame(matrix(0,nrow=1,ncol=nrow(combo1)))
  PHold1<-data.frame(matrix(0,nrow=length(gl1),ncol=nrow(combo1)))
  row.names(PHold1)<-gl1
  for (j in 1:nrow(combo1)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo1[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
    }
    T1_[,j]<-rowMeans(temp1,na.rm =TRUE)
    T10[,j]<-sum((rowMeans(temp1,na.rm =TRUE))^2)
    for (m in 1:length(gl1)){
      nacheck<-1-sum(temp1[gl1[m],]==0)/ncol(temp1)
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHold1[gl1[m],j]<-nacheck
    }
  }
  #iterate through all combinations used in within of file 2
  T_1<-data.frame(matrix(0,nrow=length(gl2),ncol=nrow(combo2)))
  row.names(T_1)<-gl2
  T01<-data.frame(matrix(0,nrow=1,ncol=nrow(combo2)))
  PHold2<-data.frame(matrix(0,nrow=length(gl2),ncol=nrow(combo2)))
  row.names(PHold2)<-gl2
  for (j in 1:nrow(combo2)){
    #separate datasets for combinations
    temp2<-data.frame(file2[,unlist(combo2[j,])]) # subset of first data
    row.names(temp2)<-glycopep2 # GP names
    if (MVCorrection!=TRUE){
      temp2[temp2==0]<-NA
    }
    T_1[,j]<-rowMeans(temp2,na.rm =TRUE)
    T01[,j]<-sum((rowMeans(temp2,na.rm =TRUE))^2)
    for (m in 1:length(gl2)){
      nacheck<-1-sum(temp2[gl2[m],]==0)/ncol(temp2)
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHold2[gl2[m],j]<-nacheck
    }
  }
  
  T__Hold<-merge(T1_,T_1,by=0,all=TRUE)
  row.names(T__Hold)<-T__Hold[,1]
  TRef<-row.names(T__Hold)
  T__Hold<-T__Hold[,-1]
  T__Hold[is.na(T__Hold)]<-0
  PHold<-merge(PHold1,PHold2,by=0,all=TRUE)
  row.names(PHold)<-PHold[,1]
  PHold<-PHold[,-1]
  PHold[is.na(PHold)]<-0
  #Bring T11 related terms together
  ncomb1<-nrow(combo1)
  ncomb2<-nrow(combo2)
  ncombt<-ncomb1+ncomb2
  ncombs<-ncomb1*ncomb2
  for (j in 1:nrow(combo1)){
    dT<-abs(T__Hold[,j]-T__Hold[,(ncomb1+1):ncombt])
    row.names(dT)<-row.names(T__Hold)
    presence<-(PHold[,j]+PHold[,(ncomb1+1):ncombt])/2
    if (mn==FALSE){
      KTerm<-1+presence
      row.names(KTerm)<-row.names(PHold)
    } else {
      KTerm<-mn
    }
    for (m in 1:length(TRef)){
      tanmathold[(((j-1)*ncomb2)+1):(j*ncomb2),TRef[m]]<-unlist(T__Hold[TRef[m],j]*T__Hold[TRef[m],(ncomb1+1):ncombt]*(KTerm[TRef[m],]^(-dT[TRef[m],])))
    }
    T11<-t(tanmathold[(((j-1)*ncomb2)+1):(j*ncomb2),])
    for (m in 1:length(TRef)){
      tanmatholdW[TRef[m],(((j-1)*ncomb2)+1):(j*ncomb2)]<-unlist(T11[TRef[m],]/(unlist(T10[j])+T01-colSums(T11,na.rm=T)))
    }
  }
  
  tan1<-colSums(tanmatholdW,na.rm=T)
  #calculate null
  colnt<-coln1+coln2
  if (bootie){
    if (colnt<5 & 2<colnt){
      Nsmplt1<-CombWRep(colnt,coln1)
      Nsmplt2<-CombWRep(colnt,coln2)
      combot1<-data.frame(matrix(1,nrow=Nsmplt1,ncol=coln1))
      combot1<-data.frame(gtools::combinations(colnt,coln1,repeats.allowed = TRUE))
      combot2<-data.frame(matrix(1,nrow=Nsmplt2,ncol=coln2))
      combot2<-data.frame(gtools::combinations(colnt,coln2,repeats.allowed = TRUE))
    } else {
      print('There be a number of samples here. The bootstrap generation will be slower.')
      combot1<-data.frame(matrix(1,nrow=150,ncol=coln1))
      combot2<-data.frame(matrix(1,nrow=150,ncol=coln2))
      for (j in 1:150){
        combot1[j,]<-sort(sample(colnt,coln1,replace = TRUE))
        combot2[j,]<-sort(sample(colnt,coln2,replace = TRUE))
      }
      if (sum(duplicated(combot1))>0){
        combot1<-combot1[-which(duplicated(combot1)),]
      }
      if (sum(duplicated(combot2))>0){
        combot2<-combot2[-which(duplicated(combot2)),]
      }
      
    }
     
  } else {
    combot<-'Error'
  }
  row.names(combot2)<-NULL
  row.names(combot1)<-NULL
  
  #remove all with less than 2 of each side
  if (coln1>3 & coln2>3){
    minsample<-1
  } else {
    minsample<-0
  }
  if (sum(rowSums(combot1<=coln1)>=(coln1-minsample))>0 | sum(rowSums(combot1>coln1)>=(coln1-minsample))>0){
    if (sum(rowSums(combot1<=coln1)>=(coln1-1))>0){
      combot1<-combot1[-which(rowSums(combot1<=coln1)>=(coln1-minsample)),]
      row.names(combot1)<-NULL
    } 
    if (sum(rowSums(combot1>coln1)>=(coln1-minsample))>0){
      combot1<-combot1[-which(rowSums(combot1>coln1)>=(coln1-minsample)),]
      row.names(combot1)<-NULL
    }
  }
  
  if (sum(rowSums(combot2>coln1)>=(coln2-minsample))>0 | sum(rowSums(combot2<=coln1)>=(coln2-minsample))>0){
    if (sum(rowSums(combot2<=coln1)>=(coln2-minsample))>0){
      combot2<-combot2[-which(rowSums(combot2<=coln1)>=(coln2-minsample)),]
      row.names(combot2)<-NULL
    } 
    if (sum(rowSums(combot2>coln1)>=(coln2-minsample))>0){
      combot2<-combot2[-which(rowSums(combot2>coln1)>=(coln2-minsample)),]
      row.names(combot2)<-NULL
    }
  }
  row.names(combot2)<-NULL
  row.names(combot1)<-NULL
  #even out null distributions
  if (coln1==coln2){
    combot<-rbind(combot1,combot2,all=T,fill=T)
    if (sum(duplicated(combot))>0){
      combot<-combot[-which(duplicated(combot)),]
    }
    combot1<-combot[1:floor(dim(combot)[1]/2),]
    combot2<-combot[(floor(dim(combot)[1]/2)+1):(dim(combot)[1]-1),]
    row.names(combot2)<-NULL
    row.names(combot1)<-NULL
  }
  row.names(combot2)<-NULL
  row.names(combot1)<-NULL
  
  
  jacN<-rep(0,nrow(combot1)*nrow(combot2)) # jaccard holder depricated
  tanN<-rep(0,nrow(combot1)*nrow(combot2)) # tanimoto holder
  tanNFinal<-rep(0,nrow(combot1)*nrow(combot2))
  tanmatholdN<-data.frame(matrix(data=0,nrow=(nrow(combot1)*nrow(combot2)) ,ncol=length(glycojoint)))
  colnames(tanmatholdN)<-glycojoint
  tanmatholdNW<-t(data.frame(matrix(data=0,nrow=(nrow(combot1)*nrow(combot2)) ,ncol=length(glycojoint))))
  row.names(tanmatholdNW)<-glycojoint
  #acquire glycopeptides
  gl1<-glycopep1
  gl2<-glycopep2
  gj<-glycojoint
  
  #first set of null distributions
  T1_N<-data.frame(matrix(0,nrow=length(gj),ncol=nrow(combot1)))
  row.names(T1_N)<-gj
  T10N<-data.frame(matrix(0,nrow=1,ncol=nrow(combot1)))
  PHoldN1<-data.frame(matrix(0,nrow=length(gj),ncol=nrow(combot1)))
  row.names(PHoldN1)<-gj
  for (j in 1:nrow(combot1)){
    #separate datasets for combinations
    temp1<-data.frame(mergedf[,unlist(combot1[j,])]) # subset of first data
    row.names(temp1)<-gj # GP names
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
    }
    T1_N[,j]<-rowMeans(temp1,na.rm =TRUE)
    T10N[,j]<-sum((rowMeans(temp1,na.rm =TRUE))^2)
    for (m in 1:length(gj)){
      nacheck<-1-sum(temp1[gj[m],]==0)/ncol(temp1)
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHoldN1[gj[m],j]<-nacheck
    }
  }
  #second set of null distributions
  T_1N<-data.frame(matrix(0,nrow=length(gj),ncol=nrow(combot2)))
  row.names(T_1N)<-gj
  T01N<-data.frame(matrix(0,nrow=1,ncol=nrow(combot2)))
  PHoldN2<-data.frame(matrix(0,nrow=length(gj),ncol=nrow(combot2)))
  row.names(PHoldN2)<-gj
  for (j in 1:nrow(combot2)){
    #separate datasets for combinations
    temp2<-data.frame(mergedf[,unlist(combot2[j,])]) # subset of first data
    row.names(temp2)<-gj # GP names
    if (MVCorrection!=TRUE){
      temp2[temp2==0]<-NA
    }
    T_1N[,j]<-rowMeans(temp2,na.rm =TRUE)
    T01N[,j]<-sum((rowMeans(temp2,na.rm =TRUE))^2)
    for (m in 1:length(gj)){
      nacheck<-1-sum(temp2[gj[m],]==0)/ncol(temp2)
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHoldN2[gj[m],j]<-nacheck
    }
  }
  
  T__HoldN<-merge(T1_N,T_1N,by=0,all=TRUE)
  row.names(T__HoldN)<-T__HoldN[,1]
  T__HoldN<-T__HoldN[,-1]
  T__HoldN[is.na(T__HoldN)]<-0
  PHoldN<-merge(PHoldN1,PHoldN2,by=0,all=TRUE)
  row.names(PHoldN)<-PHoldN[,1]
  PHoldN<-PHoldN[,-1]
  PHoldN[is.na(PHoldN)]<-0
  #Bring T11N related terms together
  ncombt1<-nrow(combot1)
  ncombt2<-nrow(combot2)
  ncombtt<-ncombt1+ncombt2
  ncombts<-ncombt1*ncombt2
  for (j in 1:nrow(combot1)){
    dTN<-abs(T__HoldN[,j]-T__HoldN[,(ncombt1+1):ncombtt])
    row.names(dTN)<-row.names(T__HoldN)
    presenceN<-(PHoldN[,j]+PHoldN[,(ncombt1+1):ncombtt])/2
    if (mn==FALSE){
      KTermN<-1+presenceN
      row.names(KTermN)<-row.names(PHoldN)
    } else {
      KTermN<-mn
    }
    for (m in 1:length(TRef)){
      tanmatholdN[(((j-1)*ncombt2)+1):(j*ncombt2),TRef[m]]<-unlist(T__HoldN[TRef[m],j]*T__HoldN[TRef[m],(ncombt1+1):ncombtt]*(KTermN[TRef[m],]^(-dTN[TRef[m],])))
    }
    T11N<-t(tanmatholdN[(((j-1)*ncombt2)+1):(j*ncombt2),])
    for (m in 1:length(TRef)){
      tanmatholdNW[TRef[m],(((j-1)*ncombt2)+1):(j*ncombt2)]<-unlist(T11N[TRef[m],]/(unlist(T10N[j])+T01N-colSums(T11N,na.rm=T)))
    }
  }
  tanN<-colSums(tanmatholdNW)
  
  #calculate actual similarity
  T1_a<-rowMeans(file1,na.rm =TRUE)
  T10a<-sum((rowMeans(file1,na.rm =TRUE))^2)
  presence1<-data.frame(matrix(0,nrow=length(gl1),ncol=1))
  row.names(presence1)<-gl1
  for (m in 1:length(gl1)){
    nacheck<-1-sum(file1[gl1[m],]==0)/ncol(file1)
    if (is.na(nacheck)){
      nacheck<-0
    }
    presence1[gl1[m],1]<-nacheck
  }
  T_1a<-rowMeans(file2,na.rm =TRUE)
  T01a<-sum((rowMeans(file2,na.rm =TRUE))^2)
  presence2<-data.frame(matrix(0,nrow=length(gl2),ncol=1))
  row.names(presence2)<-gl2
  for (m in 1:length(gl2)){
    nacheck<-1-sum(file2[gl2[m],]==0)/ncol(file2)
    if (is.na(nacheck)){
      nacheck<-0
    }
    presence2[gl2[m],1]<-nacheck
  }
  presenceA<-merge(presence1,presence2,by=0,all=TRUE)
  row.names(presenceA)<-presenceA[,1]
  presenceA<-presenceA[,-1]
  presenceA[is.na(presenceA)]<-0
  THA<-merge(T1_a,T_1a,by=0,all=TRUE)
  row.names(THA)<-THA[,1]
  THA<-THA[,-1]
  THA[is.na(THA)]<-0
  if (mn==FALSE){
    KTerm<-data.frame(matrix((1+rowMeans(presenceA,na.rm=T)),nrow=length(TRef),ncol=1))
    row.names(KTerm)<-row.names(presenceA)
  } else {
    KTerm<-rep(mn,nrow(THA))
    row.names(KTerm)<-row.names(presenceA)
  }
  dTA<-data.frame(matrix(abs(THA[,1]-THA[,2]),nrow=length(TRef),ncol=1))
  row.names(dTA)<-TRef
  tanmatA<-data.frame(matrix(0,ncol=length(TRef),nrow=1))
  tanmatAFinal<-data.frame(matrix(0,ncol=length(TRef),nrow=1))
  colnames(tanmatA)<-TRef
  colnames(tanmatAFinal)<-TRef
  T11A<-data.frame(matrix(0,nrow=length(TRef),ncol=1))
  row.names(T11A)<-TRef
  for (m in 1:length(TRef)){
    T11A[TRef[m],1]<-THA[TRef[m],1]*THA[TRef[m],2]*(KTerm[TRef[m],1]^(-dTA[TRef[m],1]))
    tanmatA[1,TRef[m]]<-unlist(T11A[TRef[m],])
  }
  for (m in 1:length(TRef)){
    tanmatAFinal[1,TRef[m]]<-unlist(T11A[TRef[m],]/(T10a+T01a-sum(T11A,na.rm=T)))
  }
  taniActualF<-sum(tanmatAFinal)
  
  
  
  Output<-data.frame(matrix(c(jac1,tan1),nrow=ncomb1*ncomb2,ncol=2))
  colnames(Output)<-c('Jaccard','Tanimoto')
  TaniOut<-data.frame(matrix(tan1,nrow=ncomb1*ncomb2[1],ncol=1))
  colnames(TaniOut)<-c('Tanimoto')
  NullOut<-data.frame(matrix(tanN,nrow=ncombts,ncol = 1))
  colnames(NullOut)<-c('NullTani')
  FinalOut<-list("Summary"=TaniOut,"Hold"=Output,"RankInfo"=tanmathold,"NullRankInfoFinal"=tanmatholdNW,"NullOut"=NullOut,'RankInfoFinal'=tanmatholdW,"RankInfoActual"=tanmatAFinal,"Actual"=taniActualF,"Boot"=list(combo1,combo2,combot1,combot2))
  return(FinalOut)
}

#Use this to convery vector of similarities into density distribution
SimDisToProbDis<-function(similarityd,breaks=100){
  incr<-1/breaks
  starter<-incr/2
  hold<-rep(0,breaks)
  mids<-seq(starter,1,by=incr)
  redge<-seq(incr,1,by=incr)
  hold[1]<-sum(similarityd<=redge[1])
  for (j in 2:breaks){
    hold[j]<-sum(similarityd<=redge[j] & similarityd>redge[(j-1)])
  }
  hold<-hold+1/breaks*length(similarityd)
  dens<-hold/sum(hold)
  return(dens)
}

EntropicInformation<-function(similarityd,breaks=20){
  dens<-SimDisToProbDis(similarityd,breaks)
  Entropy<--sum(dens*log(dens))
  return(Entropy)
}

SimKLD<-function(TestSim,RefSim,breaks=20){
  TestEntropy<-EntropicInformation(TestSim,breaks)
  TestDens<-SimDisToProbDis(TestSim,breaks)
  RefDens<-SimDisToProbDis(RefSim,breaks)
  KLD<-sum(TestDens*log(TestDens/RefDens))
  return(c(KLD,TestEntropy))
}
  
InternalSimilarity<-function(filename,BootSet,kmin=1,rel=TRUE,MVCorrection=TRUE,mn=FALSE,logopt=FALSE,normvec=c('None')){
  file1<-SimDataCleanLogA(filename,kmin,rel,normvector = normvec)
  glycopep1<-c(row.names(file1))
  glycojoint<-glycopep1
  glycopep2<-glycopep1
  coln1<-dim(file1)[2]
  cols1<-seq(coln1)
  Nsmpl1<-CombWRep(coln1,coln1-1)
  combo1<-BootSet
  ncomps<-nrow(combo1)*(nrow(combo1)-1)/2
  comps<-combn(nrow(combo1),2)
  
  jac1<-rep(0,ncomps) # jaccard holder depricated
  tan1<-rep(0,ncomps) # tanimoto holder
  tan1Final<-rep(0,ncomps)
  tanmathold<-data.frame(matrix(data=0,nrow=ncomps ,ncol=length(glycojoint)))
  colnames(tanmathold)<-glycojoint
  tanmatholdW<-t(data.frame(matrix(data=0,nrow=ncomps ,ncol=length(glycojoint))))
  row.names(tanmatholdW)<-glycojoint
  #acquire glycopeptides
  gl1<-glycopep1
  gl2<-glycopep2
  
  #iterate through all combinations used in within of file 1
  T1_<-data.frame(matrix(0,nrow=length(gl1),ncol=nrow(combo1)))
  row.names(T1_)<-gl1
  T10<-data.frame(matrix(0,nrow=1,ncol=nrow(combo1)))
  PHold1<-data.frame(matrix(0,nrow=length(gl1),ncol=nrow(combo1)))
  row.names(PHold1)<-gl1
  for (j in 1:nrow(combo1)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo1[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
    }
    T1_[,j]<-rowMeans(temp1,na.rm =TRUE)
    T10[,j]<-sum(rowMeans(temp1,na.rm =TRUE)^2)
    for (m in 1:length(gl1)){
      nacheck<-1-sum(temp1[gl1[m],]==0)/ncol(temp1)
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHold1[gl1[m],j]<-nacheck
    }
  }
  T_1<-T1_
  T01<-T10
  PHold2<-PHold1
  
  T__Hold<-merge(T1_,T_1,by=0,all=TRUE)
  row.names(T__Hold)<-T__Hold[,1]
  TRef<-row.names(T__Hold)
  T__Hold<-T__Hold[,-1]
  T__Hold[is.na(T__Hold)]<-0
  PHold<-merge(PHold1,PHold2,by=0,all=TRUE)
  row.names(PHold)<-PHold[,1]
  PHold<-PHold[,-1]
  PHold[is.na(PHold)]<-0
  #Bring T11 related terms together
  ncomb1<-nrow(combo1)
  ncombt<-2*ncomb1
  ncombs<-ncomb1*(ncomb1-1)/2
  #prevent duplicate comparisons
  tracker1<-0
  ncomb2<-ncomb1-1
  for (j in 1:ncomb2){
    dT<-abs(T__Hold[,j]-T__Hold[,(j+1):ncomb1])
    if (is.null(dim(dT))){
      dT<-data.frame(matrix(dT,nrow=length(dT),ncol=1))
      row.names(dT)<-row.names(T__Hold)
    } else {
      row.names(dT)<-row.names(T__Hold)
    }
    
    presence<-(PHold[,j]+PHold[,(j+1):ncomb1])/2
    if (mn==FALSE){
      KTerm<-1+presence
      if (is.null(dim(KTerm))){
        KTerm<-data.frame(matrix(KTerm,nrow=length(KTerm),ncol=1))
        row.names(KTerm)<-row.names(PHold)
      } else {
        row.names(KTerm)<-row.names(PHold)
      }
    } else {
      KTerm<-mn
    }
    for (m in 1:length(TRef)){
      tanmathold[(tracker1+1):(tracker1+ncomb1-j),TRef[m]]<-unlist(T__Hold[TRef[m],j]*T__Hold[TRef[m],(ncomb1+1+j):ncombt]*(KTerm[TRef[m],]^(-dT[TRef[m],])))
    }
    T11<-t(tanmathold[(tracker1+1):(tracker1+ncomb1-j),])
    for (m in 1:length(TRef)){
      tanmatholdW[TRef[m],(tracker1+1):(tracker1+ncomb1-j)]<-unlist(T11[TRef[m],]/(unlist(T10[j])+T01[-c(1:j)]-colSums(T11,na.rm=T)))
    }
    #remove possible combinations increment index
    tracker1<-tracker1+ncomb1-j
  }
  tan1<-colSums(tanmatholdW)
  
  
  OutputObj<-list("InternalTanimoto"=tan1,"InternalRankingInfo"=tanmatholdW)
  return(OutputObj)
}

#plots internal boot object[[]]: 1 is the first
InternalQuality<-function(filename,BootSet,SimilarityObj,PlotTitle,GroupName,Int=NULL,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,verbose=FALSE,legendpos='topleft',legendval=FALSE,xbound=c(0,1),logval=F){
  if (!is.null(Int)){
    #Use Int as Int object
  } else {
    Int<-InternalSimilarity(filename,BootSet,kmin,rel,MVCorrection,mn)
  }
  
  k<-100
  #build density
  TDis<-SimilarityObj$Summary$Tanimoto
  TDis[is.infinite(TDis)]<-NA
  dT<-density(TDis,from=-0.1,to=1.1,na.rm=T)
  sdTy<-k*dT$y/sum(dT$y)
  if (logval==T){
    sdTy<-log(sdTy+1)
    if (any(is.infinite(sdTy))){
      sdTy[which(is.infinite(sdTy))]<-0
    }
  }
  NDis<-SimilarityObj$NullOut$NullTani
  NDis[is.infinite(NDis)]<-NA
  dN<-density(NDis,from=-0.1,to=1.1,na.rm=T)
  sdNy<-k*dN$y/sum(dN$y)
  if (logval==T){
    sdNy<-log(sdNy+1)
    if (any(is.infinite(sdNy))){
      sdNy[which(is.infinite(sdNy))]<-0
    }
  }
  IDis<-Int$InternalTanimoto
  IDis[is.infinite(IDis)]<-NA
  dI<-density(IDis,from=-0.1,to=1.1,na.rm=T)
  sdIy<-k*dI$y/sum(dI$y)
  if (logval==T){
    sdIy<-log(sdIy+1)
    if (any(is.infinite(sdIy))){
      sdIy[which(is.infinite(sdIy))]<-0
    }
  }
  TAct<-SimilarityObj$Actual
  OverlapData<-OverlapCalculator(IDis,TDis)
  #plot densities
  if (max(TDis)==0){
    mh<-max(c(k*dI$y/sum(dI$y)))
    plot(c(0,0,0.1,0.1),c(0,mh,mh,0),xlim=c(0,1),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution')
    polygon(c(0,0,0.1,0.1),c(0,mh,mh,0),col=rgb(1,0,0,0.5))
  } else {
    mh<-max(c(k*dT$y/sum(dT$y),k*dI$y/sum(dI$y)))
    plot(dT$x,k*dT$y/sum(dT$y),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution',xlim=xbound)
    polygon(c(-0.1001,dT$x,1.1001),c(0,k*dT$y/sum(dT$y),0),col=rgb(1,0,0,0.5))
  }
  #lines(dN$x,k*dN$y/sum(dN$y))
  #polygon(dN$x,k*dN$y/sum(dN$y),col=rgb(0,0,1,0.5))
  lines(dI$x,k*dI$y/sum(dI$y))
  polygon(c(-0.1001,dI$x,1.1001),c(0,k*dI$y/sum(dI$y),0),col=rgb(0,1,0,0.5))
  lines(rep(TAct,2),c(0,100),col=1,lwd=3)
  #percentile location
  CompPerc<-ecdf(TDis)(TAct)
  JointPerc<-ecdf(NDis)(TAct)
  
  Overlap<-OverlapData$PercOverlap
  AlphaValue<-round(OverlapData$FP,3)
  BetaValue<-round(OverlapData$FN,3)
  
  if (AlphaValue>(10^-10)){
    #make alpha area
    for (j in 1:(length(OverlapData$FPi)/2)){
      jdxh<-j*2
      jdxl<-(j-1)*2+1
      fpi<-seq(OverlapData$FPi[jdxl],OverlapData$FPi[jdxh])
      xi<-mean(c(dT$x[fpi[1]],dT$x[fpi[1]-1]),na.rm=T)
      x0<-mean(c(dT$x[fpi[length(fpi)]],dT$x[fpi[length(fpi)]+1]),na.rm=T)
      polygon(c(xi,dT$x[fpi],x0,x0),c(0,sdTy[fpi],sdTy[fpi[length(fpi)]],0),col='black')
    }
  }
  if (BetaValue>(10^-10)){
    #make beta area
    for (j in 1:(length(OverlapData$FNi)/2)){
      jdxh<-j*2
      jdxl<-(j-1)*2+1
      fni<-seq(OverlapData$FNi[jdxl],OverlapData$FNi[jdxh])
      xi<-mean(c(dT$x[fni[1]],dT$x[fni[1]-1]),na.rm=T)
      x0<-mean(c(dT$x[fni[length(fni)]],dT$x[fni[length(fni)]+1]),na.rm=T)
      polygon(c(xi,dN$x[fni],x0,x0),c(0,sdIy[fni],sdIy[fni[length(fni)]],0),col='darkgrey')
    }
  }
  
  deltaIT<-1+(mean(IDis,na.rm=T)-mean(TDis,na.rm = T))
  Confidence<-round(deltaIT*log10((k*max(dI$y)/sum(dI$y))/var(IDis,na.rm=T)*10^(-AlphaValue-BetaValue)),2)
  if (legendval){
    legend(legendpos,legend=c(paste('Observed Similarity in the',round(ecdf(TDis)(TAct),2)*100,'Percentile'),paste('Internal of',GroupName,'with Score=',Confidence),'Test Distribution',paste0(AlphaValue*100,'% False Positive Rate'),paste0(BetaValue*100,'% False Negative Rate')),fill=c(NA,3,2,'darkgrey','black'),lty=c(1,rep(NA,4)),lwd=c(3,NA,NA,NA,NA),density=c(0,NA,NA,NA,NA),border=c(NA,1,1,1,1))
  }
  if (verbose==T){
    return(list(Int,Overlap,XPoint,CompPerc,JointPerc,Confidence,AlphaValue,BetaValue))
  }

  
}


#ConfidenceScoreZFormat

ConfidenceScore<-function(IntDist,SimDist){
  IntDens<-density(IntDist,na.rm=T)
  SimDens<-density(SimDist,na.rm=T)
  sdI<-sd(IntDist,na.rm=T)
  ET<-mean(SimDist,na.rm=T)
  EI<-mean(IntDist,na.rm=T)
  overlaplist<-Overlap(SimDist,IntDist)
  dSim<-abs(EI-ET)
  c<-unname(unlist(overlaplist['PercOverlap']))/100
  Score<-ScoreUnder(dSim,sdI,c)
  return(Score)
}

ScoreUnder<-function(dSim,sdSim,c){
  Score<-(dSim/sdI)*(1-c)
  return(Score)
}
#Overlap Function
Overlap<-function(Dis1,Dis2){
  density1<-density(Dis1,from=-0.1,to=1.1,na.rm=T)
  density2<-density(Dis2,from=-0.1,to=1.1,na.rm=T)
  df <- merge(
    as.data.frame(density1[c("x", "y")]),
    as.data.frame(density2[c("x", "y")]),
    by = "x", suffixes = c(".A", ".B")
  )
  df$comp <- as.numeric(df$y.A > df$y.B)
  df$cross <- c(NA, diff(df$comp))
  tempidx<-which(df$cross!=0)
  CPoint<-tempidx[which(max(df$y.A[tempidx])==df$y.A[tempidx])]
  XPoint<-df$x[CPoint]
  if (length(CPoint)==0){
    Overlap<-0
    Alpha<-0
    Beta<-0
  } else if (df$y.A[CPoint]<(10^-10)){
    Overlap<-0
    Alpha<-0
    Beta<-0
  } else {
    Area1<-(1-ecdf(Dis1)(XPoint))
    Area2<-ecdf(Dis2)(XPoint)
    Overlap1<-round(Area1+Area2,2)
    Overlap<-round(100*Overlap1/(2-Overlap1),2)
    Alpha<-round(Area1/(2-Overlap1),2)
    Beta<-round(Area2/(2-Overlap1),2)
  }
  return(list('PercOverlap'=Overlap,'Alpha'=Alpha,'Beta'=Beta))
}

#plots observed, test, null

SimPlot<-function(PlotTitle,SimilarityObj,legendpos='topleft',verbose=F,legendval=FALSE,xbound=c(0,1),logval=FALSE){
  k<-100
  #build density
  TDis<-SimilarityObj$Summary$Tanimoto
  dT<-density(TDis,from=-0.1,to=1.1,na.rm=T)
  sdTy<-k*dT$y/sum(dT$y)
  if (logval==T){
    sdTy<-log(sdTy+1)
    if (any(is.infinite(sdTy))){
      sdTy[which(is.infinite(sdTy))]<-0
    }
  }
  
  NDis<-SimilarityObj$NullOut$NullTani
  dN<-density(NDis,from=-0.1,to=1.1,na.rm=T)
  sdNy<-k*dN$y/sum(dN$y)
  if (logval==T){
    sdNy<-log(sdNy+1)
    if (any(is.infinite(sdNy))){
      sdNy[which(is.infinite(sdNy))]<-0
    }
  }
  TAct<-SimilarityObj$Actual
  
  OverlapData<-OverlapCalculator(SimilarityObj$NullOut$NullTani,SimilarityObj$Summary$Tanimoto)
  if (logval==T){
    #OverlapData<-OverlapCalculator(SimilarityObj$NullOut$NullTani,SimilarityObj$Summary$Tanimoto)
  }
  #plot densities
  if (max(TDis)==0){
    mh<-max(sdNy)
    plot(c(0,0,0.05,0.05),c(0,mh,mh,0),xlim=c(0,1),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='Density')
    polygon(c(0,0,0.05,0.05),c(0,mh,mh,0),col=rgb(1,0,0,0.5))
  } else {
    mh<-max(c(sdNy,sdTy))
    plot(dT$x,sdTy,xlim=xbound,ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='Density')
    polygon(c(-0.1001,dT$x,1.1001),c(0,sdTy,0),col=rgb(1,0,0,0.5))
  }
  lines(dN$x,sdNy)
  polygon(c(-0.1001,dN$x,1.1001),c(0,sdNy,0),col=rgb(0,0,1,0.5))
  lines(rep(TAct,2),c(0,100),col=1,lwd=3)
  
  
  #percentile location
  Overlap<-OverlapData$PercOverlap
  AlphaValue<-OverlapData$FP
  BetaValue<-OverlapData$FN
  if (legendval){
    legend(legendpos,legend=c(paste('Observed Similarity in the',round(ecdf(TDis)(TAct),2)*100,'Percentile'),'Test Distribution','Null Distribution',paste0(round(AlphaValue,3)*100,'% False Positive Rate'),paste0(round(BetaValue,3)*100,'% False Negative Rate')),fill=c(NA,rgb(1,0,0,0.5),rgb(0,0,1,0.5),'darkgray','black'),lty=c(1,rep(NA,4)),density=c(0,NA,NA,NA,NA),border=c(NA,1,1,1,1))
  }
  
  if (AlphaValue>(10^-10)){
    #make alpha area
    for (j in 1:(length(OverlapData$FPi)/2)){
      jdxh<-j*2
      jdxl<-(j-1)*2+1
      fpi<-seq(OverlapData$FPi[jdxl],OverlapData$FPi[jdxh])
      xi<-mean(c(dT$x[fpi[1]],dT$x[fpi[1]-1]),na.rm=T)
      x0<-mean(c(dT$x[fpi[length(fpi)]],dT$x[fpi[length(fpi)]+1]),na.rm=T)
      polygon(c(xi,dT$x[fpi],x0,x0),c(0,sdTy[fpi],sdTy[fpi[length(fpi)]],0),col='black')
    }
  }
  if (BetaValue>(10^-10)){
    #make beta area
    for (j in 1:(length(OverlapData$FNi)/2)){
      jdxh<-j*2
      jdxl<-(j-1)*2+1
      fni<-seq(OverlapData$FNi[jdxl],OverlapData$FNi[jdxh])
      xi<-mean(c(dT$x[fni[1]],dT$x[fni[1]-1]),na.rm=T)
      x0<-mean(c(dT$x[fni[length(fni)]],dT$x[fni[length(fni)]+1]),na.rm=T)
      polygon(c(xi,dN$x[fni],x0,x0),c(0,sdNy[fni],sdNy[fni[length(fni)]],0),col='darkgrey')
    }
  }
  
  if (verbose==T){
    return(list(Overlap,TAct,'FalsePositiveRate'=AlphaValue,'FalseNegativeRate'=BetaValue))
  }
  
}

SimPlotFromFile<-function(PlotTitle,filename1,filename2,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,bootie=TRUE){
  SimObj<-SymmetricalSimBootstrap(filename1,filename2,kmin,rel,MVCorrection,mn,bootie)
  SimPlot(PlotTitle,SimObj)
  return(SimObj)
}




SimDensFunction<-function(datavec,xlim=c(0,1),smooth=T,bins=100){
  regions<-seq(xlim[1],xlim[2],by=((xlim[2]-xlim[1])/bins))
  counts<-rep(0,bins-1)
  for (j in 1:(bins-1)){
    if (j==(bins-1)){
      counts[j]<-sum(datavec<=regions[j+1] & datavec>=regions[j])
    } else {
      counts[j]<-sum(datavec<regions[j+1] & datavec>=regions[j])
    }
    
  }
  return(counts/sum(counts))      
}



RankSimilarity<-function(RankInfo,PlotTitle){
  rInfo<-apply(RankInfo,1,order)
  cInfo<-colMeans(RankInfo)
  boxplot(RankInfo[,order(cInfo)],las=2,xlab='Ranking',ylim=c(0,1),ylab='Similarity Contribution',main=PlotTitle)
  return(list('RankOrder'=order(cInfo),'RankInfo'=RankInfo,'RankNames'=colnames(RankInfo)[order(cInfo)]))
}

OverlapCalculator<-function(ReferenceDis,TestDis){
  #Let the distributions be the needed input
  dR<-density(ReferenceDis,from=-0.1,to=1.1)
  dT<-density(TestDis,from=-.1,to=1.1)
  if (max(TestDis)==0){
    dT$y<-rep(0,length(dT$y))
    dT$y[which(dT$x>=-0.001 & dT$x<=0.001)]<-1
  }
  sdT<-dT$y/sum(dT$y)
  sdR<-dR$y/sum(dR$y)
  #FP == Ref greater than Test
  FalsePosIDX<-which((sdT<sdR) & (sdT>0.000009))
  FP_T_ecdf<-ecdf(TestDis)
  #FN == Test greater than Ref
  FalseNegIDX<-which((sdT>=sdR) & (sdR>0.000009))
  FN_R_ecdf<-ecdf(ReferenceDis)
  
  FPcont<-which(diff(FalsePosIDX)!=1)
  FNcont<-which(diff(FalseNegIDX)!=1)
  if (length(FPcont)>0){
    internalFPidx<-c(1)
    tempidx<-2
    for (j in 1:length(FPcont)){
      internalFPidx[tempidx]<-FPcont[j]
      internalFPidx[tempidx+1]<-FPcont[j]+1
      tempidx<-tempidx+2
    }
    internalFPidx[tempidx]<-length(FalsePosIDX)
  } else {
    internalFPidx<-c(1,length(FalsePosIDX))
  }
  if (length(FNcont)>0){
    internalFNidx<-c(1)
    tempidx<-2
    for (j in 1:length(FNcont)){
      internalFNidx[tempidx]<-FNcont[j]
      internalFNidx[tempidx+1]<-FNcont[j]+1
      tempidx<-tempidx+2
    }
    internalFNidx[tempidx]<-length(FalseNegIDX)
  } else {
    internalFNidx<-c(1,length(FalseNegIDX))
  }
  FP<-0
  if (length(FalsePosIDX)>0){
    for (j in 1:(length(internalFPidx)/2)){
      jdxl<-(j-1)*2+1
      jdxh<-j*2
      fptemph<-FP_T_ecdf(dT$x[FalsePosIDX[internalFPidx[jdxh]]])
      fptempl<-FP_T_ecdf(dT$x[FalsePosIDX[internalFPidx[jdxl]]])
      FP<-FP+(fptemph-fptempl)
    }
  }
  
  FN<-0
  if (length(FalseNegIDX)>0){
    for (j in 1:(length(internalFNidx)/2)){
      jdxl<-(j-1)*2+1
      jdxh<-j*2
      fntemph<-FN_R_ecdf(dR$x[FalseNegIDX[internalFNidx[jdxh]]])
      fntempl<-FN_R_ecdf(dR$x[FalseNegIDX[internalFNidx[jdxl]]])
      FN<-FN+(fntemph-fntempl)
    }
  }
  Overlap<-FN+FP
  OverlapPerc<-Overlap/(2-Overlap)
  FP<-FP/(2-Overlap)
  FN<-FN/(2-Overlap)
  FP_ibounds<-FalsePosIDX[internalFPidx]
  FN_ibounds<-FalseNegIDX[internalFNidx]
  Output<-list('FP'=FP,'FN'=FN,'Overlap'=Overlap,'PercOverlap'=OverlapPerc,'FPi'=FP_ibounds,'FNi'=FN_ibounds)
  return(Output)
  
}

RelativeRank<-function(SimilarityObj,IntObjA,IntObjB){
  #ranking info used to scale
  k<-100
  #build density
  GPs<-row.names(SimilarityObj$RankInfoFinal)
  
  Output<-data.frame(matrix(NA,nrow=length(GPs),ncol=7))
  row.names(Output)<-GPs
  colnames(Output)<-c('Overlap','Alpha','Beta','MeanZ1','ActualZ1','MeanZ2','ActualZ2')
  
  for (j in 1:length(GPs)){
    TDis<-SimilarityObj$RankInfoFinal[GPs[j],]
    dT<-density(TDis,from=-0.1,to=1.1,na.rm=T)
    NDis1<-IntObjA$InternalRankingInfo[GPs[j],]
    dN1<-density(NDis1,from=-0.1,to=1.1,na.rm=T)
    NDis2<-IntObjA$InternalRankingInfo[GPs[j],]
    dN2<-density(NDis2,from=-0.1,to=1.1,na.rm=T)
    TAct<-SimilarityObj$RankInfoActual[GPs[j]]
    #assess
    if (max(TDis)!=0){
      CompPerc<-ecdf(TDis)(TAct)
      JointPerc1<-ecdf(NDis1)(TAct)
      JointPerc2<-ecdf(NDis2)(TAct)
      CompH<-approx(dT$x,k*dT$y/sum(dT$y),TAct)
      if (JointPerc1>JointPerc2){
        dN<-dN1
        NDis<-NDis1
      } else {
        dN<-dN2
        NDis<-NDis2
      }
      
      JointH<-approx(dN$x,k*dN$y/sum(dN$y),TAct)
      if (is.na(CompH$y)){
        CompH$y<-0
      }
      if (is.na(JointH$y)){
        JointH$y<-0
      }
      df <- merge(
        as.data.frame(dN[c("x", "y")]),
        as.data.frame(dT[c("x", "y")]),
        by = "x", suffixes = c(".T", ".N")
      )
      df$comp <- as.numeric(df$y.T > df$y.N)
      df$cross <- c(NA, diff(df$comp))
      CPoint<-which(df$cross!=0)[which(max(df$y.T[which(df$cross!=0)])==df$y.T[which(df$cross!=0)])]
      XPoint<-df$x[CPoint]
    } else {
      CPoint<-1
      df <- merge(
        as.data.frame(dN[c("x", "y")]),
        as.data.frame(dT[c("x", "y")]),
        by = "x", suffixes = c(".T", ".N")
      )
      df$y.T[CPoint]<-0
    }
    
    if (df$y.T[CPoint]<(10^-10)){
      Overlap<-0
      AlphaValue=0
      BetaValue=0
    } else {
      TArea<-(1-ecdf(TDis)(XPoint))
      NArea<-ecdf(NDis)(XPoint)
      Overlap<-round(TArea+NArea,2)
      BetaValue<-round(NArea/(2-Overlap),2)
      AlphaValue<-round(TArea/(2-Overlap),2)
    }
    Output[j,1]<-BetaValue+AlphaValue
    Output[j,2]<-AlphaValue
    Output[j,3]<-BetaValue
    
    Nmean1<-mean(NDis1,na.rm=T)
    Nsd1<-sd(NDis1,na.rm=T)
    zT1<-mean((TDis-Nmean1)/Nsd1,na.rm=T)
    Output[j,4]<-zT1
    Output[j,5]<-(TAct-Nmean1)/Nsd1
    Nmean2<-mean(NDis2,na.rm=T)
    Nsd2<-sd(NDis2,na.rm=T)
    zT2<-mean((TDis-Nmean2)/Nsd2,na.rm=T)
    Output[j,6]<-zT2
    Output[j,7]<-(TAct-Nmean2)/Nsd2
  }
  
  
  return(Output)
}


##Symmetric quality of contributions
ContributionCheckInternal<-function(SimObj,IntA,IntB,ks=F){
  IdentList<-colnames(SimObj$RankInfoActual)
  InternalCheck<-data.frame(matrix(NA,nrow=length(IdentList),ncol=3),row.names = IdentList)
  colnames(InternalCheck)<-c('Both','File1','File2')
  #1 is pass 0 is fail -1 is fails to appear in file
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    if (!(l %in% row.names(IntA$InternalRankingInfo))){
      InternalCheck[l,2]=-1
    } else {
      Med<-median(IntA$InternalRankingInfo[l,])<SimObj$RankInfoActual[l]
      if (ks==T){
        KSObj<-ks.test(IntA$InternalRankingInfo[l,],SimObj$RankInfoFinal[l,],alternative='g')
        KSVal<-KSObj$p.value<0.05
        InternalCheck[l,2]<-KSVal*Med
      } else if (ks==F){
        TDis<-SimObj$RankInfoFinal[l,]
        dT<-density(TDis,from=-0.05,to=1.05,na.rm=T)
        NDis<-IntA$InternalRankingInfo[l,]
        dN<-density(NDis,from=-0.05,to=1.05,na.rm=T)
        k<-100
        df <- merge(
          as.data.frame(dN[c("x", "y")]),
          as.data.frame(dT[c("x", "y")]),
          by = "x", suffixes = c(".T", ".N")
        )
        df$comp <- as.numeric(df$y.T > df$y.N)
        df$cross <- c(NA, diff(df$comp))
        CPoint<-which(df$cross!=0)[which(max(df$y.T[which(df$cross!=0)])==df$y.T[which(df$cross!=0)])]
        XPoint<-df$x[CPoint]
        if (df$y.T[CPoint]<(10^-10)){
          Overlap<-0
          
        } else {
          TArea<-k*(1-ecdf(TDis)(XPoint))
          NArea<-k*ecdf(NDis)(XPoint)
          Overlap<-round(TArea+NArea,2)/(100+100-round(TArea+NArea,2))
        }
        
        InternalCheck[l,2]<-(Overlap<=25)*Med
      }
      
      
    }
    if (!(l %in% row.names(IntB$InternalRankingInfo))){
      InternalCheck[l,3]=-1
    } else {
      Med<-median(IntB$InternalRankingInfo[l,])<SimObj$RankInfoActual[l]
      if (ks==T){
        KSObj<-ks.test(IntB$InternalRankingInfo[l,],SimObj$RankInfoFinal[l,],alternative='g')
        KSVal<-KSObj$p.value<0.05
        InternalCheck[l,3]<-KSVal*Med
      } else if (ks==F){
        TDis<-SimObj$RankInfoFinal[l,]
        dT<-density(TDis,from=-0.05,to=1.05,na.rm=T)
        NDis<-IntB$InternalRankingInfo[l,]
        dN<-density(NDis,from=-0.05,to=1.05,na.rm=T)
        k<-100
        df <- merge(
          as.data.frame(dN[c("x", "y")]),
          as.data.frame(dT[c("x", "y")]),
          by = "x", suffixes = c(".T", ".N")
        )
        df$comp <- as.numeric(df$y.T > df$y.N)
        df$cross <- c(NA, diff(df$comp))
        CPoint<-which(df$cross!=0)[which(max(df$y.T[which(df$cross!=0)])==df$y.T[which(df$cross!=0)])]
        XPoint<-df$x[CPoint]
        if (df$y.T[CPoint]<(10^-10)){
          Overlap<-0
        } else {
          TArea<-k*(1-ecdf(TDis)(XPoint))
          NArea<-k*ecdf(NDis)(XPoint)
          Overlap<-round(TArea+NArea,2)/(100+100-round(TArea+NArea,2))
        }
        
        InternalCheck[l,3]<-(Overlap<=25)*Med
      }
    }
    if ((InternalCheck[l,2]==0)&(InternalCheck[l,3]==0)){
      InternalCheck[l,1]<-0
    } else if ((InternalCheck[l,2]==-1)|(InternalCheck[l,3]==-1)){
      InternalCheck[l,1]=-1
    } else {
      InternalCheck[l,1]<-1
    }
  }
  return(list('Raw'=InternalCheck,'WinList'=IdentList[which(InternalCheck[,1]==1)],'FailList'=IdentList[which(InternalCheck[,1]!=1)]))
}

ContributionCheckComp<-function(SimObj,ExclusionList='Default'){
  IdentList<-colnames(SimObj$RankInfoActual)
  if (ExclusionList[1]!='Default'){
    for (j in 1:length(ExclusionList)){
      IdentList<-IdentList[-which(IdentList==ExclusionList[j])]
    }
  }
  CompCheck<-data.frame(matrix(NA,nrow=length(IdentList),ncol=3),row.names = IdentList)
  colnames(CompCheck)<-c('Both','Centrality','NonZero')
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    center<-ecdf(SimObj$RankInfoFinal[l,])(SimObj$RankInfoActual[l])
    nzed<-1-sum(SimObj$RankInfoFinal[l,]==0)/length(SimObj$RankInfoFinal[l,])
    CompCheck[l,2]<-((center<=0.75)&(center>=0.25))
    CompCheck[l,3]<-(nzed>=0.85)
    CompCheck[l,1]<-(CompCheck[l,2]&CompCheck[l,3])
  }
  return(list('Raw'=CompCheck,'WinList'=IdentList[which(CompCheck[,1]==T)],'FailList'=IdentList[which(CompCheck[,1]!=T)]))
}

ContributionCheckJoint<-function(SimObj,ExclusionList='Default',ks=F){
  IdentList<-colnames(SimObj$RankInfoActual)
  if (ExclusionList[1]!='Default'){
    for (j in 1:length(ExclusionList)){
      IdentList<-IdentList[-which(IdentList==ExclusionList[j])]
    }
  }
  JointCheck<-data.frame(matrix(NA,nrow=length(IdentList),ncol=3),row.names = IdentList)
  colnames(JointCheck)<-c('Both','LQuartile','Separation')
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    vic<-which(SimObj$NullRankInfoFinal[l,]==0)
    JointCheck[l,2]<-(ecdf(SimObj$NullRankInfoFinal[l,])(SimObj$RankInfoActual[l])<=0.25)
    if (ks==T){
      JointCheck[l,3]<-(ks.test(SimObj$RankInfoFinal[l,],SimObj$NullRankInfoFinal[l,],alternative='g')$p.value<=0.05)
    } else if(ks==F) {
      TDis<-SimObj$RankInfoFinal[l,]
      dT<-density(TDis,from=-0.05,to=1.05,na.rm=T)
      NDis<-SimObj$NullRankInfoFinal[l,]
      dN<-density(NDis,from=-0.05,to=1.05,na.rm=T)
      k<-100
      df <- merge(
        as.data.frame(dN[c("x", "y")]),
        as.data.frame(dT[c("x", "y")]),
        by = "x", suffixes = c(".T", ".N")
      )
      df$comp <- as.numeric(df$y.T > df$y.N)
      df$cross <- c(NA, diff(df$comp))
      CPoint<-which(df$cross!=0)[which(max(df$y.T[which(df$cross!=0)])==df$y.T[which(df$cross!=0)])]
      XPoint<-df$x[CPoint]
      if (df$y.T[CPoint]<(10^-10)){
        Overlap<-0
      } else {
        TArea<-k*(1-ecdf(TDis)(XPoint))
        NArea<-k*ecdf(NDis)(XPoint)
        Overlap<-round(TArea+NArea,2)/(100+100-round(TArea+NArea,2))
      }
      
      JointCheck[l,3]<-(Overlap<=25)
    }
    
    JointCheck[l,1]<-(JointCheck[l,2]&JointCheck[l,3])
  }
  return(list('Raw'=JointCheck,'WinList'=IdentList[which(JointCheck[,1]==1)],'FailList'=IdentList[which(JointCheck[,1]!=1)]))
}


RankingCriterion<-function(SimObj,ExclusionList='Default'){
  IdentList<-colnames(SimObj$RankInfoActual)
  if (ExclusionList!='Default'){
    for (j in 1:length(ExclusionList)){
      IdentList<-IdentList[-which(IdentList==ExclusionList[j])]
    }
  }
  RankInfo<-data.frame(rep(NA,length(IdentList)),row.names = IdentList)
  colnames(RankInfo)<-'Identifications'
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    RankInfo[l,1]<-mean(SimObj$RankInfoFinal[l,])*var(SimObj$RankInfoFinal[l,])*(abs(ecdf(SimObj$RankInfoFinal[l,])(SimObj$RankInfoActual[l][[1]]))+1/length(SimObj$RankInfoFinal[l,]))
  }
  return(RankInfo)
}

DQAndRank<-function(ObsR,TestR,NullR,InternalR1,InternalR2,RankMethod='Else'){
  #B12 is ranked mean RankInfo of file2 compared to distribution from file 1
  #B21 is ranked mean RankInfo of file1 compared to distribution from file 2
  RR1_2Hold<-RelativeRank(WithinRank1,BetweenRank12)
  RR1_2<-RR1_2Hold$RelRank
  RR2_1Hold<-RelativeRank(WithinRank2,BetweenRank21)
  RR2_1<-RR2_1Hold$RelRank
  RRTotal<-data.frame(matrix(NA,nrow=max(c(dim(WithinRank1$Summary)[1],dim(WithinRank1$Summary)[2])),ncol=2))
  colnames(RRTotal)<-c('R1_2','R2_1')
  RRTotal$R1_2<-RR1_2Hold$Summary
  RRTotal$R2_1<-RR2_1Hold$Summary
  R1_2List<-colnames(RR1_2)
  R2_1List<-colnames(RR2_1)
  Joint<-intersect(R1_2List,R2_1List)
  DFRank<-data.frame(matrix(0,nrow=2,ncol=length(Joint)))
  MSimL<-rep(0,length(Joint))
  SDSimL<-rep(0,length(Joint))
  colnames(DFRank)<-Joint
  rownames(DFRank)<-c('Mean','SD')
  QRank<-data.frame(matrix(0,nrow=8,ncol=length(Joint)))
  colnames(QRank)<-Joint
  rownames(QRank)<-c('Lower12','Lower21','Upper12','Upper21','IndicatorL','IndicatorU','IndicatorInvRL','IndicatorInvRU')
  
  #similarity of the two means and the two variances
  # multiply similarity by inverse of the means and the variances
  if (RankMethod=='Log'){
    lRR1_2<-log(RR1_2)
    minlR12<-min(lRR1_2[!is.infinite(lRR1_2)])
    lRR1_2[is.infinite(lRR1_2)]<-minlR12
    lRR2_1<-log(RR2_1)
    minlR21<-min(lRR2_1[!is.infinite(lRR2_1)])
    lRR2_1[is.infinite(lRR2_1)]<-minlR21
    for (j in 1:length(Joint)){
      idnm<-Joint[j]
      meanR12<-mean(lRR1_2[,idnm],na.rm=T)
      sdR12<-sd(lRR1_2[,idnm],na.rm=T)
      meanR21<-mean(lRR2_1[,idnm],na.rm=T)
      sdR21<-sd(lRR2_1[,idnm],na.rm=T)
      MeanSim<-meanR12*meanR21/(meanR12^2+meanR21^2-meanR12*meanR21)
      MSimL[j]<-MeanSim
      SDSim<-sdR12*sdR21/(sdR12^2+sdR21^2-sdR12*sdR21)
      SDSimL[j]<-SDSim
      DFRank[1,idnm]<-MeanSim*(1/mean(c(lRR1_2[,idnm],lRR2_1[,idnm]),na.rm=T)) #Sim*InvMean
      DFRank[2,idnm]<-SDSim*(1/sd(c(lRR1_2[,idnm],lRR2_1[,idnm]),na.rm=T)) #Sim*InvSD
      
      AmeanR12<-mean(RR1_2[,idnm],na.rm=T)
      AsdR12<-sd(RR1_2[,idnm],na.rm=T)
      AmeanR21<-mean(RR2_1[,idnm],na.rm=T)
      AsdR21<-sd(RR2_1[,idnm],na.rm=T)
      err12<-qnorm(0.975)*AsdR12/sqrt(dim(RR1_2)[1]-sum(is.na(RR1_2[,idnm])))
      err21<-qnorm(0.975)*AsdR21/sqrt(dim(RR2_1)[1]-sum(is.na(RR2_1[,idnm])))
      QRank[1,idnm]<-AmeanR12-err12
      QRank[2,idnm]<-AmeanR21-err21
      QRank[3,idnm]<-AmeanR12+err12
      QRank[4,idnm]<-AmeanR21+err21
      QRank[5,idnm]<-(QRank[1,idnm]<1 && QRank[2,idnm]<1)
      QRank[6,idnm]<-(QRank[3,idnm]<1 && QRank[4,idnm]<1)
      QRank[7,idnm]<-(QRank[1,idnm]<mean(RRTotal$R1_2,na.rm=T)^-1 && QRank[2,idnm]<mean(RRTotal$R2_1,na.rm=T)^-1)
      QRank[8,idnm]<-(QRank[3,idnm]<mean(RRTotal$R1_2,na.rm=T)^-1 && QRank[4,idnm]<mean(RRTotal$R2_1,na.rm=T)^-1)
    }
  } else if (RankMethod!='Log'){
    for (j in 1:length(Joint)){
      idnm<-Joint[j]
      meanR12<-mean(RR1_2[,idnm],na.rm=T)
      sdR12<-sd(RR1_2[,idnm],na.rm=T)
      meanR21<-mean(RR2_1[,idnm],na.rm=T)
      sdR21<-sd(RR2_1[,idnm],na.rm=T)
      MeanSim<-meanR12*meanR21/(meanR12^2+meanR21^2-meanR12*meanR21)
      MSimL[j]<-MeanSim
      SDSim<-sdR12*sdR21/(sdR12^2+sdR21^2-sdR12*sdR21)
      SDSimL[j]<-SDSim
      DFRank[1,idnm]<-MeanSim*(1/mean(c(RR1_2[,idnm],RR2_1[,idnm]),na.rm=T)) #Sim*InvMean
      DFRank[2,idnm]<-SDSim*(1/sd(c(RR1_2[,idnm],RR2_1[,idnm]),na.rm=T)) #Sim*InvSD
      err12<-qnorm(0.975)*sdR12/sqrt(dim(RR1_2)[1]-sum(is.na(RR1_2[,idnm])))
      err21<-qnorm(0.975)*sdR21/sqrt(dim(RR2_1)[1]-sum(is.na(RR2_1[,idnm])))
      QRank[1,idnm]<-meanR12-err12
      QRank[2,idnm]<-meanR21-err21
      QRank[3,idnm]<-meanR12+err12
      QRank[4,idnm]<-meanR21+err21
      QRank[5,idnm]<-(QRank[1,idnm]<1 && QRank[2,idnm]<1)
      QRank[6,idnm]<-(QRank[3,idnm]<1 && QRank[4,idnm]<1)
      QRank[7,idnm]<-(QRank[1,idnm]<mean(RRTotal$R1_2,na.rm=T)^-1 && QRank[2,idnm]<mean(RRTotal$R2_1,na.rm=T)^-1)
      QRank[8,idnm]<-(QRank[3,idnm]<mean(RRTotal$R1_2,na.rm=T)^-1 && QRank[4,idnm]<mean(RRTotal$R2_1,na.rm=T)^-1)
      
    }
  }
  
  
  
  Output<-list('RankData'=DFRank,'QualityCheck'=QRank,'MeanSim'=MSimL,'SDSim'=SDSimL,'RankTotal'=RRTotal)
  return(Output)
  
}

RankListMaker<-function(WtihinRank1,WithinRank2,BetweenRank12,BetweenRank21){
  #Take Aggregate Object and Produce list of IDs:
  ## are too low quality
  ## Ranking
  ## 
  
  AggregateObj<-RankAggregate(WithinRank1,WithinRank2,BetweenRank12,BetweenRank21)
  RRT12<-mean(AggregateObj$RankTotal[,1],na.rm = T)
  RRT21<-mean(AggregateObj$RankTotal[,2],na.rm = T)
  
  if (mean(AggregateObj$RankTotal[,1],na.rm=T)>1 & mean(AggregateObj$RankTotal[,2],na.rm=T)>1){
    print('Datasets 1 & 2 are both Below Overall Quality Threshold.')
  } else if (mean(AggregateObj$RankTotal[,1],na.rm=T)>1 & mean(AggregateObj$RankTotal[,2],na.rm=T)<=1){
    print('Dataset 2 is Below Quality Threshold; Dataset 1 may be a subset of the more broadly defined Dataset 2.')
  } else if (mean(AggregateObj$RankTotal[,1],na.rm=T)<=1 & mean(AggregateObj$RankTotal[,2],na.rm=T)>1){
    print('Dataset 1 is Below Quality Threshold; Dataset 2 may be a subset of the more broadly defined Dataset 1.')
  } else {
    print('Datasets 1 & 2 are both Above Overall Quality Threshold.')
    NumIDs<-dim(AggregateObj$QualityCheck)[2]
    LBound12<-AggregateObj$QualityCheck['Lower12',]<RRT12
    LBound21<-AggregateObj$QualityCheck['Lower21',]<RRT21
    for (j in 1:NumIDs){
      
    }
  }
  
  
}






#MissingValue Counter and Plotter
MVCountPlot<-function(filename,PlotTitle,HistTitle,AbunTitle,kmin=1,plotV=TRUE){
  file1<-SimDataClean(filename,kmin) #load file
  jlen<-dim(file1)[1] # glycopeptides
  jtot<-dim(file1)[2] # samples
  mvhold<-rep(0,jlen) # glycopeptides missing value holder
  meanhold<-rep(0,jlen)
  for (j in seq(jlen)){
    mvhold[j]<-sum(file1[j,]==0)/jtot
    meanhold[j]<-sum(file1[j,which(file1[j,]!=0)])/length(which(file1[j,]!=0))
  }
  if (plotV){
    par(mfrow=c(1,1))
    hist(mvhold,xlim=c(0,1),breaks = 20, main=HistTitle, xlab='Missing Value %',ylab='# Glycopeptides')
    glycopep<-row.names(file1)
    par(mfrow=c(2,1))
    plot(meanhold,1-mvhold,main=AbunTitle,xlab='Avg Relative Abundance',ylab='Presence',xlim=c(0,1),ylim=c(0,1))
    plot(file1[,1],1-mvhold,main=AbunTitle,xlab='Relative Abundance',ylab='Presence',xlim=c(0,1),ylim=c(0,1))
    for (j in 2:jtot){
      points(file1[,j],1-mvhold)
    }
    barplot(mvhold,main=PlotTitle,ylab='MV%',ylim=c(0,1), names.arg=glycopep)
    par(las=2)
  }
  par(mfrow=c(1,1))
  return(mvhold)
}

#MVPlottervAbundance
MVAbun<-function(filename,AbunTitle,kmin=1,rel=TRUE){
  file1<-SimDataClean(filename,kmin,rel)
  jlen<-dim(file1)[1]
  jtot<-dim(file1)[2]
  mvhold<-rep(0,jlen)
  meanhold<-rep(0,jlen)
  for (j in seq(jlen)){
    mvhold[j]<-sum(file1[j,]==0)/jtot
    meanhold[j]<-sum(file1[j,which(file1[j,]!=0)])/length(which(file1[j,]!=0))
  }
  if (rel){
    xl<-c(0,1)
  } else {
    xl<-c(0,max(file1))
  }
  yl<-c(0,1)
  par(mfrow=c(2,1))
  plot(meanhold,1-mvhold,main=AbunTitle,xlab='Avg Relative Abundance',ylab='Presence',xlim=xl,ylim=yl)
  plot(file1[,1],1-mvhold,main=AbunTitle,xlab='Relative Abundance',ylab='Presence',xlim=xl,ylim=yl)
  for (j in 2:jtot){
    points(file1[,j],1-mvhold)
  }
  par(las=2)
  par(mfrow=c(1,1))
}

RAMZISMain<-function(){
  #Sym comparison
  #Evaluation
  #internal quality
  #plots
  #Individual Quality
  #Ranking
  #Plots
}

#Function to separate glycopeptides by same peptide backbone
PeptideSegment<-function(filename1,filename2,kmin=2,simtype="Tanimoto",rel=TRUE,OverallName,SampleName1,SampleName2,mn=FALSE,logopt=FALSE){
  datfile1<-SimDataCleanLog(filename1,kmin=kmin,rel = rel)
  datfile2<-SimDataCleanLog(filename2,kmin=kmin,rel = rel)
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

#ModelComparison

