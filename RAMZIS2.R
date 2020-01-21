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
GlycReRead<-function(Filelist,NormVector){
  GFiles<-read.table(Filelist,stringsAsFactors = FALSE)
  GLen<-dim(GFiles)[1]
  NameHold<-rep(0,GLen)
  HoldTable<-data.frame(matrix(0,nrow = 1,ncol=(GLen+1)))
  for (j in 1:GLen){
    FileG<-GFiles[GLen,1]
    FNameA<-strsplit(FileG,'\\/')[length(strsplit(FileG,'\\/'))]
    FNameB<-strsplit(FNameA,'\\.csv')[length(strsplit(FNameA,'\\.csv'))]
    NameHold[j]<-FNameB
    GFile<-read.table(FileG,header=TRUE,row.names=1,stringsAsFactors = FALSE)
  }
}

GlyLineRead<-function(Filelist){}

#SimDataClean
SimDataClean<-function(filename,kmin=2,rel=TRUE,normvector='Default'){
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

TheoreticalDataGenerator<-function(n,g,p,w,alp=2,bet=2,maxim=TRUE,w0=2){
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
        
        FractionObj<-MASS:::.rat(round(datavals[j],w0))$rat
        FractionObj<-FractionObj/(10^(w[l]-6))
        bet0<-FractionObj[2]-FractionObj[1]
        alp0<-FractionObj[1]
        datagen[j,l]<-rbeta(1,alp0,bet0)
      }
    } else if (length(w)==g) {
      FractionObj<-MASS:::.rat(round(datavals[j],w0))$rat
      FractionObj<-FractionObj/(10^(w[j]-6))
      bet0<-FractionObj[2]-FractionObj[1]
      alp0<-FractionObj[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    } else {
      FractionObj<-MASS:::.rat(round(datavals[j],w0))$rat
      FractionObj<-FractionObj/(10^(w[1]-6))
      bet0<-FractionObj[2]-FractionObj[1]
      alp0<-FractionObj[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    }
    padj<-p*exp(mean(unlist(datagen[j,])-datmean))
    padj[padj>=1]<-0.99999
    padj[padj<=0]<-0.00001
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


#Bootstrap with n less loops
WithinSimBootstrap2<-function(filename,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,bootie=TRUE){
  # filename is file or matrix: if matrix be sure to normalize first
  # kmim, rel are simdataclean variables
  # MVCorrection only to be touched if you don't want to use A*GPPresence(A)
  # mn leave alone unless you make custom plotting function. delivers all withins rather than summary
  library('gtools')
  file1<-SimDataClean(filename,kmin,rel) #open file
  #acquire all glycopeptides in file
  glycopep1<-c(row.names(file1)) # list of glycopeptides
  #set up possible combinations
  coln<-ncol(file1) # sample num
  cols<-seq(coln) # list of sample nums
  Nsmpl<-CombWRep(coln,coln-1) #get number of n-1 samples with replacements
  combo<-data.frame(matrix(1,nrow=Nsmpl,ncol=(coln-1)))
  if (bootie){
    if (coln<8 & 2<coln){
      combo<-data.frame(gtools::combinations(coln,(coln-1),repeats.allowed = TRUE))
      rmlist<-c()
      for (j in 1:coln){
        rmlist<-c(rmlist,min(which(combo[,1]==j))) #removes first combination of each number which is singularities
      }
      combo<-combo[-rmlist,]
    } else if (coln<=2) {
      print(paste(filename,' has too few samples. This will not work.'))
      stop()
    } else {
      print('There be a number of samples here. The bootstrap generation will be slower.')
      combo<-data.frame(matrix(1,nrow=4290,ncol=(coln-1)))
      for (j in 1:4290){
        combo[j,]<-sort(sample(coln,(coln-1),replace = TRUE))
      }
      combo<-combo[-duplicated(combo),]
    }
  } else {
    combo<-combn(cols,(coln-1))
  }
  
  
  #jaccard and tanimoto trackers
  jac1<-rep(0,dim(combo)[1]*coln) # jaccard holder
  tan1<-rep(0,dim(combo)[1]*coln) # tanimoto holder
  tan1Final<-rep(0,dim(combo)[1])
  tanmathold<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycopep1)))
  colnames(tanmathold)<-glycopep1
  tanmatholdW<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycopep1)))
  colnames(tanmatholdW)<-glycopep1
  
  #iterate through combinations
  for (j in 1:nrow(combo)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    
    temp2<-file1 # subset of non combo data
    row.names(temp2)<-glycopep1 # GP names
    
    #acquire glycopeptides
    gl1<-row.names(temp1)
    gl2<-row.names(temp2)
    #compare for jaccard
    M11<-sum(gl1 %in% gl2)
    jac1[((j-1)*coln)]<-M11/length(glycopep1)
    
    #turns 0s into NAs if want to reduce impact of missing values
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
      temp2[temp2==0]<-NA
    }
    
    #generate T10, T01, T1_ and T_1
    if (ncol(temp1)==1){
      T1_<-data.frame(matrix(temp1[,1],nrow=length(gl1),ncol=1))
      row.names(T1_)<-gl1
      T10<-sum(temp1[,1]^2,na.rm = TRUE)
    } else {
      T1_<-data.frame(matrix(rowMeans(temp1,na.rm =TRUE),nrow=length(gl1),ncol=1))
      row.names(T1_)<-gl1
      T10<-sum(rowMeans(temp1,na.rm =TRUE)^2)
    }
    T_1<-temp2
    row.names(T_1)<-gl2
    T01<-colSums(temp2^2,na.rm=TRUE)
    
    
    
    #Bring T11 related terms together
    T__Hold<-merge(T1_,T_1,by=0,all=TRUE)
    row.names(T__Hold)<-T__Hold[,1]
    T__Hold<-T__Hold[,-1]
    T__Hold[is.na(T__Hold)]<-0
    
    dT<-sqrt((T__Hold[,1]-T__Hold[,-1])^2)
    
    PHold1<-data.frame(matrix(0,nrow=length(gl1),ncol=1))
    row.names(PHold1)<-gl1
    colnames(PHold1)<-'P1'
    PHold2<-data.frame(matrix(0,nrow=length(gl2),ncol=coln))
    row.names(PHold2)<-gl2
    colnames(PHold2)<-paste0(1:coln)
    for (m in 1:length(gl1)){
      nacheck<-1-sum(temp1[gl1[m],]==0)/length(temp1[gl1[m],])
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHold1[gl1[m],]<-nacheck
    }
    for (m in 1:length(gl2)){
      nacheck<-temp2[gl2[m],]>0
      nacheck[is.na(nacheck)]<-0
      PHold2[gl2[m],]<-nacheck
    }
    PHold<-merge(PHold1,PHold2,by=0,all=TRUE)
    row.names(PHold)<-PHold[,1]
    PHold<-PHold[,-1]
    PHold[is.na(PHold)]<-0
    presence<-(PHold[,1]+PHold[,-1])/2
    if (mn==FALSE){
      KTerm<-1+presence
    } else {
      KTerm<-mn
    }
    for (l in 1:coln){
      T11<-T__Hold[,1]*T__Hold[,2:(coln+1)]*(KTerm^(-dT))
      TRef<-row.names(T__Hold)
      for (m in 1:length(TRef)){
        tanmathold[((j-1)*coln+l),TRef[m]]<-T11[TRef[m],l]
        tanmatholdW[((j-1)*coln+l),TRef[m]]<-T11[TRef[m],l]/(T10+T01[l]-sum(T11[,l]))
      }
      tan1[((j-1)*coln+l)]<-sum(T11[,l])/(T10+T01[l]-sum(T11[,l]))
    }
    
    
  }
  
  Output<-data.frame(matrix(c(jac1,tan1),nrow=nrow(combo)*coln,ncol=2))
  colnames(Output)<-c('Jaccard','Tanimoto')
  tanmatendw<-data.frame(matrix(data=0,nrow=(dim(combo)[1]) ,ncol=length(glycopep1)))
  colnames(tanmatendw)<-glycopep1
  for (j in 1:dim(combo)[1]){
    tanmatendw[j,]<-colMeans(tanmatholdW[((j-1)*coln+1):(j*coln),],na.rm=TRUE)
    tan1Final[j]<-mean(tan1[((j-1)*coln+1):(j*coln)])
  }
  TaniOut<-data.frame(matrix(tan1Final,nrow=dim(combo)[1],ncol=1))
  colnames(TaniOut)<-c('Tanimoto')
  
  FinalOut<-list("Summary"=TaniOut,"Hold"=Output,"RankInfo"=tanmathold,"Boot"=combo,'RankInfoW'=tanmatholdW,'RankInfoFinal'=tanmatendw)
  return(FinalOut)
}

BetweenSimBootstrap2<-function(filename1,filename2,combo,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE){
  #load data and acquire glycopeptides
  file1<-SimDataClean(filename1,kmin,rel)
  file2<-SimDataClean(filename2,kmin,rel)
  glycopep1<-c(row.names(file1))
  glycopep2<-c(row.names(file2))
  glycojoint<-unique(c(glycopep1,glycopep2))
  mergedf<-merge(file1,file2,by=0,all=FALSE)
  row.names(mergedf)<-mergedf[,1]
  mergedf<-mergedf[,-1]
  coln<-dim(file2)[2]
  
  jac1<-rep(0,nrow(combo)*coln) # jaccard holder
  tan1<-rep(0,nrow(combo)*coln) # tanimoto holder
  tan1Final<-rep(0,dim(combo)[1])
  coln<-dim(file2)[2]
  tanmathold<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycojoint)))
  colnames(tanmathold)<-glycojoint
  tanmatholdW<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycojoint)))
  colnames(tanmatholdW)<-glycojoint
  
  
  #iterate through all combinations used in within of file 1
  for (j in 1:nrow(combo)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    temp2<-file2 # subset of non combo data
    row.names(temp2)<-glycopep2 # GP names
    
    #acquire glycopeptides
    gl1<-glycopep1
    gl2<-glycopep2
    #compare for jaccard
    M11<-sum(gl1 %in% gl2)
    jac1[((j-1)*coln)]<-M11/length(glycojoint)
    
    #turns 0s into NAs if want to reduce impact of missing values
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
      temp2[temp2==0]<-NA
    }
    
    #generate T10, T01, T1_ and T_1
    if (ncol(temp1)==1){
      T1_<-data.frame(matrix(temp1[,1],nrow=length(gl1),ncol=1))
      row.names(T1_)<-gl1
      T10<-sum(temp1[,1]^2,na.rm = TRUE)
    } else {
      T1_<-data.frame(matrix(rowMeans(temp1,na.rm =TRUE),nrow=length(gl1),ncol=1))
      row.names(T1_)<-gl1
      T10<-sum(rowMeans(temp1,na.rm =TRUE)^2)
    }
    T_1<-temp2
    row.names(T_1)<-gl2
    T01<-colSums(temp2^2,na.rm=TRUE)
    
    
    
    #Bring T11 related terms together
    T__Hold<-merge(T1_,T_1,by=0,all=TRUE)
    row.names(T__Hold)<-T__Hold[,1]
    T__Hold<-T__Hold[,-1]
    T__Hold[is.na(T__Hold)]<-0
    
    rowpres<-sum(T1_>0)
    distMod<-1/sum() #adjusts for presence
    dT<-sqrt((T__Hold[,1]-T__Hold[,-1])^2)
    
    PHold1<-data.frame(matrix(0,nrow=length(gl1),ncol=1))
    row.names(PHold1)<-gl1
    colnames(PHold1)<-'P1'
    PHold2<-data.frame(matrix(0,nrow=length(gl2),ncol=coln))
    row.names(PHold2)<-gl2
    colnames(PHold2)<-paste0(1:coln)
    for (m in 1:length(gl1)){
      nacheck<-1-sum(temp1[gl1[m],]==0)/length(temp1[gl1[m],])
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHold1[gl1[m],]<-nacheck
    }
    for (m in 1:length(gl2)){
      nacheck<-temp2[gl2[m],]>0
      nacheck[is.na(nacheck)]<-0
      PHold2[gl2[m],]<-nacheck
    }
    PHold<-merge(PHold1,PHold2,by=0,all=TRUE)
    row.names(PHold)<-PHold[,1]
    PHold<-PHold[,-1]
    PHold[is.na(PHold)]<-0
    presence<-(PHold[,1]+PHold[,-1])/2
    if (mn==FALSE){
      KTerm<-1+presence
    } else {
      KTerm<-mn
    }
    for (l in 1:coln){
      T11<-T__Hold[,1]*T__Hold[,2:(coln+1)]*(KTerm^(-dT))
      TRef<-row.names(T__Hold)
      for (m in 1:length(TRef)){
        tanmathold[((j-1)*coln+l),TRef[m]]<-T11[TRef[m],l]
        tanmatholdW[((j-1)*coln+l),TRef[m]]<-T11[TRef[m],l]/(T10+T01[l]-sum(T11[,l]))
      }
      tan1[((j-1)*coln+l)]<-sum(T11[,l])/(T10+T01[l]-sum(T11[,l]))
    }
    
    
  }
  
  Output<-data.frame(matrix(c(jac1,tan1),nrow=nrow(combo)*coln,ncol=2))
  colnames(Output)<-c('Jaccard','Tanimoto')
  tanmatendw<-data.frame(matrix(data=0,nrow=(dim(combo)[1]) ,ncol=length(glycojoint)))
  colnames(tanmatendw)<-glycojoint
  for (j in 1:dim(combo)[1]){
    tanmatendw[j,]<-colMeans(tanmatholdW[((j-1)*coln+1):(j*coln),],na.rm=TRUE)
    tan1Final[j]<-mean(tan1[((j-1)*coln+1):(j*coln)])
  }
  TaniOut<-data.frame(matrix(tan1Final,nrow=dim(combo)[1],ncol=1))
  colnames(TaniOut)<-c('Tanimoto')
  
  FinalOut<-list("Summary"=TaniOut,"Hold"=Output,"RankInfo"=tanmathold,"Boot"=combo,'RankInfoW'=tanmatholdW,'RankInfoFinal'=tanmatendw)
  return(FinalOut)
}

###Preserved Symmetry Comparison
SymmetricalSimBootstrap<-function(filename1,filename2,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,bootie=TRUE){
  #load data and acquire glycopeptides
  file1<-SimDataCleanLog(filename1,kmin,rel)
  file2<-SimDataCleanLog(filename2,kmin,rel)
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
    if (coln2<5 & 2<coln2){
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
  
InternalSimilarity<-function(filename,BootSet,kmin=1,rel=TRUE,MVCorrection=TRUE,mn=FALSE){
  file1<-SimDataClean(filename,kmin,rel)
  glycopep1<-c(row.names(file1))
  glycojoint<-glycopep1
  glycopep2<-glycopep1
  coln1<-dim(file1)[2]
  cols1<-seq(coln1)
  Nsmpl1<-CombWRep(coln1,coln1-1)
  combo1<-BootSet
  ncomps<-nrow(combo1)*(nrow(combo1)-1)/2
  
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
  #iterate through all combinations used in within of file 2
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
    tracker1<-tracker1+ncomb1-j
  }
  tan1<-colSums(tanmatholdW)
  
  
  OutputObj<-list("InternalTanimoto"=tan1,"InternalRankingInfo"=tanmatholdW)
  return(OutputObj)
}

#plots internal boot object[[]]: 1 is the first
InternalQuality<-function(filename,BootSet,SimilarityObj,PlotTitle,GroupName,Int=NULL,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,verbose=FALSE,legendpos='topleft'){
  if (!is.null(Int)){
    #Use Int as Int object
  } else {
    Int<-InternalSimilarity(filename,BootSet,kmin,rel,MVCorrection,mn)
  }
  
  k<-100
  #build density
  TDis<-SimilarityObj$Summary$Tanimoto
  dT<-density(TDis,from=-0.1,to=1.1,na.rm=T)
  NDis<-SimilarityObj$NullOut$NullTani
  dN<-density(NDis,from=-0.1,to=1.1,na.rm=T)
  IDis<-Int$InternalTanimoto
  dI<-density(IDis,from=-0.2,to=1.2,na.rm=T)
  TAct<-SimilarityObj$Actual
  #plot densities
  if (max(TDis)==0){
    mh<-max(c(k*dN$y/sum(dN$y)))
    plot(c(0,0,0.1,0.1),c(0,mh,mh,0),xlim=c(0,1),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution')
    polygon(c(0,0,0.1,0.1),c(0,mh,mh,0),col=rgb(1,0,0,0.5))
  } else {
    mh<-max(c(k*dT$y/sum(dT$y),k*dN$y/sum(dN$y)))
    plot(dT$x,k*dT$y/sum(dT$y),xlim=c(0,1),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='% of Distribution')
    polygon(dT$x,k*dT$y/sum(dT$y),col=rgb(1,0,0,0.5))
  }
  lines(dN$x,k*dN$y/sum(dN$y))
  polygon(dN$x,k*dN$y/sum(dN$y),col=rgb(0,0,1,0.5))
  lines(dI$x,k*dI$y/sum(dI$y))
  polygon(dI$x,k*dI$y/sum(dI$y),col=rgb(0,1,0,0.5))
  lines(rep(TAct,2),c(0,100),col=1,lwd=3)
  #percentile location
  CompPerc<-ecdf(TDis)(TAct)
  JointPerc<-ecdf(NDis)(TAct)
  CompH<-approx(dT$x,k*dT$y/sum(dT$y),TAct)
  JointH<-approx(dN$x,k*dN$y/sum(dN$y),TAct)
  if (is.na(CompH$y)){
    CompH$y<-0
  }
  if (is.na(JointH$y)){
    JointH$y<-0
  }
  #text(TAct,CompH$y,labels = paste0(round(k*CompPerc,2),'% of Test'),pos=4)
  #text(TAct,JointH$y,labels = paste0(round(k*JointPerc,2),'% of Null'),pos=2)
  #coverage overlap
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
    Overlap<-round(TArea+NArea,2)
  }
  legend(legendpos,legend=c('Observed Similarity',paste('Internal of',GroupName),'Test Distribution','Null Distribution',paste0(Overlap,'% Overlap of Distributions')),fill=c(1,3,2,4,rgb(1,0,1)))
  if (verbose==T){
    return(list(Int,Overlap,XPoint,CompPerc,JointPerc))
  }

  
}

#plots observed, test, null

SimPlot<-function(PlotTitle,SimilarityObj,legendpos='topleft',verbose=F){
  k<-100
  #build density
  TDis<-SimilarityObj$Summary$Tanimoto
  dT<-density(TDis,from=-0.1,to=1.1,na.rm=T)
  NDis<-SimilarityObj$NullOut$NullTani
  dN<-density(NDis,from=-0.1,to=1.1,na.rm=T)
  TAct<-SimilarityObj$Actual
  #plot densities
  if (max(TDis)==0){
    mh<-max(c(k*dN$y/sum(dN$y)))
    plot(c(0,0,0.05,0.05),c(0,mh,mh,0),xlim=c(0,1),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='Density')
    polygon(c(0,0,0.05,0.05),c(0,mh,mh,0),col=rgb(1,0,0,0.5))
  } else {
    mh<-max(c(k*dT$y/sum(dT$y),k*dN$y/sum(dN$y)))
    plot(dT$x,k*dT$y/sum(dT$y),xlim=c(0,1),ylim=c(0,mh),main=PlotTitle,type='l',xlab='Similarity',ylab='Density')
    polygon(dT$x,k*dT$y/sum(dT$y),col=rgb(1,0,0,0.5))
  }
  lines(dN$x,k*dN$y/sum(dN$y))
  polygon(dN$x,k*dN$y/sum(dN$y),col=rgb(0,0,1,0.5))
  lines(rep(TAct,2),c(0,100),col=1,lwd=3)
  #percentile location
  if (max(TDis)!=0){
    CompPerc<-ecdf(TDis)(TAct)
    JointPerc<-ecdf(NDis)(TAct)
    CompH<-approx(dT$x,k*dT$y/sum(dT$y),TAct)
    JointH<-approx(dN$x,k*dN$y/sum(dN$y),TAct)
    if (is.na(CompH$y)){
      CompH$y<-0
    }
    if (is.na(JointH$y)){
      JointH$y<-0
    }
    # text(TAct,CompH$y,labels = paste0(round(k*CompPerc,2),'% of Test'),pos=4)
    # text(TAct,JointH$y,labels = paste0(round(k*JointPerc,2),'% of Null'),pos=2)
    #coverage overlap
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
    df$Y.T[CPoint]<-0
  }
  
  if (df$y.T[CPoint]<(10^-10)){
    Overlap<-0
    legend(legendpos,legend=c('Observed Similarity','Test Distribution','Null Distribution','0 = Alpha','0 = Beta'),fill=c(1,rgb(1,0,0,0.5),rgb(0,0,1,0.5),'darkgray','black'))
    AlphaValue=0
    BetaValue=0
  } else {
    TArea<-(1-ecdf(TDis)(XPoint))
    NArea<-ecdf(NDis)(XPoint)
    Overlap<-round(TArea+NArea,2)
    AlphaValue<-round(NArea/(2-Overlap),2)
    BetaValue<-round(TArea/(2-Overlap),2)
    #lines(rep(XPoint,2),c(0,df$y.T[CPoint]))
    #make alpha area
    polygon(dN$x[dN$x<=XPoint],c(k*dN$y[dN$x<XPoint]/sum(dN$y),0),col='darkgray')
    #make beta area
    polygon(dT$x[dT$x>=XPoint],c(0,k*dT$y[dT$x>XPoint]/sum(dT$y)),col='black')
    legend(legendpos,legend=c('Observed Similarity','Test Distribution','Null Distribution',paste0(AlphaValue,'= Alpha'),paste0(BetaValue,'= Beta')),fill=c(1,rgb(1,0,0,0.5),rgb(0,0,1,0.5),'darkgray','black'))
  }
  
  
  if (verbose==T){
    return(list(Overlap,XPoint,CompPerc,JointPerc,AlphaValue,BetaValue))
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

RelativeRank<-function(WithinRankObj,BetweenRankObj){
  #ranking info used to scale
  WR<-WithinRankObj$RankInfoFinal
  BR<-BetweenRankObj$RankInfoFinal
  WList<-colnames(WR)
  BList<-colnames(BR)
  Joint<-intersect(WList,BList)
  RRank<-data.frame(matrix(0,nrow=dim(WR)[1],ncol=length(Joint)))
  colnames(RRank)<-Joint
  for (j in 1:length(Joint)){
    idnm<-Joint[j]
    RRank[,idnm]<-BR[,idnm]/WR[,idnm]
  }
  Sumry<-BetweenRankObj$Summary$Tanimoto/WithinRankObj$Summary$Tanimoto
  Output<-list("RelRank"=RRank,"WNotInB"=WList[!WList %in% BList],"BNotInW"=BList[!BList %in% WList],"Summary"=Sumry)
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
PeptideSegment<-function(filename1,filename2,kmin=2,simtype="Tanimoto",rel=TRUE,OverallName,SampleName1,SampleName2,mn=FALSE){
  datfile1<-SimDataClean(filename1,kmin=kmin,rel = rel)
  datfile2<-SimDataClean(filename2,kmin=kmin,rel = rel)
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

