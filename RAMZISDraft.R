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
#RAMZIS
#  df1=GlycReRead(FileList1,NormVector1) #load files from list 1 and format
#  df2=GlycReRead(FileList2,NormVector2) #load files from list 2 and format
#  #Comparative Similarity By Group
#  PeptideSegment(df1,df2,kmin,simtype,rel,OverallName,SampleName1,SampleName2)
#    LoopThroughUniquePeptideBackbones:
#      WithinS1=WithinSim(Subset(df1),kmin)
#        BetweenSim(sample(Subset(df1)),sample(Subset(df1)),kmin) 
#      WithinS2=WithinSim(Subset(df2),kmin)
#        BetweenSim(sample(Subset(df2)),sample(Subset(df2)),kmin)
#      Between12=BetweenSim(Subset(df1),Subset(df2),kmin)
#      Simplot(WithinS1,WithinS2,Between12) #outputs plot of comparative similarity
#        Boxplot(WithinS1,WithinS2,main=paste(OverallName,Backbone),names=c(SampleName1,SampleName2))
#		 Abline(Between12)
#        Labeling
#      EndPlot
#      RAMZISPeptideBackbone(WithinS1,WithinS2,Between12) #outputs three ranking lists
#      Storage<-WithinS1,WithinS2,Between12
#	   StoragePBLists<-RAMZISPBLists
#    EndLoop
#  RAMZISProtein(StoragePBListsW1,StoragePBListsW2,StoragePBListsB)
#  StoragePrLists<-RAMZISPrLists
#  RAMZISProteome(StoragePrListsW1,StoragePrListsW2,StoragePrListsB)
#  StoragePmLists<-RAMZISPmLists
#  Output: StoragePBLists, StoragePrLists, StoragePmLists
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

#SimDataClean
SimDataClean<-function(filename,kmin=2,rel=TRUE){
  #Read in Data and measure dataframe size
  if (typeof(filename)=='character'){
    file1<-read.csv(filename,header=TRUE, row.names=1,stringsAsFactors = FALSE)
    data1<-file1
    #normalize data by max abundance per sample
    if (rel){
      for (i in 1:ncol(file1)){
        data1[,i]<- file1[,i]/max(file1[,i])
      }
    }
  } else if(typeof(filename)=='list'){
    file1<-filename
    data1<-filename
    #print('Program Detected List for filename. Attempting use as dataframe')
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

TheoreticalDataGenerator<-function(n,g,p,w,alp=2,bet=2){
  #n samples
  #g identifications
  #p average presence
  ##p determines the presence distribution that appears across the sample
  ## binomial distribution rate of observation is proportional to signal intensity
  ## eg P(1,n,p_adj)=#observations. p_adj=p*exp(x-mean(x))
  #w term controlling the variance in the data
  #alp and bet are terms for general data distribution
  datagen<-data.frame(matrix(0,nrow=g,ncol=n))
  rowsampname<-paste0('Identification',1:g)
  row.names(datagen)<-rowsampname
  
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
  
  #use rounding to generate betas for each value so they all are within 0-1
  datmean<-mean(datavals)
  for (j in 1:g){
    if (length(w)==n){
      for (l in 1:n){
        bet0<-10^w[l]
        alp0<-round(datavals[j],w[l])*10^w[l]
        datagen[j,l]<-rbeta(1,alp0,bet0)
        
      }
    } else if (length(w)==g) {
      bet0<-10^w[j]
      alp0<-round(datavals[j],w[j])*10^w[j]
      datagen[j,]<-rbeta(n,alp0,bet0)
    } else {
      bet0<-10^w[1]
      alp0<-round(datavals[j],w[1])*10^w[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    }
    padj<-p*exp(mean(unlist(datagen[j,])-datmean))
    keepnum<-rbinom(1,(n-1),padj)+1
    datagen[j,-match(sample(datagen[j,],keepnum),datagen[j,])]<-0
  }
  datagen[which(rowMeans(datagen)==max(rowMeans(datagen))),]<-rep(1,n)
  FinalOut<-list("Dataset"=datagen,"parameters"=c(n,g,p),"w"=w,"alp"=alp,"bet"=bet,"DistributionCenters"=datavals)
  return(FinalOut)
}

# Function not worth it:
# DataPairGenerator<-function(a1,b1,a2,b2,w1,w2,s1,s2,p1,p2,n1,n2,g1,g2,sharerate){
#   #a1b1,w1 and a2b2,w2 are parameters controlling the internal similarity of the data
#   #s1 and s2 are parameters controlling the relatedness of the groups
#   #p1 and p2 are the average presence for the group
#   #n is the number of samples
#   #g is the number of identifications
#   #sharerate is the rate of identifications in common between the two groups
#   
#   g1adj<-ceiling(g1*(1/sharerate))
#   g2adj<-ceiling(g2*(1/sharerate))
#   
#   dat1<-TheoreticalDataGenerator(n1,g1adj,p1,w1,a1,b1)
#   
#   
# }

#WithinSim
WithinSim<-function(filename,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE){
  # filename is file or matrix: if matrix be sure to normalize first
  # kmim, rel are simdataclean variables
  # MVCorrection only to be touched if you don't want to use A*GPPresence(A)
  # mn leave alone unless you make custom plotting function. delivers all withins rather than summary
  file1<-SimDataClean(filename,kmin,rel) #open file
  #acquire all glycopeptides in file
  glycopep1<-c(row.names(file1)) # list of glycopeptides
  #set up possible combinations
  coln<-ncol(file1) # sample num
  cols<-seq(coln) # list of sample nums
  combo<-combn(cols,floor(coln/2)) #possible combination
  #jaccard and tanimoto trackers
  jac1<-rep(0,ncol(combo)) # jaccard holder
  tan1<-rep(0,ncol(combo)) # tanimoto holder
  
  #iterate through combinations
  for (j in 1:ncol(combo)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,combo[,j]]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    temp2<-data.frame(file1[,-combo[,j]]) # subset of non combo data
    row.names(temp2)<-glycopep1 # GP names
    #reduce to only present glycopeptides
    if (length(which(apply(temp1,1,sum)==0))>0){
      ghold<-row.names(temp1)[-which(apply(temp1,1,sum)==0)]
      temp1<-data.frame(temp1[-which(apply(temp1,1,sum)==0),])
      row.names(temp1)<-ghold
    }
    if (length(which(apply(temp2,1,sum)==0))>0){
      ghold<-row.names(temp2)[-which(apply(temp2,1,sum)==0)]
      temp2<-data.frame(temp2[-which(apply(temp2,1,sum)==0),])
      row.names(temp2)<-ghold
    }
    #acquire glycopeptides
    gl1<-row.names(temp1)
    gl2<-row.names(temp2)
    #compare for jaccard
    M11<-sum(gl1 %in% gl2)
    jac1[j]<-M11/length(glycopep1)
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
      temp2[temp2==0]<-NA
    }
    if (ncol(temp1)==1){
      T1_<-temp1[which(gl1 %in% gl2),1]
      T10<-sum(temp1[,1]^2,na.rm = TRUE)
    } else {
      T1_<-rowMeans(temp1[which(gl1 %in% gl2),],na.rm =TRUE)
      T10<-sum(rowMeans(temp1,na.rm =TRUE)^2)
    }
    if (ncol(temp2)==1){
      T_1<-temp2[which(gl2 %in% gl1),1]
      T01<-sum(temp2[,1]^2,na.rm = TRUE)
    } else {
      T_1<-rowMeans(temp2[which(gl2 %in% gl1),],na.rm =TRUE)
      T01<-sum(rowMeans(temp2,na.rm =TRUE)^2)
    }
    if (length(which(gl1 %in% gl2))!=0){
      for (k in 1:length(T1_)){
        if (is.na(T1_)[k]){
          T1_[k]<-0
        }
      }
      for (k in 1:length(T_1)){
        if (is.na(T_1)[k]){
          T_1[k]<-0
        }
      }
      if (is.na(T10)){
        T10<-0
      }
      if (is.na(T01)){
        T01<-0
      }
    } else {
      T_1<-0
      T1_<-0
    }
    
    T11<-sum(T1_*T_1)
    dT1_<-sqrt(sum(T1_^2))
    dT_1<-sqrt(sum(T_1^2))
    tan1[j]<-T11/(T10+T01-T11)
    if (dT1_>=dT_1){
      dT<-sqrt(sum((T_1-T1_)^2))/dT1_
    } else if (dT1_<dT_1){
      dT<-sqrt(sum((T_1-T1_)^2))/dT_1
    }
    pres1<-GPPresence(temp1,kmin=1,rel=rel)
    pres2<-GPPresence(temp2,kmin=1,rel=rel) #change presence function to tag with row names to line up,
    #fix program to make sure it is comparing gp1 to gp1
    presm<-mean(c(pres1,pres2))
    if (tan1[j]!=0){
      tan1[j]<-tan1[j]/((1+presm)^dT)
    }
    
  }
  if (mn==FALSE){
    Output<-data.frame(matrix(c(jac1,tan1),nrow=ncol(combo),ncol=2))
    colnames(Output)<-c('Jaccard','Tanimoto SD')
    
  } else {
    Output<-data.frame(matrix(c(mean(jac1),sd(jac1),mean(tan1),sd(tan1)),nrow=1,ncol=4))
    colnames(Output)<-c('Jaccard Mean','Jaccard SD','Tanimoto Mean','Tanimoto SD')
    
  }
  return(Output)
}

CombWRep<-function(n,r){
  return(factorial(n+r-1)/(factorial(n-1)*factorial(r)))
}

#WithinSimBootstrapSingleComp
WithinSimBootstrap1<-function(filename,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,bootie=TRUE){
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
    if (coln<10 & 2<coln){
      combo<-data.frame(gtools::combinations(coln,(coln-1),repeats.allowed = TRUE))
      rmlist<-c()
      for (j in 1:coln){
        rmlist<-c(rmlist,min(which(combo[,1]==j))) #removes first combination of each number which is singularities
      }
      combo<-combo[-rmlist,]
    } else if (coln<2) {
      print(paste(filename,' has too few samples. This will not work.'))
      stop()
    } else {
      print('There be a number of samples here. The bootstrap generation will be slower.')
      combo<-data.frame(matrix(1,nrow=4290,ncol=(coln-1)))
      for (j in 1:12870){
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
  tanmathold<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycopep1)))
  colnames(tanmathold)<-glycopep1
    
  #iterate through combinations
  for (j in 1:nrow(combo)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    for (l in 1:coln){
      temp2<-data.frame(file1[,l]) # subset of non combo data
      row.names(temp2)<-glycopep1 # GP names
      
      #acquire glycopeptides
      gl1<-row.names(temp1)
      gl2<-row.names(temp2)
      #compare for jaccard
      M11<-sum(gl1 %in% gl2)
      jac1[((j-1)*coln+l)]<-M11/length(glycopep1)
      
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
      if (ncol(temp2)==1){
        T_1<-data.frame(matrix(temp2[,1],nrow=length(gl2),ncol=1))
        row.names(T_1)<-gl2
        T01<-sum(temp2[,1]^2,na.rm = TRUE)
      } else {
        T_1<-data.frame(matrix(rowMeans(temp2,na.rm =TRUE),nrow=length(gl2),ncol=1))
        row.names(T_1)<-gl2
        T01<-sum(rowMeans(temp2,na.rm =TRUE)^2)
      }
      ### Replace above code with commented code below
      # T1_<-data.frame(matrix(rowMeans(temp1,na.rm =TRUE),nrow=length(gl1),ncol=1))
      # row.names(T1_)<-gl1
      # T10<-sum(rowMeans(temp1,na.rm =TRUE)^2)
      # 
      # T_1<-data.frame(matrix(temp2[,1],nrow=length(gl2),ncol=1))
      # row.names(T_1)<-gl2
      # T01<-sum(temp2[,1]^2,na.rm = TRUE)
      ### Can cut runtime in fourth if compare all samples at once. eg vecA vs matrixB
      ### Use colsums in T01
      ### use temp2 to hold all values
      ### 
      
      
      #if something in common, remove nas, else drop values to zero
      if (length(which(gl1 %in% gl2))!=0){
        for (k in 1:length(T1_)){
          if (is.na(T1_)[k]){
            T1_[k]<-0
          }
        }
        for (k in 1:length(T_1)){
          if (is.na(T_1)[k]){
            T_1[k]<-0
          }
        }
        if (is.na(T10)){
          T10<-0
        }
        if (is.na(T01)){
          T01<-0
        }
      } else {
        T_1<-0
        T1_<-0
      }
      
      #Bring T11 related terms together
      T__Hold<-merge(T1_,T_1,by=0,all=TRUE)
      row.names(T__Hold)<-T__Hold[,1]
      T__Hold<-T__Hold[,-1]
      T__Hold[is.na(T__Hold)]<-0
      
      
      
      dT1_<-sqrt(sum(T__Hold[,1]^2))
      dT_1<-sqrt(sum(T__Hold[,2]^2))
      if (dT1_>=dT_1){
        dT<-sqrt(sum((T__Hold[,1]-T__Hold[,2])^2))/dT1_
      } else if (dT1_<dT_1){
        dT<-sqrt(sum((T__Hold[,1]-T__Hold[,2])^2))/dT_1
      }
      PHold1<-data.frame(matrix(0,nrow=length(gl1),ncol=1))
      row.names(PHold1)<-gl1
      PHold2<-data.frame(matrix(0,nrow=length(gl2),ncol=1))
      row.names(PHold2)<-gl2
      for (m in 1:length(gl1)){
        nacheck<-1-sum(temp1[gl1[m],]==0)/length(temp1[gl1[m],])
        if (is.na(nacheck)){
          nacheck<-0
        }
        PHold1[gl1[m]]<-nacheck
      }
      for (m in 1:length(gl2)){
        nacheck<-1-sum(temp2[gl2[m],]==0)/length(temp2[gl2[m],])
        if (is.na(nacheck)){
          nacheck<-0
        }
        PHold2[gl2[m]]<-nacheck
      }
      PHold<-merge(PHold1,PHold2,by=0,all=TRUE)
      row.names(PHold)<-PHold[,1]
      PHold<-PHold[,-1]
      PHold[is.na(PHold)]<-0
      presence<-mean(rowMeans(PHold,na.rm = FALSE))
      if (mn==FALSE){
        KTerm<-1+presence
      } else {
        KTerm<-mn
      }
      T11<-T__Hold[,1]*T__Hold[,2]*(KTerm^(-dT))
      TRef<-row.names(T__Hold)
      for (m in 1:length(T11)){
        tanmathold[((j-1)*coln+l),TRef[m]]<-T11[m]
      }
      
      tan1[((j-1)*coln+l)]<-sum(T11)/(T10+T01-sum(T11))
    }
  }

  Output<-data.frame(matrix(c(jac1,tan1),nrow=nrow(combo)*coln,ncol=2))
  colnames(Output)<-c('Jaccard','Tanimoto')

  FinalOut<-list("Summary"=Output,"RankInfo"=tanmathold,"Boot"=combo)
  return(FinalOut)
}

#Bootstrap with n less loops
WithinSimBootstrap2<-function(filename,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,bootie=TRUE){
  # filename is file or matrix: if matrix be sure to normalize first
  # kmim, rel are simdataclean variables
  # MVCorrection only to be touched if you don't want to use A*GPPresence(A)
  # mn now controls K term 
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
    if (coln<10 & 2<coln){
      combo<-data.frame(gtools::combinations(coln,(coln-1),repeats.allowed = TRUE))
      rmlist<-c()
      for (j in 1:coln){
        rmlist<-c(rmlist,min(which(combo[,1]==j))) #removes first combination of each number which is singularities
      }
      combo<-combo[-rmlist,]
    } else if (coln<2) {
      print(paste(filename,' has too few samples. This will not work.'))
      stop()
    } else {
      print('There be a number of samples here. The bootstrap generation will be slower.')
      combo<-data.frame(matrix(1,nrow=4290,ncol=(coln-1)))
      for (j in 1:12870){
        combo[j,]<-sort(sample(coln,(coln-1),replace = TRUE))
      }
      combo<-combo[-duplicated(combo),]
    }
  } else {
    combo<-combn(cols,(coln-1))
  }
  
  
  #tanimoto trackers
  tan1<-rep(0,dim(combo)[1]*coln) # tanimoto holder
  tanmathold<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycopep1)))
  colnames(tanmathold)<-glycopep1
  
  
  #iterate through combinations
  for (j in 1:nrow(combo)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    temp2<-data.frame(file1) # temp datafile
    row.names(temp2)<-glycopep1 # GP names
    gl1<-row.names(temp1)
    gl2<-row.names(temp2)
    
    
    #turns 0s into NAs if want to reduce impact of missing values
    if (MVCorrection!=TRUE){
      temp1[temp1==0]<-NA
      temp2[temp2==0]<-NA
    }
    
    ## Replace above code with commented code below
    T1_<-rowMeans(temp1,na.rm =TRUE)
    T10<-sum(rowMeans(temp1,na.rm =TRUE)^2)

    T11pre<-temp2*T1_
    row.names(T11pre)<-glycopep1
    T11pre[is.na(T11pre)]<-0
    T01<-colSums(temp2^2,na.rm = TRUE)
    
    
    
    dT1_<-sqrt(sum(T1_^2))
    dT_1<-sqrt(colSums(temp2^2))
    dT<-dT_1
    dT[dT<dT1_]<-dT1_
    PHold1<-data.frame(matrix(0,nrow=length(glycopep1),ncol=coln))
    PHold2<-data.frame(matrix(0,nrow=length(glycopep1),ncol=coln))
    row.names(PHold2)<-gl2
    for (m in 1:length(gl1)){
      nacheck<-1-sum(temp1[gl1[m],]==0)/length(temp1[gl1[m],])
      if (is.na(nacheck)){
        nacheck<-0
      }
      PHold1[m,]<-rep(nacheck,coln)
    }
    PHold2[temp2==0]<-0
    PHold2[temp2!=0]<-1
    PHold<-(PHold2+PHold1)/2
    
    PHold[is.na(PHold)]<-0
    
    
    presence<-colMeans(PHold,na.rm = FALSE)
    if (mn==FALSE){
      KTerm<-rep(1,coln)+presence
    } else {
      KTerm<-rep(mn,coln)
    }
    for (l in 1:coln){
      T11<-T11pre[,l]*(KTerm[l]^-dT[l])
      TRef<-row.names(T11pre)
      for (m in 1:length(T11)){
        tanmathold[((j-1)*coln+l),TRef[m]]<-T11[m]
      }
      tan1[((j-1)*coln+l)]<-sum(T11)/(T10+T01[l]-sum(T11))
    }
    
    
    
  }
  
  Output<-data.frame(matrix(c(rep(0,length(tan1)),tan1),nrow=nrow(combo)*coln,ncol=2))
  colnames(Output)<-c('Tanimoto')
  
  FinalOut<-list("Summary"=Output,"RankInfo"=tanmathold,"Boot"=combo)
  return(FinalOut)
}

BetweenSimBootstrap<-function(filename1,filename2,combo,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE){
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
  coln<-dim(file2)[2]
  tanmathold<-data.frame(matrix(data=0,nrow=(dim(combo)[1]*coln) ,ncol=length(glycojoint)))
  colnames(tanmathold)<-glycojoint
  
  #iterate through all combinations used in within of file 1
  for (j in 1:nrow(combo)){
    #separate datasets for combinations
    temp1<-data.frame(file1[,unlist(combo[j,])]) # subset of first data
    row.names(temp1)<-glycopep1 # GP names
    for (l in 1:coln){
      temp2<-data.frame(file2[,l]) # subset of non combo data
      row.names(temp2)<-glycopep2 # GP names
      
      #acquire glycopeptides
      gl1<-row.names(temp1)
      gl2<-row.names(temp2)
      #compare for jaccard
      M11<-sum(gl1 %in% gl2)
      jac1[((j-1)*coln+l)]<-M11/length(glycojoint)
      
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
      if (ncol(temp2)==1){
        T_1<-data.frame(matrix(temp2[,1],nrow=length(gl2),ncol=1))
        row.names(T_1)<-gl2
        T01<-sum(temp2[,1]^2,na.rm = TRUE)
      } else {
        T_1<-data.frame(matrix(rowMeans(temp2,na.rm =TRUE),nrow=length(gl2),ncol=1))
        row.names(T_1)<-gl2
        T01<-sum(rowMeans(temp2,na.rm =TRUE)^2)
      }
      
      #if something in common, remove nas, else drop values to zero
      if (length(which(gl1 %in% gl2))!=0){
        for (k in 1:length(T1_)){
          if (is.na(T1_)[k]){
            T1_[k]<-0
          }
        }
        for (k in 1:length(T_1)){
          if (is.na(T_1)[k]){
            T_1[k]<-0
          }
        }
        if (is.na(T10)){
          T10<-0
        }
        if (is.na(T01)){
          T01<-0
        }
      } else {
        T_1<-0
        T1_<-0
      }
      
      #Bring T11 related terms together
      T__Hold<-merge(T1_,T_1,by=0,all=TRUE)
      row.names(T__Hold)<-T__Hold[,1]
      T__Hold<-T__Hold[,-1]
      T__Hold[is.na(T__Hold)]<-0
      
      dT1_<-sqrt(sum(T__Hold[,1]^2))
      dT_1<-sqrt(sum(T__Hold[,2]^2))
      if (dT1_>=dT_1){
        dT<-sqrt(sum((T__Hold[,1]-T__Hold[,2])^2))/dT1_
      } else if (dT1_<dT_1){
        dT<-sqrt(sum((T__Hold[,1]-T__Hold[,2])^2))/dT_1
      }
      PHold1<-data.frame(matrix(0,nrow=length(gl1),ncol=1))
      row.names(PHold1)<-gl1
      PHold2<-data.frame(matrix(0,nrow=length(gl2),ncol=1))
      row.names(PHold2)<-gl2
      for (m in 1:length(gl1)){
        nacheck<-1-sum(temp1[gl1[m],]==0)/length(temp1[gl1[m],])
        if (is.na(nacheck)){
          nacheck<-0
        }
        PHold1[gl1[m]]<-nacheck
      }
      for (m in 1:length(gl2)){
        nacheck<-1-sum(temp2[gl2[m],]==0)/length(temp2[gl2[m],])
        if (is.na(nacheck)){
          nacheck<-0
        }
        PHold2[gl2[m]]<-nacheck
      }
      PHold<-merge(PHold1,PHold2,by=0,all=TRUE)
      row.names(PHold)<-PHold[,1]
      PHold<-PHold[,-1]
      PHold[is.na(PHold)]<-0
      presence<-mean(rowMeans(PHold,na.rm = FALSE))
      if (mn==FALSE){
        KTerm<-1+presence
      } else {
        KTerm<-mn
      }
      T11<-T__Hold[,1]*T__Hold[,2]*(KTerm^(-dT))
      TRef<-row.names(T__Hold)
      for (m in 1:length(T11)){
        tanmathold[((j-1)*coln+l),TRef[m]]<-T11[m]
      }
      
      tan1[((j-1)*coln+l)]<-sum(T11)/(T10+T01-sum(T11))
    }
  }

  Output<-data.frame(matrix(c(jac1,tan1),nrow=nrow(combo)*coln,ncol=2))
  colnames(Output)<-c('Jaccard','Tanimoto')

  FinalOut<-list("Summary"=Output,"RankInfo"=tanmathold,"Boot"=combo)
  return(FinalOut)
}

RankSimilarity<-function(RankInfo,PlotTitle){
  rInfo<-apply(RankInfo,1,order)
  cInfo<-colMeans(RankInfo)
  boxplot(RankInfo[,order(cInfo)],las=2,xlab='Ranking',ylim=c(0,1),ylab='Similarity Contribution',main=PlotTitle)
  return(sort(cInfo))
}

RankComparison<-function(WithinRank1,WithinRank2,BetweenRank12,BetweenRank21){
  #B12 is ranked mean RankInfo of file2 compared to distribution from file 1
  #B21 is ranked mean RankInfo of file1 compared to distribution from file 2
  WN1<-colnames(WithinRank1)
  WN2<-colnames(WithinRank2)
  BN12<-colnames(BetweenRank12)
  BN21<-colnames(BetweenRank21)
  DF12<-data.frame(matrix(0,nrow=1,ncol=length(BN12)))
  colnames(DF12)<-BN12
  
}



#Similarity Plotting Function
SimPlot<-function(WithinM1,WithinSD1,WithinM2,WithinSD2,Between,PlotTitle,GName1,GName2){
  d1<-rnorm(1000,mean=WithinM1,sd=WithinSD1)
  d2<-rnorm(1000,mean=WithinM2,sd=WithinSD2)
  counter<-1
  while (sum(d1>1)!=0 | sum(d2>1)!=0 | counter<100){
    if (sum(d1>1)!=0){
      d1[d1>1]<-rnorm(sum(d1>1),mean=WithinM1,sd=WithinSD1)
    }
    if (sum(d2>1)!=0){
      d2[d2>1]<-rnorm(sum(d2>1),mean=WithinM2,sd=WithinSD2)
    }
    if (counter==99){
      if (sum(d1>1)!=0){
        d1[d1>1]<-1
      }
      if (sum(d2>1)!=0){
        d2[d2>1]<-1
      }
    }
    counter<-counter+1
  }
  d1[d1>1]<-1
  d2[d2>1]<-1
  boxplot(d1,d2,ylim=c(0,1),col=3,main=PlotTitle,names=c(GName1,GName2),xlab='Samples',ylab='Similarity Ratio',outline=FALSE)
  abline(h=Between,col=2,lty=5,lwd=2)
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

AdvancedSimPlot<-function(filename1,filename2,kmin=2,simtype="Tanimoto",PlotTitle,GName1,GName2,rel=TRUE,mn=FALSE,gg=FALSE,gg2=FALSE){
  #Change Within Sim to give all values, not just mean and sd
  #change Between to boot strap comparisons
  library(ggplot2)
  if (simtype=='Tanimoto'){
    WithinM1<-WithinSimBootstrap1(filename1,kmin,rel=rel,mn=mn)
    WithinM2<-WithinSimBootstrap1(filename2,kmin,rel=rel,mn=mn)
    Between1<-BetweenSimBootstrap(filename1,filename2,WithinM1$Boot,kmin,rel,mn=mn)$Summary[,2]
    Between2<-BetweenSimBootstrap(filename2,filename1,WithinM2$Boot,kmin,rel,mn=mn)$Summary[,2]
  } else if (simtype=='Jaccard'){
    WithinM1<-WithinSim(filename1,kmin,rel)[,1]
    WithinM2<-WithinSim(filename2,kmin,rel)[,1]
    Between<-BetweenSim(filename1,filename2,kmin,rel)[1]
  } else {
    print('simtype must equal Jaccard or Tanimoto')
  }
  
  d1<-WithinM1$Summary[,2]
  d2<-WithinM2$Summary[,2]
  b1<-mean(Between1)
  b1s<-format(round(sd(Between1),3),nsmall=3)[1]
  b2<-mean(Between2)
  b2s<-format(round(sd(Between2),3),nsmall=3)[1]
  
  if (is.na(b1)){
    b1<-0
    b1s<-0
  }
  if (is.na(b2)){
    b2<-0
    b2s<-0
  }
  
  if (gg){
    if (gg2){
      dframe<-data.frame(c(d1,d2,Between1,Between2),c(rep(GName1,length(d1)),rep(GName2,length(d2)),rep(paste(GName2,'in',GName1),length(Between1)),rep(paste(GName1,'in',GName2),length(Between2))))
      colnames(dframe)<-c('Similarity','Distribution')
      
    } else {
      dframe<-data.frame(c(d1,d2),c(rep(GName1,length(d1)),rep(GName2,length(d2))))
      colnames(dframe)<-c('Similarity','Distribution')
      
    }
    p<-ggplot(dframe,aes(x=Distribution,y=Similarity,color=Distribution))+geom_violin(trim=FALSE)
    p+stat_summary(fun.y=mean,geom="point",shape=23,size=2)
    #p+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    
  } else {
    boxplot(d1,d2,ylim=c(0,1),col=3,main=PlotTitle,names=c(GName1,GName2),xlab='Samples',ylab='Similarity Ratio',outline=FALSE)
    abline(h=b1,col=2,lty=5,lwd=2)
    if (b1>.1){
      text(1,b1-.05,labels=paste(GName2,'in',GName1,format(round(b1,3),nsmall=3),'+-',b1s))
    } else {
      text(1,b1+.05,labels=paste(GName2,'in',GName1,format(round(b1,3),nsmall=3),'+-',b1s))
    }
    
    abline(h=b2,col=4,lty=5,lwd=2)
    if (b2>.1){
      text(2,b2-.05,labels=paste(GName1,'in',GName2,format(round(b2,3),nsmall=3),'+-',b2s))
    } else {
      text(2,b2+.05,labels=paste(GName1,'in',GName2,format(round(b2,3),nsmall=3),'+-',b2s))
    }
    
    text(1,1,labels=paste('N=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[,1]),':',format(round(mean(d1),3),nsmall=3),'+-',format(round(sd(d1),3),nsmall=3)[1]))
    text(2,1,labels=paste('N=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[,1]),':',format(round(mean(d2),3),nsmall=3),'+-',format(round(sd(d2),3),nsmall=3)[1]))
  }
    
  
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
  failvec1<-c()
  failvec2<-c()
  if (length(UniCle2)>=length(UniCle1)){
    for (j in seq(length(UniCle2))){
      if (UniCle2[j] %in% UniCle1){
        k<-which(UniCle1 %in% UniCle2[j])
        temp1<-datfile1[which(NameVec1==k),]
        temp2<-datfile2[which(NameVec2==j),]
        AdvancedSimPlot(temp1,temp2,kmin=1,simtype = simtype,paste0(OverallName,': ',UniCle2[j]),SampleName1,SampleName2)
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
        AdvancedSimPlot(temp1,temp2,kmin=1,simtype = simtype,paste0(OverallName,': ',UniCle1[j]),SampleName1,SampleName2)
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
  
  
}

ModelComparison<-function(filename1,filename2,filename3,PlotTitle,NameList,rel=TRUE,kmin=2,mn=FALSE){
  w1<-WithinSimBootstrap1(filename1,kmin,rel,mn)
  w2<-WithinSimBootstrap1(filename2,kmin,rel,mn)
  w3<-WithinSimBootstrap1(filename3,kmin,rel,mn)
  b12<-BetweenSimBootstrap(filename1,filename2,kmin,rel,mn)
  b21<-BetweenSimBootstrap(filename2,filename1,kmin,rel,mn)
  b13<-BetweenSimBootstrap(filename1,filename3,kmin,rel,mn)
  b31<-BetweenSimBootstrap(filename3,filename1,kmin,rel,mn)
  b23<-BetweenSimBootstrap(filename2,filename3,kmin,rel,mn)
  b32<-BetweenSimBootstrap(filename3,filename2,kmin,rel,mn)
  

  boxplot(w1$Summary[,2],w2$Summary[,2],w3$Summary[,2],ylim=c(0,1),col=3,main=PlotTitle,names=NameList,xlab='Samples',outline = FALSE)
  abline(h=b12$Summary[,2],col=2,lty=5,lwd=2)
  abline(h=b13$Summary[,2],col=1,lty=5,lwd=2)
  abline(h=b23$Summary[,2],col=4,lty=5,lwd=2)
  abline(h=b21$Summary[,2],col=2,lty=5,lwd=2)
  abline(h=b31$Summary[,2],col=1,lty=5,lwd=2)
  abline(h=b32$Summary[,2],col=4,lty=5,lwd=2)
  legend('bottomright',legend=c(paste(NameList[1],'V',NameList[2],'=',format(round(b12[2],3),nsmall = 3)),paste(NameList[1],'V',NameList[3],'=',format(round(b13[2],3),nsmall=3)),paste(NameList[2],'V',NameList[3],'=',format(round(b23[2],3),nsmall=3))),col=c(2,1,4),lty=5)
  ###FIX THIS MONSTER
  #text(1,1,labels=paste('N=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[,1]),':',format(round(w1[1,3],3),nsmall=3),'+-',format(round(w1[1,4],3),nsmall=3)[1]))
  #text(2,1,labels=paste('N=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[,1]),':',format(round(w2[1,3],3),nsmall=3),'+-',format(round(w2[1,4],3),nsmall=3)[1]))
  #text(3,1,labels=paste('N=',length(SimDataClean(filename3,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename3,kmin=kmin,rel=rel)[,1]),':',format(round(w3[1,3],3),nsmall=3),'+-',format(round(w3[1,4],3),nsmall=3)[1]))
}

