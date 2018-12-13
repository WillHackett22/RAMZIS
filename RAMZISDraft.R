#JacTanIterFunction Script
#Will Hackett
#9-6-18
#Functions for Within and Between Similarity Measures of Glycopeptides
#Expected input: CSV files that have columns of total signal values from glycresoft output
#Expected input Cont: Files should be glycopeptides for a single protein
#Expected input cont: eg File1 and File2 should only have values from AGP1 to compare AGP1
#Expected input cont: The first column should be the glycopeptide, all others are signal columns from different samples
#
#Functions:
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

WithinSimCascade<-function(filename,kmin=2,rel=TRUE,MVCorrection=TRUE){
  file1<-SimDataClean(filename,kmin,rel)
  #acquire all glycopeptides in file
  glycopep1<-c(row.names(file1))
  #set up possible combinations
  coln<-ncol(file1)
  cols<-seq(coln)
  #figure out total possible number of combos
  comblt<-0
  combl<-rep(0,floor(coln/2))
  for (k in 1:floor(coln/2)){
    comblt<-comblt+length(combn(cols,k)[1,])
    combl[k]<-ncol(combn(cols,k))
  }
  #run through all combos
  jac1<-rep(0,comblt)
  tan1<-rep(0,comblt)
  for (k in 1:floor(coln/2)){
    combos<-combn(cols,k)
    #remove redundant combinations
    combo<-combos
    if (k>1){
      priorc<-sum(combl[1:(k-1)])
    } else {
      priorc<-0
    }
    
    #iterate through combinations
    for (j in 1:ncol(combo)){
      #separate datasets for combinations
      temp1<-data.frame(file1[,combo[,j]])
      row.names(temp1)<-glycopep1
      temp2<-data.frame(file1[,-combo[,j]])
      row.names(temp2)<-glycopep1
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
      jac1[priorc+j]<-M11/length(glycopep1)
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
      T11<-sum(T1_*T_1)
      tan1[priorc+j]<-T11/(T10+T01-T11)
      dT1_<-sqrt(sum(T1_^2))
      dT_1<-sqrt(sum(T_1^2))
      tan1[priorc+j]<-T11/(T10+T01-T11)
      if (dT1_>=dT_1){
        dT<-sqrt(sum((T_1-T1_)^2))/dT1_
      } else if (dT1_<dT_1){
        dT<-sqrt(sum((T_1-T1_)^2))/dT_1
      }
      pres1<-GPPresence(temp1,kmin=1,rel=rel)
      pres2<-GPPresence(temp2,kmin=1,rel=rel)
      presm<-mean(c(pres1,pres2))
      if (tan1[j]!=0){
        tan1[priorc+j]<-tan1[priorc+j]/((1+presm)^dT)
      }
      
    }
  }
  Output<-data.frame(matrix(c(mean(jac1),sd(jac1),mean(tan1),sd(tan1)),nrow=1,ncol=4))
  colnames(Output)<-c('Jaccard Mean','Jaccard SD','Tanimoto Mean','Tanimoto SD')
  return(Output)
}

BetweenSim<-function(filename1,filename2,kmin=2,rel=TRUE,MVCorrection=TRUE){
  #load data and acquire glycopeptides
  file1<-SimDataClean(filename1,kmin,rel)
  file2<-SimDataClean(filename2,kmin,rel)
  glycopep1<-c(row.names(file1))
  glycopep2<-c(row.names(file2))
  glycojoint<-unique(c(glycopep1,glycopep2))
  mergedf<-merge(file1,file2,by=0,all=FALSE)
  row.names(mergedf)<-mergedf[,1]
  mergedf<-mergedf[,-1]
  jacbet<-dim(mergedf)[1]/length(glycojoint)
  if (MVCorrection!=TRUE){
    file1[file1==0]<-NA
    file2[file2==0]<-NA
  }
  if (ncol(file1)==1){
    T1_<-file1[which(glycopep1 %in% glycopep2),1]
    T10<-sum(file1[,1]^2,na.rm = TRUE)
  } else {
    T1_<-rowMeans(file1[which(glycopep1 %in% glycopep2),],na.rm =TRUE)
    T10<-sum(rowMeans(file1,na.rm =TRUE)^2)
  }
  if (ncol(file2)==1){
    T_1<-file2[which(glycopep2 %in% glycopep1),1]
    T01<-sum(file2[,1]^2,na.rm = TRUE)
  } else {
    T_1<-rowMeans(file2[which(glycopep2 %in% glycopep1),],na.rm =TRUE)
    T01<-sum(rowMeans(file2,na.rm =TRUE)^2)
  }
  T11<-sum(T1_*T_1)
  tanbet<-T11/(T10+T01-T11)
  dT1_<-sqrt(sum(T1_^2))
  dT_1<-sqrt(sum(T_1^2))
  if (dT1_>=dT_1){
    dtan<-sqrt(sum((T_1-T1_)^2))/dT1_
  } else if (dT1_<dT_1){
    dtan<-sqrt(sum((T_1-T1_)^2))/dT_1
  }
  pres1<-GPPresence(file1,kmin=1,rel=rel)
  pres2<-GPPresence(file2,kmin=1,rel=rel)
  presm<-mean(c(pres1,pres2))
  if (tanbet!=0){
    tanbet<-tanbet/((1+presm)^dtan)
  }
  return(c(jacbet,tanbet))
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

AdvancedSimPlot<-function(filename1,filename2,kmin=2,simtype="Tanimoto",PlotTitle,GName1,GName2,rel=TRUE){
  #Change Within Sim to give all values, not just mean and sd
  #change Between to boot strap comparisons
  if (simtype=='Tanimoto'){
    WithinM1<-WithinSim(filename1,kmin,rel=rel)[,2]
    WithinM2<-WithinSim(filename2,kmin,rel=rel)[,2]
    Between<-BetweenSim(filename1,filename2,kmin,rel)[2]
  } else if (simtype=='Jaccard'){
    WithinM1<-WithinSim(filename1,kmin,rel)[,1]
    WithinM2<-WithinSim(filename2,kmin,rel)[,1]
    Between<-BetweenSim(filename1,filename2,kmin,rel)[1]
  } else {
    print('simtype must equal Jaccard or Tanimoto')
  }
  
  d1<-WithinM1
  d2<-WithinM2
  # counter<-1
  # while (sum(d1>1)!=0 | sum(d2>1)!=0 | counter<100){
  #   if (sum(d1>1)!=0){
  #     d1[d1>1]<-rnorm(sum(d1>1),mean=WithinM1,sd=WithinSD1)
  #   }
  #   if (sum(d2>1)!=0){
  #     d2[d2>1]<-rnorm(sum(d2>1),mean=WithinM2,sd=WithinSD2)
  #   }
  #   if (counter==99){
  #     if (sum(d1>1)!=0){
  #       d1[d1>1]<-1
  #     }
  #     if (sum(d2>1)!=0){
  #       d2[d2>1]<-1
  #     }
  #   }
  #   counter<-counter+1
  # }
  # counter<-1
  # while (sum(d1<0)!=0 | sum(d2<0)!=0 | counter<100){
  #   if (sum(d1<0)!=0){
  #     d1[d1<0]<-rnorm(sum(d1<0),mean=WithinM1,sd=WithinSD1)
  #   }
  #   if (sum(d2<0)!=0){
  #     d2[d2<0]<-rnorm(sum(d2<0),mean=WithinM2,sd=WithinSD2)
  #   }
  #   if (counter==99){
  #     if (sum(d1<0)!=0){
  #       d1[d1<0]<-0
  #     }
  #     if (sum(d2<0)!=0){
  #       d2[d2<0]<-0
  #     }
  #   }
  #   counter<-counter+1
  # }
  # d1[d1>1]<-1
  # d2[d2>1]<-1
  # d1[d1<0]<-0
  # d2[d2<0]<-0
  boxplot(d1,d2,ylim=c(0,1),col=3,main=PlotTitle,names=c(GName1,GName2),xlab='Samples',ylab='Similarity Ratio',outline=FALSE)
  abline(h=Between,col=2,lty=5,lwd=2)
  text(1.5,Between,labels=paste(Between))
  text(1,1,labels=paste('N=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[,1]),':',format(round(mean(WithinM1),3),nsmall=3),'+-',format(round(sd(WithinM1),3),nsmall=3)[1]))
  text(2,1,labels=paste('N=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[,1]),':',format(round(mean(WithinM2),3),nsmall=3),'+-',format(round(sd(WithinM2),3),nsmall=3)[1]))
  
}

#Function to separate glycopeptides by same peptide backbone
PeptideSegment<-function(filename1,filename2,kmin=2,simtype="Tanimoto",rel=TRUE,OverallName,Sample1,Sample2){
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
        AdvancedSimPlot(temp1,temp2,kmin=1,simtype = simtype,paste0(OverallName,': ',UniCle2[j]),Sample1,Sample2)
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
        AdvancedSimPlot(temp1,temp2,kmin=1,simtype = simtype,paste0(OverallName,': ',UniCle1[j]),Sample1,Sample2)
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

ModelComparison<-function(filename1,filename2,filename3,PlotTitle,NameList,rel=TRUE,kmin=2){
  w1<-WithinSimCascade(filename1,kmin,rel)
  w2<-WithinSimCascade(filename2,kmin,rel)
  w3<-WithinSimCascade(filename3,kmin,rel)
  b12<-BetweenSim(filename1,filename2,kmin,rel)
  b13<-BetweenSim(filename1,filename3,kmin,rel)
  b23<-BetweenSim(filename2,filename3,kmin,rel)
  d1<-rnorm(1000,mean=w1[1,3],sd=w1[1,4])
  d2<-rnorm(1000,mean=w2[1,3],sd=w2[1,4])
  d3<-rnorm(1000,mean=w3[1,3],sd=w3[1,4])
  counter<-1
  while (sum(d1>1)!=0 | sum(d2>1)!=0 | sum(d3>1)!=0  | counter<100){
    if (sum(d1>1)!=0){
      d1[d1>1]<-rnorm(sum(d1>1),mean=w1[1,3],sd=w1[1,4])
    }
    if (sum(d2>1)!=0){
      d2[d2>1]<-rnorm(sum(d2>1),mean=w2[1,3],sd=w2[1,4])
    }
    if (sum(d3>1)!=0){
      d3[d3>1]<-rnorm(sum(d3>1),mean=w3[1,3],sd=w3[1,4])
    }
    if (counter==99){
      if (sum(d1>1)!=0){
        d1[d1>1]<-1
      }
      if (sum(d2>1)!=0){
        d2[d2>1]<-1
      }
      if (sum(d3>1)!=0){
        d3[d3>1]<-1
      }
    }
    counter<-counter+1
  }
  boxplot(d1,d2,d3,ylim=c(0,1),col=3,main=PlotTitle,names=NameList,xlab='Samples',outline = FALSE)
  abline(h=b12[2],col=2,lty=5,lwd=2)
  abline(h=b13[2],col=1,lty=5,lwd=2)
  abline(h=b23[2],col=4,lty=5,lwd=2)
  legend('bottomright',legend=c(paste(NameList[1],'V',NameList[2],'=',format(round(b12[2],3),nsmall = 3)),paste(NameList[1],'V',NameList[3],'=',format(round(b13[2],3),nsmall=3)),paste(NameList[2],'V',NameList[3],'=',format(round(b23[2],3),nsmall=3))),col=c(2,1,4),lty=5)
  text(1,1,labels=paste('N=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename1,kmin=kmin,rel=rel)[,1]),':',format(round(w1[1,3],3),nsmall=3),'+-',format(round(w1[1,4],3),nsmall=3)[1]))
  text(2,1,labels=paste('N=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename2,kmin=kmin,rel=rel)[,1]),':',format(round(w2[1,3],3),nsmall=3),'+-',format(round(w2[1,4],3),nsmall=3)[1]))
  text(3,1,labels=paste('N=',length(SimDataClean(filename3,kmin=kmin,rel=rel)[1,]),', G=',length(SimDataClean(filename3,kmin=kmin,rel=rel)[,1]),':',format(round(w3[1,3],3),nsmall=3),'+-',format(round(w3[1,4],3),nsmall=3)[1]))
}

