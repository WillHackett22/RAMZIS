#' SymmetricalSimBootstrap compares two datasets of glycopeptide identifications
#'
#' @param filename1 first sample dataset
#' @param filename2 second sample dataset
#' @param kmin minimum number of observations Default=2
#' @param rel Boolean designating relativity Default=TRUE
#' @param MVCorrection Boolean missing value correction, if disabled it turns 0s to NAs. Default=TRUE
#' @param mn Default=FALSE. In default settings it adjusts the distance scaling by average presence. Setting it to a numeric value will use that as a constant instead. (Recommended 1:2)
#' @param logopt Boolean indicating the use of the log transform before relative scaling of abundance. Default=TRUE
#' @param normvec Optional: the normalization vectors to adjust signal abundance between samples. Default=c('None','None')
#'
#' @return Similarity Object
#' @export
#'
#' @examples #
SymmetricalSimBootstrap<-function(filename1,filename2,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,logopt=TRUE,normvec=c('None','None')){
  #load data and acquire glycopeptides
  file1<-SimDataClean(filename1,kmin,rel,normvector = normvec[1])
  file2<-SimDataClean(filename2,kmin,rel,normvector = normvec[2])
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
  TempSimObject<-SimilarityCalculation(T__Hold,T10,T01,PHold,ncomb1,ncomb2,mn)
  tanmathold<-TempSimObject$Numerator
  tanmatholdW<-TempSimObject$Contribution
  tan1<-TempSimObject$Similarity

  #calculate null
  colnt<-coln1+coln2

  if (colnt<5 & 2<colnt){
    Nsmplt1<-CombWRep(colnt,coln1)
    Nsmplt2<-CombWRep(colnt,coln2)
    combot1<-data.frame(matrix(1,nrow=Nsmplt1,ncol=coln1))
    combot1<-data.frame(gtools::combinations(colnt,coln1,repeats.allowed = TRUE))
    combot2<-data.frame(matrix(1,nrow=Nsmplt2,ncol=coln2))
    combot2<-data.frame(gtools::combinations(colnt,coln2,repeats.allowed = TRUE))
  } else {
    print('There be a number of samples here. The bootstrap generation will be slower.')
    nbootsize<-200
    combot1<-data.frame(matrix(1,nrow=nbootsize,ncol=coln1))
    combot2<-data.frame(matrix(1,nrow=nbootsize,ncol=coln2))
    for (j in 1:nbootsize){
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
  TempSimObject<-SimilarityCalculation(T__HoldN,T10N,T01N,PHoldN,ncombt1,ncombt2,mn)
  tanmatholdN<-TempSimObject$Numerator
  tanmatholdNW<-TempSimObject$Contribution
  tanN<-TempSimObject$Similarity

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

  #Output list set up
  Output<-data.frame(matrix(c(jac1,tan1),nrow=ncombs,ncol=2))
  colnames(Output)<-c('Jaccard','Tanimoto')
  TaniOut<-data.frame(matrix(tan1,nrow=ncombs[1],ncol=1))
  colnames(TaniOut)<-c('Tanimoto')
  NullOut<-data.frame(matrix(tanN,nrow=ncombts,ncol = 1))
  colnames(NullOut)<-c('NullTani')
  FinalOut<-list("Summary"=TaniOut,"Hold"=Output,"UnadjustedRankInfo"=tanmathold,"NullRankInfoFinal"=tanmatholdNW,"NullOut"=NullOut,'RankInfoFinal'=tanmatholdW,"RankInfoActual"=tanmatAFinal,"Actual"=taniActualF,"Boot"=list(combo1,combo2,combot1,combot2))
  return(FinalOut)
}
