#' InternalSimilarity calculates the similarity of an sample group against itself using bootstraps derived
#'
#' @param filename abundance dataset of sample group
#' @param BootSet Sampling of sample group used in Test comparisons
#' @param kmin minimum number of observations needed to be considered real. Default=1
#' @param rel Boolean indicating relative abundance transformation
#' @param MVCorrection Boolean missing value correction, if disabled it turns 0s to NAs. Default=TRUE
#' @param mn Default=FALSE. In default settings it adjusts the distance scaling by average presence. Setting it to a numeric value will use that as a constant instead. (Recommended 1:2)
#' @param logopt Boolean indicating the use of the log transform before relative scaling of abundance. Default=TRUE
#' @param normvec Optional: the normalization vectors to adjust signal abundance between samples. Default=c('None')
#'
#' @return Internal Similarity Object
#' @export
#'
#' @examples #
InternalSimilarity<-function(filename,BootSet,kmin=1,rel=TRUE,MVCorrection=TRUE,mn=FALSE,logopt=TRUE,normvec=c('None')){
  file1<-SimDataCleanLogA(filename,kmin,rel,normvector = normvec,logoption=logopt)
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
  TempSimObject<-SimilarityCalculation(T__Hold,T10,T01,PHold,ncomb1,ncomb1,mn)
  tanmatholdW<-TempSimObject$Contribution
  tan1<-TempSimObject$Similarity
  keepidx<-c()
  for (j in 1:ncomb1){
    keepidx<-c(keepidx,(ncomb1*(j-1)+j):(ncomb1*j))
    ident<-ncomb1*(j-1)+j
    keepidx<-keepidx[-which(ident==keepidx)]
  }
  tanmatholdW<-tanmatholdW[,keepidx]
  tan1<-tan1[keepidx]


  OutputObj<-list("InternalTanimoto"=tan1,"InternalRankingInfo"=tanmatholdW)
  return(OutputObj)
}
