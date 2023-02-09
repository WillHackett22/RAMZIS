SymmetricalSimBootstrapV2<-function(filename1,filename2,kmin=2,rel=TRUE,MVCorrection=TRUE,mn=FALSE,logopt=TRUE,normvec=list('None','None')){
  #load data and acquire glycopeptides
  file1<-SimDataClean(filename1,kmin,rel,normvector = normvec[[1]],logoption = logopt)
  file2<-SimDataClean(filename2,kmin,rel,normvector = normvec[[2]],logoption = logopt)
  gl1<-c(row.names(file1))
  gl2<-c(row.names(file2))
  gj<-unique(c(gl1,gl2))
  mergedf<-MatrixMerge_Helper(file1,file2)
  coln1<-dim(file1)[2]
  coln2<-dim(file2)[2]
  cols1<-seq(coln1)
  cols2<-seq(coln2)
  #generate test samples
  Nsmpl2<-CombWRep(coln2,coln2)
  Nsmpl1<-CombWRep(coln1,coln1)
  combo1<-SAMPLER_Helper(coln1)
  combo2<-SAMPLER_Helper(coln2)
  jac1<-rep(0,nrow(combo1)*nrow(combo2)) # jaccard holder depricated
  #send to test similarity function
  TestTempSimObject<-TestSimilarityFunction(file1,combo1,gl1,file2,combo2,gl2,MVCorrection,mn)
  tanmathold<-TestTempSimObject$Numerator
  tanmatholdW<-TestTempSimObject$Contribution
  tan1<-TestTempSimObject$Similarity

  #generate null samples
  comboN<-NULLSAMPLER_Helper(coln1,coln2)
  #send to test similarity function with mergedf
  NullTempSimObject<-TestSimilarityFunction(mergedf,comboN$NDis1,gj,mergedf,comboN$NDis2,gj,MVCorrection,mn)
  tanmatholdN<-NullTempSimObject$Numerator
  tanmatholdNW<-NullTempSimObject$Contribution
  tanN<-NullTempSimObject$Similarity
  jacN<-rep(0,nrow(comboN$NDis1)*nrow(comboN$NDis2)) # jaccard holder depricated

  #calculate actual similarity
  ActualTempSimObject<-SimilarityCalculation_Singular(file1,file2)
  tanmatAFinal<-ActualTempSimObject$Contribution
  taniActualF<-ActualTempSimObject$Similarity

  #Output list set up
  ncombs<-nrow(combo1)*nrow(combo2)
  TaniOut<-data.frame(matrix(tan1,nrow=ncombs,ncol=1))
  colnames(TaniOut)<-c('Tanimoto')
  ncombns<-length(tanN)
  NullOut<-data.frame(matrix(tanN,nrow=ncombns,ncol = 1))
  colnames(NullOut)<-c('NullTani')
  FinalOut<-list("Summary"=TaniOut,"UnadjustedRankInfo"=tanmathold,"NullRankInfoFinal"=tanmatholdNW,"NullOut"=NullOut,'RankInfoFinal'=tanmatholdW,"RankInfoActual"=tanmatAFinal,"Actual"=taniActualF,"Boot"=list(combo1,combo2,comboN$NDis1,comboN$NDis2))
  return(FinalOut)
}
