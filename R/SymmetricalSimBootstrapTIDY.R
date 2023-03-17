#' SymmetricalSimBootstrapTIDY compares two datasets of glycopeptide identifications. It is a temporary updated version
#'
#' @param filename1 first sample dataset
#' @param filename2 second sample dataset
#' @param kmin minimum number of observations Default=2
#' @param rel Boolean designating standardization method. Default='Joint'
#' @param MVCorrection Boolean missing value correction, if disabled it turns 0s to NAs. Default=TRUE
#' @param mn Default=FALSE. In default settings it adjusts the distance scaling by average presence. Setting it to a numeric value will use that as a constant instead. (Recommended 1:2)
#' @param logopt Boolean indicating the use of the log transform before relative scaling of abundance. Default=TRUE
#' @param normvec Optional: the normalization vectors to adjust signal abundance between samples. Default=list('None','None')
#'
#' @return Similarity Object
#' @export
#'
#' @examples #
SymmetricalSimBootstrapTIDY<-function(filename1,filename2,kmin=2,rel='Joint',MVCorrection=TRUE,mn=FALSE,logopt=TRUE,normvec=list('None','None')){
  #load data and acquire glycopeptides
  datalist<-SimDataCleanJoint(filename1,filename2,kmin,rel,normvector = normvec,logoption = logopt)
  file1<-datalist$DF1
  file2<-datalist$DF2
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
  combo1<-SAMPLER_HelperTIDY(coln1,GroupID="TestA")
  combo2<-SAMPLER_HelperTIDY(coln2,GroupID="TestB")
  jac1<-rep(0,nrow(combo1)*nrow(combo2)) # jaccard holder depricated
  #send to test similarity function
  TestTempSimObject<-TestSimilarityFunctionTIDY(file1,combo1,gl1,file2,combo2,gl2,MVCorrection,mn)
  tanmat<-TestTempSimObject$Details
  tan1<-TestTempSimObject$Similarity

  #generate null samples
  comboN<-NULLSAMPLER_HelperTIDY(coln1,coln2)
  #send to test similarity function with mergedf
  NullTempSimObject<-TestSimilarityFunction(mergedf,comboN$NDis1,gj,mergedf,comboN$NDis2,gj,MVCorrection,mn)
  tanmatN<-NullTempSimObject$Details
  tanN<-NullTempSimObject$Similarity

  #calculate actual similarity
  ActualTempSimObject<-SimilarityCalculation_Singular(file1,file2)
  tanmatAFinal<-ActualTempSimObject$Contribution
  taniActualF<-ActualTempSimObject$Similarity

  #Output list set up
  ncombs<-nrow(combo1)*nrow(combo2)
  TaniOut<-data.frame(matrix(unlist(tan1$SimTot),nrow=ncombs,ncol=1))
  colnames(TaniOut)<-c('TestTanimoto')
  ncombns<-length(tanN$SimTot)
  NullOut<-data.frame(matrix(tanN$SimTot,nrow=ncombns,ncol = 1))
  colnames(NullOut)<-c('NullTanimoto')
  FinalOut<-list("Summary"=TaniOut,"TestRankInfo"=tanmathold,"NullRankInfo"=tanmatN,"NullOut"=NullOut,"RankInfoActual"=tanmatAFinal,"Actual"=taniActualF,"Boot"=list(combo1,combo2,comboN$NDis1,comboN$NDis2))
  return(FinalOut)
}
