#' SymmetricalSimBootstrapV2 compares two datasets of glycopeptide identifications. It is a temporary updated version
#'
#' @param filename1 first sample dataset
#' @param filename2 second sample dataset
#' @param kmin minimum number of observations Default=2
#' @param rel Boolean designating standardization method. Default='Joint'
#' @param MVCorrection Boolean missing value correction, if disabled it turns 0s to NAs. Default=TRUE
#' @param mn Default=FALSE. In default settings it adjusts the distance scaling by average presence. Setting it to a numeric value will use that as a constant instead. (Recommended 1:2)
#' @param logopt Boolean indicating the use of the log transform before relative scaling of abundance. Default=TRUE
#' @param normvec Optional: the normalization vectors to adjust signal abundance between samples. Default=list('None','None')
#' @param rel_force Boolean to force relativization in addition to TIC normalization. Default=FALSE
#'
#' @return Similarity Object
#' @export
#'
#' @examples #
SymmetricalSimBootstrapV2<-function(filename1,filename2,kmin=2,rel='Joint',MVCorrection=TRUE,mn=FALSE,logopt=TRUE,normvec=list('None','None'),rel_force=FALSE){
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
  combo1<-SAMPLER_Helper(coln1)
  combo2<-SAMPLER_Helper(coln2)
  jac1<-rep(0,nrow(combo1)*nrow(combo2)) # jaccard holder depricated
  #send to test similarity function
  TestTempSimObject<-TestSimilarityFunction(file1,combo1,gl1,file2,combo2,gl2,MVCorrection,mn)
  tanmathold<-TestTempSimObject$Numerator
  tanmatholdW<-TestTempSimObject$WeightedContributions
  tan1<-TestTempSimObject$Similarity

  #generate null samples
  comboN<-NULLSAMPLER_Helper(coln1,coln2)
  #send to test similarity function with mergedf
  NullTempSimObject<-TestSimilarityFunction(mergedf,comboN$NDis1,gj,mergedf,comboN$NDis2,gj,MVCorrection,mn)
  tanmatholdN<-NullTempSimObject$Numerator
  tanmatholdNW<-NullTempSimObject$WeightedContributions
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
  FinalOut<-list("Summary"=TaniOut,"RankInfo"=list("Test"=tanmathold,"Null"=tanmatholdN),"NullOut"=NullOut,'WeightedContributions'=list("Test"=tanmatholdW,"Null"=tanmatholdNW),"RankInfoActual"=tanmatAFinal,"Actual"=taniActualF,"Boot"=list(combo1,combo2,comboN$NDis1,comboN$NDis2))
  return(FinalOut)
}
