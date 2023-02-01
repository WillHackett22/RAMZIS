#' RAMZISMain is intended to be the hassle free way to run RAMZIS. It should produce all the similarity comparisons necessary for a given set of files
#'
#' @param filename1 The first file or dataset of the similarity comparison
#' @param filename2 The second file or dataset of the similarity comparison
#' @param alpha Default=0.05
#' @param beta Default=0.20
#' @param conf_thresh Default=2
#' @param int_thresh Default=0.25
#' @param observedbounds Default= c("ZScore",3)
#' @param kmin Default=2
#' @param kmin_int Default=1
#' @param rel Default=TRUE
#' @param MVCorrection default=TRUE
#' @param mn Default=FALSE
#' @param verbose default=TRUE
#' @param logopt Default=FALSE
#' @param normvec Default=c('None','None')
#'
#' @return
#' @export
#'
#' @examples
RAMZISMain<-function(filename1,filename2,alpha=0.05,beta=0.20,conf_thresh=2,int_thresh=0.25,observedbounds=c("ZScore",3),
                     kmin=2,kmin_int=1,rel=TRUE,MVCorrection=TRUE,mn=FALSE,verbose=T,logopt=FALSE,normvec=c('None','None')){
  #load data and clean
  df1<-SimDataClean(filename1,kmin,rel,normvector = normvec[1])
  df2<-SimDataClean(filename2,kmin,rel,normvector = normvec[2])
  mergedf<-MatrixMerge_Helper(df1,df2)
  coln1<-dim(df1)[2]
  coln2<-dim(df2)[2]
  cols1<-seq(coln1)
  cols2<-seq(coln2)


  #acquire glycopeptides
  gl1<-c(row.names(df1))
  gl2<-c(row.names(df2))
  glycojoint<-unique(c(gl1,gl2))
  #sample for Internal and Test
  combo1<-SAMPLER_Helper(coln1,cols1,filename1)
  combo2<-SAMPLER_Helper(coln2,cols2,filename2)
  #Generate Internal Similarity
  IntSimObj1<-InternalSimilarityFunction(filename1,combo1,kmin=kmin_int,rel=rel,MVCorrection=MVCorrection,mn=mn,logopt=logopt,normvec=normvec[1])
  IntSimObj2<-InternalSimilarityFunction(filename2,combo2,kmin=kmin_int,rel=rel,MVCorrection=MVCorrection,mn=mn,logopt=logopt,normvec=normvec[2])
  #Generate Test Similarity
  TestSimObject<-TestSimilarityFunction(df1,combo1,gl1,df2,combo2,gl2,MVCorrection,mn)
  InternalOverlap1<-OverlapCalculator(IntSimObj1$Similarity,TestSimObject$Similarity)
  InternalOverlap2<-OverlapCalculator(IntSimObj2$Similarity,TestSimObject$Similarity)
  ConfInt1<-InternalConfidenceScore(IntSimObj1$Similarity,TestSimObject$Similarity,InternalOverlap1)
  ConfInt2<-InternalConfidenceScore(IntSimObj2$Similarity,TestSimObject$Similarity,InternalOverlap2)

  #calculate actual similarity
  Actual_Temp<-SimilarityCalculation_Singular(df1,df2,mn)
  #check if actual similarity is within bounds of Test Similarity
  if (!is.numeric(observedbounds)){
    if (observedbounds[1]=="ZScore"){
      TMean<-mean(TestSimObject$Similarity,na.rm=T)
      Tsd<-sd(TestSimObject$Similarity,na.rm=T)*as.numeric(observedbounds[2])
      observedb<-c(TMean-Tsd,TMean+Tsd)
      TestBoundCheck<-Actual_Temp$Similarity
    } else if(observedbounds[1]=="Quartile"){
      observedb<-c(0.25,0.75)
      TestBoundCheck<-ecdf(TestSimObject$Similarity)(Actual_Temp$Similarity)
    }
  } else{
    observedb<-observedbounds
    TestBoundCheck<-ecdf(TestSimObject$Similarity)(Actual_Temp$Similarity)
  }
  if (verbose==T){
    print(paste0("Data Quality Thresholds:: Confidence<=",conf_thresh,"   OverlapTotal<",int_thresh))
    print(paste0(filename1," has an Internal Confidence of ",ConfInt1," and total alpha+beta=",InternalOverlap1$PercOverlap))
    if (InternalOverlap1$PercOverlap>int_thresh | ConfInt1<=conf_thresh){
      print(paste0(filename1, "has Insufficient Data Quality for Valid Similarity Assessment between ",filename1," and ",filename2))
    } else {
      print(paste0(filename1," has Sufficient Data Quality."))
    }
    print(paste0(filename2," has an Internal Confidence of ",ConfInt2," and total alpha+beta=",InternalOverlap2$PercOverlap))
    if (InternalOverlap2$PercOverlap>int_thresh | ConfInt2<=conf_thresh){
      print(paste0(filename2, "has Insufficient Data Quality for Valid Similarity Assessment between ",filename1," and ",filename2))
    } else {
      print(paste0(filename2," has Sufficient Data Quality."))
    }
    print(paste0("Test Similarity Modeling Quality:: Observed Similarity should be between ",round(observedb[1],3)," and ",round(observedb[2],3)))
    if (TestBoundCheck<observedb[1] | TestBoundCheck>observedb[2]){
      print(paste0("RAMZIS Simulation does not acceptably model comparison of ",filename1," and ",filename2,"."))
      print("Consistent failures may indicate the existence of subgroups.")
    } else {
      print(paste0("RAMZIS Simulation acceptably modeled the comparison of ",filename1," and ",filename2,"."))
    }
  }
  #Quality Checks Done
  #Begin Null Similarity Generation for comparison
  #sample null
  ncombos<-NULLSAMPLER_Helper(coln1,coln2,sampn = 200)
  ncn1dim<-dim(ncombos$NDis1)
  ncn2dim<-dim(ncombos$NDis2)
  ncrr<-ncn1dim[1]*ncn2dim[1]
  #initialize output objects
  jacN<-rep(0,ncrr) # jaccard holder depricated
  NullSimObject<-NullSimilarityFunction(df1,df2,ncombos,MVCorrection,mn)
  NullOverlap<-OverlapCalculator(NullSimObject$Similarity,TestSimObject$Similarity)
  if (verbose==T){
    print("Comparison of Test and Null Similarity Distributions Finds:")
    NTAlpha<-round(NullOverlap$FP,3)
    NTBeta<-round(NullOverlap$FN,3)
    print(paste0("Alpha=",NTAlpha))
    print(paste0("Beta=",NTBeta))
    if (NullOverlap$FP<=alpha & NullOverlap$FN<=beta){
      print(paste0("Significant Separation of Test and Null Distributions observed with alpha<=",alpha," and beta<=",beta,"in ",filename1," and ",filename2,"."))
    } else if (NullOverlap$FP>alpha & NullOverlap$FN>beta){
      print(paste0("Observed alpha>",alpha," and beta>",beta,"in ",filename1," and ",filename2,"."))
    }
  }

  #Output
  # SimilarityTotals<-list("InternalSimilarity1"=IntSimObj1$Similarity,"InternalSimilarity2"=IntSimObj2$Similarity,"ActualSimilarity"=Actual_Temp$Similarity,"TestSimilarity"=TestSimObject$Similarity,"NullSimilarity"=NullSimObject$Similarity)
  # SimilarityContributions<-list("InternalContributions1"=IntSimObj1$InternalRankingInfo,"InternalContributions2"=IntSimObj2$InternalRankingInfo,"ActualContributions"=Actual_Temp$SimilarityContributions,"TestContributions"=TestSimObject$Contribution,"NullContributions"=NullSimObject$Contribution)
  # FinalOut<-list("SimValues"=SimilarityTotals,"SimContributions"=SimilarityContributions,"Bootstraps"=list(combo1,combo2,ncombos))
  FinalOut<-MAIN_OUTPUT_GATHERER(IntSimObj1,IntSimObj2,TestSimObject,NullSimObject,Actual_Temp,InternalOverlap1,InternalOverlap2,NullOverlap,combo1,combo2,ncombos)
  return(FinalOut)
}
