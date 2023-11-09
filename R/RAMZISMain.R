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
#' @param logopt Default=TRUE
#' @param normvec Default=c('None','None')
#' @param rel_force Default=FALSE
#' @param QualityInfo Default="WeightedContributions" Alternative="Numerator"
#' @param RankingInfo Default="Numerator" Alternative="WeightedContributions"
#'
#' @return Automated RAMZIS Output Currently underdevelopment
#' @export
#'
#' @examples ###
#' #
RAMZISMain<-function(filename1,filename2,alpha=0.05,beta=0.20,conf_thresh=2,int_thresh=0.25,observedbounds=c("ZScore",2),
                     kmin=2,kmin_int=1,rel="Joint",MVCorrection=TRUE,mn=FALSE,verbose=T,logopt=TRUE,normvec=c('None','None'),
                     rel_force=FALSE,QualityInfo="WeightedContributions",RankingInfo="Numerator"){
  #load data and clean
  datalist<-SimDataCleanJoint(filename1,filename2,kmin,rel,normvector = normvec,logoption = logopt)
  file1<-datalist$DF1
  file2<-datalist$DF2
  if ((dim(file1)[1]<2) | (dim(file2)[1]<2)){
    print("There is only one glycopeptide in one dataset. Use other metrics.")
  }
  gl1<-c(row.names(file1))
  gl2<-c(row.names(file2))
  gj<-unique(c(gl1,gl2))
  mergedf<-MatrixMerge_Helper(file1,file2)
  coln1<-dim(file1)[2]
  coln2<-dim(file2)[2]
  cols1<-seq(coln1)
  cols2<-seq(coln2)
  #generate test samples
  combo1<-SAMPLER_Helper(coln1)
  combo2<-SAMPLER_Helper(coln2)
  #send to test similarity function
  TestTempSimObject<-TestSimilarityFunction(file1,combo1,gl1,file2,combo2,gl2,MVCorrection,mn)
  tanmathold<-TestTempSimObject$Numerator
  tanmatholdW<-TestTempSimObject$Contribution
  tan1<-TestTempSimObject$Similarity
  #Produce Internal Similarity data
  InternalSimObj1<-InternalSimilarity(filename=filename1,BootSet=combo1,kmin=kmin_int,rel="Within",MVCorrection=MVCorrection,mn=mn,logopt=logopt,normvec=normvec[[1]],rel_force=rel_force)
  InternalSimObj2<-InternalSimilarity(filename=filename2,BootSet=combo2,kmin=kmin_int,rel="Within",MVCorrection=MVCorrection,mn=mn,logopt=logopt,normvec=normvec[[2]],rel_force=rel_force)
  #generate null samples
  comboN<-NULLSAMPLER_Helper(coln1,coln2)
  #send to test similarity function with mergedf
  NullTempSimObject<-TestSimilarityFunction(mergedf,comboN$NDis1,gj,mergedf,comboN$NDis2,gj,MVCorrection,mn)
  tanmatholdN<-NullTempSimObject$Numerator
  tanmatholdNW<-NullTempSimObject$Contribution
  tanN<-NullTempSimObject$Similarity
  #calculate actual similarity
  ActualTempSimObject<-SimilarityCalculation_Singular(file1,file2)
  tanmatAFinal<-ActualTempSimObject$Contribution
  taniActualF<-ActualTempSimObject$Similarity

  ObsSim<-ObservedSimilarityStats(TestTempSimObject$Similarity,ActualTempSimObject$Similarity)
  if ((dim(file1)[1]>1) & (dim(file2)[1]>1)){
    RankingData<-RankingQualityFunctionV3(TestTempSimObject,NullTempSimObject,ActualTempSimObject,InternalSimObj1,InternalSimObj2,QualityInfo,RankingInfo)
  } else {
    RankingData<-NULL
  }
  InternalOverlap1<-OverlapCalculator(InternalSimObj1$InternalSimilarity,TestTempSimObject$Similarity)
  InternalOverlap2<-OverlapCalculator(InternalSimObj2$InternalSimilarity,TestTempSimObject$Similarity)
  ConfInt1<-InternalConfidenceScore(InternalSimObj1$InternalSimilarity,TestTempSimObject$Similarity,InternalOverlap1)
  ConfInt2<-InternalConfidenceScore(InternalSimObj2$InternalSimilarity,TestTempSimObject$Similarity,InternalOverlap2)
  ModalityObj<-ModalityTest(list(combo1,combo2),InternalSimObj1$InternalSimilarity,InternalSimObj2$InternalSimilarity)
  #check if actual similarity is within bounds of Test Similarity
  if (typeof(filename1)=="list"){
    filename1<-"File_1"
  }
  if (typeof(filename2)=="list"){
    filename2<-"File_2"
  }

  if (!is.numeric(observedbounds)){
    if (observedbounds[1]=="ZScore"){
      TMean<-mean(TestTempSimObject$Similarity,na.rm=T)
      Tsd<-sd(TestTempSimObject$Similarity,na.rm=T)*as.numeric(observedbounds[2])
      observedb<-c(TMean-Tsd,TMean+Tsd)
      TestBoundCheck<-ActualTempSimObject$Similarity
    } else if(observedbounds[1]=="Quartile"){
      observedb<-c(0.25,0.75)
      TestBoundCheck<-ecdf(TestTempSimObject$Similarity)(ActualTempSimObject$Similarity)
    }
  } else{
    observedb<-observedbounds
    TestBoundCheck<-ecdf(TestTempSimObject$Similarity)(ActualTempSimObject$Similarity)
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
      print(paste0(filename2, " has Insufficient Data Quality for Valid Similarity Assessment between ",filename1," and ",filename2))
    } else {
      print(paste0(filename2," has Sufficient Data Quality."))
    }
    if (!is.null(ModalityObj$Int1)){
      print(paste0(filename1," likely has outliers due to multimodality."))
    }
    if (!is.null(ModalityObj$Int2)){
      print(paste0(filename2," likely has outliers due to multimodality."))
    }
    print(paste0("Test Similarity Modeling Quality:: Observed Similarity should be between ",round(observedb[1],3)," and ",round(observedb[2],3)))
    print(paste0("Observed Similarity is ",round(ActualTempSimObject$Similarity,3)))
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
  NullTempSimObject<-TestSimilarityFunction(mergedf,ncombos$NDis1,gj,mergedf,ncombos$NDis2,gj,MVCorrection,mn)
  NullOverlap<-OverlapCalculator(NullTempSimObject$Similarity,TestTempSimObject$Similarity)

  #General Comparison assessment
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
  #Output list set up
  SimilarityOut<-list("Test"=TestTempSimObject$Similarity,"Null"=NullTempSimObject$Similarity,
                      "Internal1"=InternalSimObj1$InternalSimilarity,"Internal2"=InternalSimObj2$InternalSimilarity,
                      "Actual"=ActualTempSimObject$Similarity)
  RankingDataOut<-list("RankMatrix"=RankingData,
                       "RankData"=list("Test"=TestTempSimObject$Numerator,"Null"=NullTempSimObject$Numerator,
                                       "Internal1"=InternalSimObj1$InternalRankingInfo,"Internal2"=InternalSimObj2$InternalRankingInfo,
                                       "Actual"=ActualTempSimObject$Contribution))
  WeightedContributions<-list("Test"=TestTempSimObject$WeightedContributions,"Null"=NullTempSimObject$WeightedContributions,
                              "Internal1"=InternalSimObj1$WeightedContributions,"Internal2"=InternalSimObj2$WeightedContributions)
  BootOut<-list(combo1,combo2,comboN$NDis1,comboN$NDis2)
  QualityOut<-list("Internal1"=list("Confidence"=ConfInt1,"Overlap"=InternalOverlap1$PercOverlap,"Alpha"=InternalOverlap1$FP,"Beta"=InternalOverlap1$FN,"Modality"=ModalityObj$Int1),
                   "Internal2"=list("Confidence"=ConfInt2,"Overlap"=InternalOverlap2$PercOverlap,"Alpha"=InternalOverlap2$FP,"Beta"=InternalOverlap2$FN,"Modality"=ModalityObj$Int2),
                   "Test"=list("ZScore"=ObsSim$ZScore,"Percentile"=ObsSim$Percentile))
  GenComOut<-list("Alpha"=NullOverlap$FP,"Beta"=NullOverlap$FN)
  # Summarize Ranking Data Matrix
  if ((dim(file1)[1]>1) & (dim(file2)[1]>1)){
    RankTrunc<-RankingData[,c("ZScore","PassOverall")]
    FailCause<-rep("Passed",dim(RankTrunc)[1])
    gpRanked<-row.names(RankingData)
    for (j in 1:dim(RankingData)[1]){
      if (RankingData$PassOverall[j]==FALSE){
        f1p<-RankingData$`Overlap%_1`[j]
        f2p<-RankingData$`Overlap%_2`[j]
        if (!is.na(f1p) & !is.na(f2p)){
          if (f1p==0 & f2p==1){
            if (gpRanked[j] %in% row.names(file2)){
              FailCause[j]<-"LowQualityIn_File2"
            } else {
              FailCause[j]<-"UnseenIn_File2"
            }
          } else if (f1p==1 & f2p==0){
            if (gpRanked[j] %in% row.names(file1)){
              FailCause[j]<-"LowQualityIn_File1"
            } else {
              FailCause[j]<-"UnseenIn_File1"
            }
          } else {
            if (f1p>0.25){
              if (f2p>0.25){
                FailCause[j]<-"LowQualityIn_Both"
              } else {
                FailCause[j]<-"LowQualityIn_File1"
              }
            } else {
              FailCause[j]<-"LowQualityIn_File2"
            }
          }
        } else{
          FailCause[j]<-"BelowMinimumObs"
        }
      }
    }
  } else {
    RankTrunc<-NULL
  }
  Int1PlotData<-GGFormatterInternal(TestTempSimObject$Similarity,InternalSimObj1$InternalSimilarity,1)
  Int2PlotData<-GGFormatterInternal(TestTempSimObject$Similarity,InternalSimObj2$InternalSimilarity,2)
  NullPlotData<-GGFormatterGeneral(TestTempSimObject$Similarity,NullTempSimObject$Similarity)
  PlotData<-list("Internal1"=Int1PlotData,"Internal2"=Int2PlotData,"General"=NullPlotData)
  RankTrunc$FailCause<-FailCause
  SummaryOut<-list("QualityChecks"=QualityOut,"GeneralComparison"=GenComOut,"Rankings"=RankTrunc)
  FinalOut<-list("Summary"=SummaryOut,"Similarity"=SimilarityOut,"PlotData"=PlotData,"RankInfo"=RankingDataOut,
            'WeightedContributions'=WeightedContributions,"Boot"=BootOut)

  return(FinalOut)
}
