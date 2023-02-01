#' MAIN_OUTPUT_GATHERER is used by RAMZISMain() to collect and organize all outputs of the similarity comparisons
#'
#' @param Internal1 PlaceHolder
#' @param Internal2 PlaceHolder
#' @param Test PlaceHolder
#' @param Null PlaceHolder
#' @param Observed PlaceHolder
#' @param ITO1 PlaceHolder
#' @param ITO2 PlaceHolder
#' @param NTO PlaceHolder
#' @param Boot1 PlaceHolder
#' @param Boot2 PlaceHolder
#' @param NullBoot1 PlaceHolder
#' @param NullBoot2 PlaceHolder
#'
#' @return PlaceHolder
#' @export
#'
#' @examples #
MAIN_OUTPUT_GATHERER<-function(Internal1,Internal2,Test,Null,Observed,ITO1,ITO2,NTO,Boot1,Boot2,NullBoot1,NullBoot2){
  il1<-length(Internal1$Similarity)
  il2<-length(Internal2$Similarity)
  sl<-length(Test$Similarity)
  nl<-length(Null$Similarity)
  tl<-il1+il2+sl+nl
  #Gather Similarities
  OutL<-data.frame(matrix(0,nrow=tl+1,ncol=2))
  colnames(OutL)<-c("Similarity","Distribution")
  OutL$Similarity<-c(Internal1$Similarity,Internal2$Similarity,Test$Similarity,Null$Similarity,Observed$Similarity)
  OutL$Distribution<-c(rep("Internal1",il1),rep("Internal2",il2),rep("Test",sl),rep("Null",nl),"Observed")
  #Gather Glycopeptides
  gp_all<-unique(c(row.names(Internal1$Contribution),row.names(Internal2$Contribution),row.names(Test$Contribution),row.names(Null$Contribution)))
  GPID<-data.frame(matrix(0,nrow=length(gp_all),ncol=2))
  colnames(GPID)<-c("GP","GPID")
  GPID$GP<-gp_all
  GPID$GPID<-paste0('GP_',seq(1,length(gp_all)))
  #Gather Contributions
  OutC_I1<-TidyContributions(Internal1,Boot1,Boot1,CompID ='I1_',IntBool=T)
  OutC_I2<-TidyContributions(Internal2,Boot2,Boot2,CompID ='I2_',IntBool=T)
  OutC_T<-TidyContributions(Test,Boot1,Boot2,CompID ='T_',IntBool=F)
  OutC_N<-TidyContributions(Null,NBoot1,NBoot2,CompID ='N_',IntBool=F)
  OutC<-rbind(OutC_I1,OutC_I2,OutC_T,OutC_N)
  #Gather Densities for easier plotting
  I1Dens<-DENSITY_RD(Internal1$Similarity)
  I2Dens<-DENSITY_RD(Internal2$Similarity)
  TDens<-DENSITY_RD(Test$Similarity)
  NDens<-DENSITY_RD(Null$Similarity)
  xdenslen<-length(I1Dens$x)
  OutD<-data.frame(matrix(0,nrow=xdenslen*4,ncol=4))
  colnames(OutD)<-c("Similarity","Density","Distribution")
  OutD$Similarity<-rep(I1Dens$x,4)
  OutD$Density<-c(I1Dens$y,I2Dens$y,TDens$y,NDens$y)
  OutD$Distribution<-rep(c("Internal1","Internal2","Test","Null"),each=xdenslen)
  # Store Summary Data of Overlaps
  OutO<-data.frame(matrix(0,nrow=3,ncol=4))
  colnames(OutO)<-c("Comparison","Alpha","Beta","ABTotal","Confidence","ReferenceMean","TestMean")
  TestMean<-mean(Test$Similarity,na.rm=T)
  OutO[1,]<-c("Internal1",ITO1$FP,ITO1$FN,ITO1$PercOverlap,InternalConfidenceScore(Internal1$Similarity,Test$Similarity,ITO1),mean(Internal1$Similarity,na.rm=T),TestMean)
  OutO[2,]<-c("Internal2",ITO2$FP,ITO2$FN,ITO2$PercOverlap,InternalConfidenceScore(Internal2$Similarity,Test$Similarity,ITO2),mean(Internal2$Similarity,na.rm=T),TestMean)
  OutO[3,]<-c("Null",NTO$FP,NTO$FN,NTO$PercOverlap,NA,mean(Null$Similarity,na.rm=T),TestMean)
  #Final Out List
  Out<-list("Distributions"=OutL,"Contribution"=OutC,"Densities"=OutD,"GPID"=GPID,"OverlapSummary"=OutO,"Samplings"=list("Internal1"=combo1,"Internal2"=combo2,"Null"=ncombos))
  return(Out)
}
