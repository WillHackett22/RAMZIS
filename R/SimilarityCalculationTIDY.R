#' SimilarityCalculationTIDY calculates the adjusted tanimoto similarity coefficient for two samplings
#'
#' @param ABMatrix joint data set matrix
#' @param AATerm vector for the A denominator terms
#' @param BBTerm vector for the B denominator terms
#' @param PresenceMatrix presence weights for comparisons of size equal to ABMatrix
#' @param AMatSize number of A samplings
#' @param BMatSize number of B samplings
#' @param DistScale distance scaling term
#'
#' @return Similarity coefficients and the composite parts
#' @export
#'
#' @examples #
SimilarityCalculationTIDY<-function(A_Obj,B_Obj,AMatSize,BMatSize,IDVec,DistScale){
  CombSize<-AMatSize*BMatSize
  ABtemp<-lapply(IDVec,PerIdent,dfA=A_Obj$Partial,dfB=B_Obj$Partial)
  ABMatrix<-do.call("rbind",ABtemp)
  if (DistScale==FALSE){
    ABMatrix$KTerm<-1+ABMatrix$Presence
  } else {
    ABMatrix$KTerm<-DistScale
  }
  ABMatrix$T11<-ABMatrix$AB*(ABMatrix$KTerm)^-ABMatrix$dAB
  AABB<-data.frame(matrix(0,nrow=CombSize,ncol=3))
  colnames(AABB)<-c("SampleIDA","SampleIDB","AABB")
  AABB$SampleIDA<-rep(A_Obj$Self$SampleID,each=BMatSize)
  AABB$SampleIDB<-rep(B_Obj$Self$SampleID,AMatSize)
  AABB$AABB<-rep(A_Obj$Self$Value,each=BMatSize)+rep(B_Obj$Self$Value,AMatSize)
  ABK<-ABMatrix %>%
    group_by(SampleIDA,SampleIDB) %>%
    summarise(ABk=sum(T11))
  TermMatrix<-merge(AABB,ABK)
  TermMatrix$Denom<-TermMatrix$AABB+TermMatrix$ABk
  ABMatrix<-merge(ABMatrix,TermMatrix)
  ABMatrix$Sim<-ABMatrix$T11/ABMatrix$Denom
  tan1<-ABMatrix %>% group_by(SampleIDA,SampleIDB) %>%
    summarise(SimTotal=sum(Sim))
  Outlist<-list("Details"=ABMatrix,"Similarity"=tan1)
  return(Outlist)
}


PerIdent<-function(gi,dfA,dfB){
  ixAg<-which(dfA$Identification==gi)
  ixBg<-which(dfB$Identification==gi)
  ila<-length(ixAg)
  ilb<-length(ixBg)
  out<-data.frame(matrix(0,nrow=ila*ilb,ncol=6))
  colnames(out)<-c("Identification","AB","Presence","dAB","SampleIDA","SampleIDB")
  out$Identification<-rep(gi,ila*ilb)
  vec1<-rep(dfA$Value[ixAg],each=ilb)
  vec2<-rep(dfB$Value[ixBg],ila)
  out$dAB<-abs(vec1-vec2)
  out$AB<-vec1*vec2
  out$Presence<-(rep(dfA$Presence[ixAg],each=ilb)+rep(dfB$Presence[ixBg],ila))/2
  out$SampleIDA<-rep(dfA$SampleID[ixAg],each=ilb)
  out$SampleIDB<-rep(dfB$SampleID[ixBg],ila)
  return(out)
}
