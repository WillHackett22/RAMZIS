#' SimilarityCalculation calculates the adjusted tanimoto similarity coefficient for two samplings
#'
#' @param ABMatrix joint data set matrix
#' @param AATerm vector for the A denominator terms
#' @param BBTerm vector for the B denominator terms
#' @param PresenceMatrix presence weights for comparisons of size equal to ABMatrix
#' @param AMatSize number of A samplings
#' @param BMatSize number of B samples
#' @param DistScale distance scaling term
#'
#' @return Similarity coefficients and the composite parts
#' @export
#'
#' @examples #
SimilarityCalculation<-function(ABMatrix,AATerm,BBTerm,PresenceMatrix,AMatSize,BMatSize,DistScale){
  TotalSize<-AMatSize+BMatSize
  dTtemp<-apply(ABMatrix[,1:AMatSize],2,ManhattanVectorDistance,dataframe1=ABMatrix[,(AMatSize+1):TotalSize])
  dT<-data.frame(matrix(unlist(dTtemp),nrow=nrow(ABMatrix)),row.names = row.names(ABMatrix))
  prestemp<-apply(PresenceMatrix[,1:AMatSize],2,VectorMatrixMean,dataframe1=PresenceMatrix[,(AMatSize+1):TotalSize])
  presence<-data.frame(matrix(unlist(prestemp),nrow=nrow(PresenceMatrix)),row.names = row.names(PresenceMatrix))
  T11temp<-apply(ABMatrix[,1:AMatSize],2,VectorMatrixMultiplication,dataframe1=ABMatrix[,(AMatSize+1):TotalSize])
  T11mod<-data.frame(matrix(unlist(T11temp),nrow=nrow(ABMatrix)),row.names = row.names(ABMatrix))
  if (DistScale==FALSE){
    KTerm<-1+presence
    row.names(KTerm)<-row.names(PresenceMatrix)
  } else {
    KTerm<-DistScale
  }
  T11<-T11mod*(KTerm)^-dT
  tanmathold<-t(T11)
  Denominator<-rep(0,dim(T11)[2])
  for (j in 1:length(AATerm)){
    Denominator[((j-1)*length(BBTerm)+1):(j*length(BBTerm))]<-unlist(AATerm[j])+unlist(BBTerm)
  }
  tanmatholdW<-t(t(T11)/(Denominator-colSums(T11,na.rm=T)))
  tan1<-colSums(tanmatholdW,na.rm=T)
  Outlist<-list("Numerator"=tanmathold,"Contribution"=tanmatholdW,"Similarity"=tan1)
}
