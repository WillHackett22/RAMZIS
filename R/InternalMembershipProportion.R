#' InternalMembershipProportion is for determining the membership of samples to peaks to the Internal similarity distribution
#'
#' @param IntDist Internal Similarity Distribution
#' @param BootDis1 Sampling
#' @param ModalityList Output from Modality()
#'
#' @return Sample Proportionship
#' @export
#'
#' @examples #
InternalMembershipProportion<-function(IntDist,BootDis1,ModalityList){
  BootCombinations<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis1)[1],ncol=2))
  BootCombinations[,1]<-rep(seq(1,dim(BootDis1)[1]),each=dim(BootDis1)[1])
  BootCombinations[,2]<-rep(seq(1,dim(BootDis1)[1]),dim(BootDis1)[1])
  BootCombinations<-BootCombinations[-which(BootCombinations[,1]==BootCombinations[,2]),]

  #make density
  SimDens<-density(IntDist)
  #find minima
  LocalMinimaIdx<-ModalityList$Minima
  idxL<-length(LocalMinimaIdx)+1
  LBound<-c(0,SimDens$x[LocalMinimaIdx])
  UBound<-c(SimDens$x[LocalMinimaIdx],1)
  #determine the membership of each comparison
  BootMembers<-data.frame(matrix(NA,nrow=dim(BootCombinations)[1],ncol=dim(BootDis1)[2]))

  for (j in 1:dim(BootMembers)[2]){
    BootMembers[,j]<-rowSums(BootDis1[BootCombinations[,1],]==j)+rowSums(BootDis1[BootCombinations[,2],]==j)
  }
  MemberProp<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))

  #find members by

  Members<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))
  PropList<-paste(seq(1,max(BootDis1)))
  colnames(Members)<-PropList
  rownames(Members)<-paste(UBound)
  for (j in 1:idxL){
    Members[j,]<-colSums(BootMembers[which((IntDist>=LBound[j]&IntDist<=UBound[j])),])
  }
  MemberProp<-t(t(Members)/colSums(Members))
  return(MemberProp)
}
