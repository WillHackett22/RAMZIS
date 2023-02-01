#' TestMembershipProportion is for determining the membership of samples to peaks to the Internal similarity distribution
#'
#' @param SimDist Similarity Distribution
#' @param BootDis1 Bootstrapping of the first sample
#' @param BootDis2 Bootstrapping of the second sample
#' @param ModalityList Output from Mo
#'
#' @return Sample Proportionship
#' @export
#'
#' @examples #
TestMembershipProportion<-function(SimDist,BootDis1,BootDis2,ModalityList){
  BootCombinations<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=2))
  BootCombinations[,1]<-rep(seq(1,dim(BootDis1)[1]),each=dim(BootDis2)[1])
  BootCombinations[,2]<-rep(seq(1,dim(BootDis2)[1]),dim(BootDis1)[1])
  #make density
  SimDens<-DENSITY_RD(SimDist)
  #find minima
  LocalMinimaIdx<-ModalityList$Minima
  idxL<-length(LocalMinimaIdx)+1
  LBound<-c(0,SimDens$x[LocalMinimaIdx])
  UBound<-c(SimDens$x[LocalMinimaIdx],1)
  #determine the membership of each comparison
  BootMembers<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=dim(BootDis1)[2]+dim(BootDis2)[2]))
  for (j in 1:dim(BootMembers)[2]){
    if (j<=dim(BootDis1)[2]){
      BootMembers[,j]<-rowSums(BootDis1[BootCombinations[,1],]==j)
    } else{
      BootMembers[,j]<-rowSums(BootDis2[BootCombinations[,1],]==(j-dim(BootDis1)[2]))
    }
  }
  MemberProp<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))


  Members<-data.frame(matrix(NA,nrow=idxL,ncol=max(BootDis1)))
  PropList<-paste(seq(1,max(BootDis1)))
  colnames(Members)<-PropList
  rownames(Members)<-paste(UBound)
  for (j in 1:idxL){
    Members[j,]<-colSums(BootMembers[which((SimDist>=LBound[j]&SimDist<=UBound[j])),])
  }
  MemberProp<-t(t(Members)/colSums(Members))
  return(MemberProp)
}
