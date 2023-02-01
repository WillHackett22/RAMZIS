#' NullMembershipProportion is for determining the membership of samples to peaks to the Null similarity distribution
#'
#' @param SimDist Similarity Distribution
#' @param BootDis1 The first of the null samplings
#' @param BootDis2 The second of the null samplings
#' @param ModalityList The data from Modality about
#'
#' @return The degree of proportion to a peak that each sample contributes to
#' @export
#'
#' @examples #
NullMembershipProportion<-function(SimDist,BootDis1,BootDis2,ModalityList){
  BootCombinations<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=2))
  BootCombinations[,1]<-rep(seq(1,dim(BootDis1)[1]),each=dim(BootDis2)[1])
  BootCombinations[,2]<-rep(seq(1,dim(BootDis2)[1]),dim(BootDis1)[1])
  #make density
  SimDens<-density(SimDist)
  #find minima
  LocalMinimaIdx<-ModalityList$Minima
  idxL<-length(LocalMinimaIdx)+1
  LBound<-c(0,SimDens$x[LocalMinimaIdx])
  UBound<-c(SimDens$x[LocalMinimaIdx],1)
  #determine the membership of each comparison
  BootMembers<-data.frame(matrix(NA,nrow=dim(BootDis1)[1]*dim(BootDis2)[1],ncol=dim(BootDis1)[2]+dim(BootDis2)[2]))
  MemberProp<-data.frame(matrix(NA,nrow=dim(BootCombinations)[1],ncol=3))
  MemberProp[,1]<-rowSums(BootDis1[BootCombinations[,1],]<=dim(BootDis1)[2])/dim(BootDis1)[2]
  MemberProp[,2]<-rowSums(BootDis2[BootCombinations[,2],]<=dim(BootDis2)[2])/dim(BootDis2)[2]
  MemberProp[,3]<-MemberProp[,1]*(1-MemberProp[,2])
  for (j in 1:dim(BootMembers)[2]){
    BootMembers[,j]<-rowSums(BootDis1[BootCombinations[,1],]==j)+rowSums(BootDis2[BootCombinations[,2],]==j)
  }

  #find members by
  Members<-data.frame(matrix(NA,nrow=idxL,ncol=dim(BootDis1)[2]+dim(BootDis2)[2]))
  MemberPropCount<-data.frame(matrix(NA,nrow=idxL,ncol=length(unique(unlist(MemberProp)))))
  PropList<-unique(unlist(MemberProp))
  colnames(MemberPropCount)<-PropList
  MemberDegree<-data.frame(matrix(NA,nrow=idxL,ncol=ceiling(dim(BootDis1)[2]+dim(BootDis2)[2])/2))
  colnames(MemberDegree)<-paste0(seq(1,dim(MemberDegree)[2])/dim(MemberDegree)[2])
  for (j in 1:idxL){
    Members[j,]<-colSums(BootMembers[which((SimDist>=LBound[j]&SimDist<=UBound[j])),])
    Origin<-rowSums(BootMembers[which((SimDist>=LBound[j]&SimDist<=UBound[j])),1:dim(MemberDegree)[2]])/(dim(BootDis1)[2]+dim(BootDis2)[2])
    MemberDegree[j,1]<-sum(Origin<=1/dim(MemberDegree)[2])

    for (l in 1:dim(MemberPropCount)[2]){
      MemberPropCount[j,l]<-sum(MemberProp[which((SimDist>=LBound[j]&SimDist<=UBound[j])),]==PropList[l])
    }
    for (l in 2:dim(MemberDegree)[2]){
      MemberDegree[j,l]<-sum(Origin<=(l/dim(MemberDegree)[2])&Origin>((l-1)/dim(MemberDegree)[2]))/(length(Origin)*idxL)
    }
  }
  MemberProp<-t(t(Members)/colSums(Members))
  return(MemberProp)
}
