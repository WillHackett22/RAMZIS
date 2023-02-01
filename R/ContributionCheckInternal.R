#' ContributionCheckInternal
#'
#' @param SimObj SimilarityObject
#' @param IntA Internal Similarity Object of group 1
#' @param IntB Internal Similarity Object of group 2
#' @param ks Use of KS test to determine separation. Default=FALSE. Highly Permissive
#'
#' @return contributions
#' @export
#'
#' @examples #
ContributionCheckInternal<-function(SimObj,IntA,IntB,ks=F){
  IdentList<-colnames(SimObj$RankInfoActual)
  InternalCheck<-data.frame(matrix(NA,nrow=length(IdentList),ncol=3),row.names = IdentList)
  colnames(InternalCheck)<-c('Both','File1','File2')
  #1 is pass 0 is fail -1 is fails to appear in file
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    if (!(l %in% row.names(IntA$InternalRankingInfo))){
      InternalCheck[l,2]=-1
    } else {
      Med<-median(IntA$InternalRankingInfo[l,])<SimObj$RankInfoActual[l]
      if (ks==T){
        KSObj<-ks.test(IntA$InternalRankingInfo[l,],SimObj$RankInfoFinal[l,],alternative='g')
        KSVal<-KSObj$p.value<0.05
        InternalCheck[l,2]<-KSVal*Med
      } else if (ks==F){
        TDis<-SimObj$RankInfoFinal[l,]
        dT<-density(TDis,from=-0.05,to=1.05,na.rm=T)
        NDis<-IntA$InternalRankingInfo[l,]
        dN<-density(NDis,from=-0.05,to=1.05,na.rm=T)
        k<-100
        df <- merge(
          as.data.frame(dN[c("x", "y")]),
          as.data.frame(dT[c("x", "y")]),
          by = "x", suffixes = c(".T", ".N")
        )
        df$comp <- as.numeric(df$y.T > df$y.N)
        df$cross <- c(NA, diff(df$comp))
        CPoint<-which(df$cross!=0)[which(max(df$y.T[which(df$cross!=0)])==df$y.T[which(df$cross!=0)])]
        XPoint<-df$x[CPoint]
        if (df$y.T[CPoint]<(10^-10)){
          Overlap<-0

        } else {
          TArea<-k*(1-ecdf(TDis)(XPoint))
          NArea<-k*ecdf(NDis)(XPoint)
          Overlap<-round(TArea+NArea,2)/(100+100-round(TArea+NArea,2))
        }

        InternalCheck[l,2]<-(Overlap<=25)*Med
      }


    }
    if (!(l %in% row.names(IntB$InternalRankingInfo))){
      InternalCheck[l,3]=-1
    } else {
      Med<-median(IntB$InternalRankingInfo[l,])<SimObj$RankInfoActual[l]
      if (ks==T){
        KSObj<-ks.test(IntB$InternalRankingInfo[l,],SimObj$RankInfoFinal[l,],alternative='g')
        KSVal<-KSObj$p.value<0.05
        InternalCheck[l,3]<-KSVal*Med
      } else if (ks==F){
        TDis<-SimObj$RankInfoFinal[l,]
        dT<-density(TDis,from=-0.05,to=1.05,na.rm=T)
        NDis<-IntB$InternalRankingInfo[l,]
        dN<-density(NDis,from=-0.05,to=1.05,na.rm=T)
        k<-100
        df <- merge(
          as.data.frame(dN[c("x", "y")]),
          as.data.frame(dT[c("x", "y")]),
          by = "x", suffixes = c(".T", ".N")
        )
        df$comp <- as.numeric(df$y.T > df$y.N)
        df$cross <- c(NA, diff(df$comp))
        CPoint<-which(df$cross!=0)[which(max(df$y.T[which(df$cross!=0)])==df$y.T[which(df$cross!=0)])]
        XPoint<-df$x[CPoint]
        if (df$y.T[CPoint]<(10^-10)){
          Overlap<-0
        } else {
          TArea<-k*(1-ecdf(TDis)(XPoint))
          NArea<-k*ecdf(NDis)(XPoint)
          Overlap<-round(TArea+NArea,2)/(100+100-round(TArea+NArea,2))
        }

        InternalCheck[l,3]<-(Overlap<=25)*Med
      }
    }
    if ((InternalCheck[l,2]==0)&(InternalCheck[l,3]==0)){
      InternalCheck[l,1]<-0
    } else if ((InternalCheck[l,2]==-1)|(InternalCheck[l,3]==-1)){
      InternalCheck[l,1]=-1
    } else {
      InternalCheck[l,1]<-1
    }
  }
  return(list('Raw'=InternalCheck,'WinList'=IdentList[which(InternalCheck[,1]==1)],'FailList'=IdentList[which(InternalCheck[,1]!=1)]))
}
