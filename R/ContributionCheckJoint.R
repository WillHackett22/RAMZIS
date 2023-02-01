#' ContributionCheckJoint QC Pipeline component
#'
#' @param SimObj Test-Null similarity object
#' @param ExclusionList failures
#' @param ks KSTest indicator for separation. DEFAULT=FALSE
#'
#' @return QC data
#' @export
#'
#' @examples #
ContributionCheckJoint<-function(SimObj,ExclusionList='Default',ks=F){
  IdentList<-colnames(SimObj$RankInfoActual)
  if (ExclusionList[1]!='Default'){
    for (j in 1:length(ExclusionList)){
      IdentList<-IdentList[-which(IdentList==ExclusionList[j])]
    }
  }
  JointCheck<-data.frame(matrix(NA,nrow=length(IdentList),ncol=3),row.names = IdentList)
  colnames(JointCheck)<-c('Both','LQuartile','Separation')
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    vic<-which(SimObj$NullRankInfoFinal[l,]==0)
    JointCheck[l,2]<-(ecdf(SimObj$NullRankInfoFinal[l,])(SimObj$RankInfoActual[l])<=0.25)
    if (ks==T){
      JointCheck[l,3]<-(ks.test(SimObj$RankInfoFinal[l,],SimObj$NullRankInfoFinal[l,],alternative='g')$p.value<=0.05)
    } else if(ks==F) {
      TDis<-SimObj$RankInfoFinal[l,]
      dT<-density(TDis,from=-0.05,to=1.05,na.rm=T)
      NDis<-SimObj$NullRankInfoFinal[l,]
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

      JointCheck[l,3]<-(Overlap<=25)
    }

    JointCheck[l,1]<-(JointCheck[l,2]&JointCheck[l,3])
  }
  return(list('Raw'=JointCheck,'WinList'=IdentList[which(JointCheck[,1]==1)],'FailList'=IdentList[which(JointCheck[,1]!=1)]))
}
