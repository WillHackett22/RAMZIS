#' ContributionCheckComp part of the QC set of functions
#'
#' @param SimObj similarity object
#' @param ExclusionList list of already failed
#'
#' @return More QC data
#' @export
#'
#' @examples #
ContributionCheckComp<-function(SimObj,ExclusionList='Default'){
  IdentList<-colnames(SimObj$RankInfoActual)
  if (ExclusionList[1]!='Default'){
    for (j in 1:length(ExclusionList)){
      IdentList<-IdentList[-which(IdentList==ExclusionList[j])]
    }
  }
  CompCheck<-data.frame(matrix(NA,nrow=length(IdentList),ncol=3),row.names = IdentList)
  colnames(CompCheck)<-c('Both','Centrality','NonZero')
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    center<-ecdf(SimObj$RankInfoFinal[l,])(SimObj$RankInfoActual[l])
    nzed<-1-sum(SimObj$RankInfoFinal[l,]==0)/length(SimObj$RankInfoFinal[l,])
    CompCheck[l,2]<-((center<=0.75)&(center>=0.25))
    CompCheck[l,3]<-(nzed>=0.85)
    CompCheck[l,1]<-(CompCheck[l,2]&CompCheck[l,3])
  }
  return(list('Raw'=CompCheck,'WinList'=IdentList[which(CompCheck[,1]==T)],'FailList'=IdentList[which(CompCheck[,1]!=T)]))
}
