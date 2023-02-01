#' RankingCriterion
#'
#' @param SimObj similarity object
#' @param ExclusionList failures
#'
#' @return Ranking and QC data
#' @export
#'
#' @examples #
RankingCriterion<-function(SimObj,ExclusionList='Default'){
  IdentList<-colnames(SimObj$RankInfoActual)
  if (ExclusionList!='Default'){
    for (j in 1:length(ExclusionList)){
      IdentList<-IdentList[-which(IdentList==ExclusionList[j])]
    }
  }
  RankInfo<-data.frame(rep(NA,length(IdentList)),row.names = IdentList)
  colnames(RankInfo)<-'Identifications'
  for (j in 1:length(IdentList)){
    l<-IdentList[j]
    RankInfo[l,1]<-mean(SimObj$RankInfoFinal[l,])*var(SimObj$RankInfoFinal[l,])*(abs(ecdf(SimObj$RankInfoFinal[l,])(SimObj$RankInfoActual[l][[1]]))+1/length(SimObj$RankInfoFinal[l,]))
  }
  return(RankInfo)
}
