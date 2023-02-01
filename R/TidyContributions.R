#' TidyContributions reformats the similarity contribution data into TIDY friendly format
#'
#' @param SimObject Similarity Object
#' @param Boot1 Sampling of dataset 1
#' @param Boot2 Sampling of dataset 2
#' @param GPID Glycopeptide IDs. Default=NULL
#' @param CompID Comparison ID. Default='X_'
#' @param IntBool Boolean indicating if an internal similarity comparison or not. Default assumes internal. Default=T
#'
#' @return TIDY formatted contribution data
#'
#' @examples
TidyContributions<-function(SimObject,Boot1,Boot2,GPID=NULL,CompID='X_',IntBool=T){
  Temp<-SimObject$Contribution
  Tcs<-dim(Temp)[2]
  Trs<-dim(Temp)[1]
  if (!is.null(GPID)){
    for (l in 1:Trs){
      if (any(GPID$GP==row.names(Temp)[l])){
        row.names(Temp)[l]<-GPID$GPID[which(GPID$GP==row.names(Temp)[l])]
      }
    }
  }
  OutC<-data.frame(matrix(0,nrow=Trs*Tcs,ncol=6))
  colnames(OutC)<-c("GPID","Contribution","Dist","Sampling1","Sampling2","Comparison")
  OutC$GPID<-rep(row.names(Temp),each=Tcs)
  OutC$Comparison<-paste0(CompID,rep(seq(1,Tcs),Trs))
  #correlate samplings, ramzis produces order [1*[1:j],[2*[1:j],...[i*[1:j]]]
  B1List<-rep(seq(1,dim(Boot1)[1]),each=dim(Boot2)[1])
  B2List<-rep(seq(1,dim(Boot2)[1]),dim(Boot1)[1])
  OutC$Sampling1<-rep(B1List,Tcs)
  OutC$Sampling2<-rep(B2List,Tcs)
  if (IntBool==T){
    OutC<-OutC[-which(OutC$Sampling1==OutC$Sampling2),]
  }
  constemp<-c()
  for (l in 1:Trs){
    constemp<-c(constemp,Temp[l,])
  }
  OutC$Contribution<-constemp
  return(OutC)
}
