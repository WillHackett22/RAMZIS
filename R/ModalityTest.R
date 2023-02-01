#' ModalityTest tests for over or under-represented membership in peaks
#'
#' @param SimObject Test-Null Similarity Object
#' @param IntObject1 Internal Similarity Object
#' @param IntObject2 Internal Similarity Object
#'
#' @return Rejection list of samples based on Modality criterion
#' @export
#'
#' @examples #
ModalityTest<-function(SimObject,IntObject1,IntObject2){
  Int1Mode<-Modality(IntObject1$InternalTanimoto)
  Int2Mode<-Modality(IntObject2$InternalTanimoto)
  IMP1<-InternalMembershipProportion(IntObject1$InternalTanimoto,SimObject$Boot[[1]],Int1Mode)
  IMP2<-InternalMembershipProportion(IntObject2$InternalTanimoto,SimObject$Boot[[2]],Int2Mode)
  TZ<-NULL
  NZ<-NULL
  IZ1<-ModalityZ(Int1Mode,IMP1)
  IZ2<-ModalityZ(Int2Mode,IMP2)
  RejectList<-list('Test'=TZ,'Null'=NZ,'Int1'=IZ1,'Int2'=IZ2)
  return(RejectList)
}
