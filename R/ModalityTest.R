#' ModalityTest tests for over or under-represented membership in peaks
#'
#' @param BootSet List of samplings produced for the Test and Internal comparisons
#' @param IntObject1 Internal Similarity
#' @param IntObject2 Internal Similarity
#'
#' @return Rejection list of samples based on Modality criterion
#' @export
#'
#' @examples #
ModalityTest<-function(BootSet,IntSim1,IntSim2){
  Int1Mode<-Modality(IntSim1)
  Int2Mode<-Modality(IntSim2)
  IMP1<-InternalMembershipProportion(IntSim1,BootSet[[1]],Int1Mode)
  IMP2<-InternalMembershipProportion(IntSim2,BootSet[[2]],Int2Mode)
  IZ1<-ModalityZ(Int1Mode,IMP1)
  IZ2<-ModalityZ(Int2Mode,IMP2)
  RejectList<-list('Int1'=IZ1,'Int2'=IZ2)
  return(RejectList)
}
