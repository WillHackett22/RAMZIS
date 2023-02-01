#' ModalityInterpreter
#'
#' @param RejectList Rejection list from modality data
#' @param IntObject1 Internal similarity object 1
#' @param Mode1 Modality data of internal 1
#' @param IntObject2 Internal similarity object 2
#' @param Mode2 Modality data of internal 2
#' @param SimObject Test-Null Similarity Object
#'
#' @return interpretation of modality results
#' @export
#'
#' @examples #
ModalityInterpreter<-function(RejectList,IntObject1,Mode1,IntObject2,Mode2,SimObject=NULL){
  Int1Out<-NULL
  Int2Out<-NULL
  if (!is.null(RejectList$Int1)){
    PrimaryPeak1<-Mode1[[1]]$Peak_x[which(Mode1[[1]]$N_Members==max(Mode1[[1]]$N_Members))]
    for (j in 1:length(RejectList$Int1)){
      if (!is.null(RejectList$Int1[[j]])){
        for (l in 1:dim(RejectList$Int1[[j]])[2]){
          if ((RejectList$Int1[[j]][1,l]==PrimaryPeak1 & RejectList$Int1[[j]][2,l]<0)|(RejectList$Int1[[j]][1,l]!=PrimaryPeak1 & RejectList$Int1[[j]][2,l]>0)){
            Int1Out<-c(Int1Out,j)
          }
        }
      }
    }
  }
  if (!is.null(RejectList$Int2)){
    PrimaryPeak2<-Mode2[[1]]$Peak_x[which(Mode2[[1]]$N_Members==max(Mode2[[1]]$N_Members))]
    for (j in 1:length(RejectList$Int2)){
      if (!is.null(RejectList$Int2[[j]])){
        for (l in 1:dim(RejectList$Int2[[j]])[2]){
          if ((RejectList$Int2[[j]][1,l]==PrimaryPeak2 & RejectList$Int2[[j]][2,l]<0)|(RejectList$Int2[[j]][1,l]!=PrimaryPeak2 & RejectList$Int2[[j]][2,l]>0)){
            Int2Out<-c(Int2Out,j)
          }
        }
      }
    }
  }
  if (!is.null(Int1Out)){
    Int1Out<-unique(Int1Out)
  }
  if (!is.null(Int2Out)){
    Int2Out<-unique(Int2Out)
  }
  if (verbose){
    print('Due to over-representation in non-primary peaks, dataset 1 should see the removal of replicates:')
    print(Int1Out)
    print('Due to over-representation in non-primary peaks, dataset 2 should see the removal of replicates:')
    print(Int2Out)
  }
  return(list(Int1Out,Int2Out))

}
