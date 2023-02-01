#' OverlapIdxs is used by AlphaIdxs and BetaIdxs to get the indices from the overlap data for the alpha and beta regions respectively
#'
#' @param OverlapData PlaceHolder
#' @param AlpBet PlaceHolder
#'
#' @return PlaceHolder
#'
#' @examples #
OverlapIdxs<-function(OverlapData,AlpBet){
  fi<-c()
  for (j in 1:(length(OverlapData[[AlpBet]])/2)){
    jdxh<-j*2
    jdxl<-(j-1)*2+1
    fi<-c(fi,seq(OverlapData[[AlpBet]][jdxl],OverlapData[[AlpBet]][jdxh]))
  }
  return(fi)
}

#' AlphaIdxs finds the indices from the overlap data for the alpha region
#'
#' @param OverlapData PlaceHolder
#' @param AlpBet PlaceHolder
#'
#' @return PlaceHolder
#'
#' @examples  #
AlphaIdxs<-function(OverlapData){
  fi<-OverlapIdxs(OverlapData,'FPi')
  return(fi)
}

#' BetaIdxs finds the indices from the overlap data for the beta region
#'
#' @param OverlapData PlaceHolder
#' @param AlpBet PlaceHolder
#'
#' @return PlaceHolder
#'
#' @examples #
BetaIdxs<-function(OverlapData){
  fi<-OverlapIdxs(OverlapData,'FNi')
  return(fi)
}
