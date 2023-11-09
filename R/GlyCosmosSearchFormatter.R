#' GlyCosomosSearchFormatter takes a vector of of GlycReSoft formatted glycopeptides and turns it into a .tsv for input into glycosmos composition search engine at https://glycosmos.org/glycans/composition
#'
#' @param glycopeptides vector of glycopeptides
#' @param outfile name and location of file to be written to. adds .tsv to given input
#' @param verbose Default False. Set to true to get table in console
#'
#' @return .tsv file for glycosmos composition search
#' @export
#' @import utils
#'
#' @examples #
GlyCosmosSearchFormatter<-function(glycopeptides,outfile,verbose=F){
  decomp<-as.data.frame(do.call(rbind,lapply(glycopeptides,GlycanSplitter)))
  out<-data.frame(matrix(0,nrow=length(glycopeptides,ncol=8)))
  colnames(out)<-c("hex","hexnac","dhex","neu5ac","neu5gc","P","S","Ac")
  out$hex<-unlist(decomp$Hex)
  out$hexnac<-unlist(decomp$HexNAc)
  out$dhex<-unlist(decomp$Fuc)
  out$neu5ac<-unlist(decomp$Neu5Ac)
  out$S<-unlist(decomp$Sulfate)
  utils::write.table(out,paste0(outfile,'.tsv'),row.names=F,sep="\t")
  if (verbose==T){
    return(out)
  }
}
