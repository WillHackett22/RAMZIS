#' GlycanSplitter takes a GlycReSoft formatted glycopeptide and produces a list of its components
#'
#' @param glypep Glycopeptide to be split into list
#'
#' @return list with format "GP" is input, "GP_simp" is the simplified form, "pep" is the peptide backbone. Remaining entries are all glycan compositions: Hex, HexNAc, Fuc, Neu5Ac, Sulfate
#' @export
#'
#' @examples #
GlycanSplitter<-function(glypep){
  gp<-GlycanSimplifier(glypep)
  pep<-gsub('\\{.*','',gp)
  gly<-as.numeric(strsplit(gsub('.*\\{|\\}','',gp),';')[[1]])
  return(list("GP"=glypep,"GP_simp"=gp,"pep"=pep,
              "Hex"=gly[1],"HexNAc"=gly[2],"Fuc"=gly[3],
              "Neu5Ac"=gly[4],"Sulfate"=gly[5]))
}
