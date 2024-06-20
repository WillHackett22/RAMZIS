#' GlycanSimplifier reduces GlycReSoft glycan notation into abbreviated form
#'
#' @param glypep The gylcopeptide to be shortened
#'
#' @return The abbreviated glycopeptide in format PEPTIDE[Hex;HexNAc;Fuc;Neu5Ac;Sulfation]
#' @export
#'
#' @examples
#' GlycanSimplifier('AMIN(N-Glycosylation)OACIDS{Hex:7; HexNAc:6; Neu5Ac:3}')
GlycanSimplifier<-function(glypep,glysplit='\\{'){
  list<-c(0,0,0,0,0) #hex, hexnac, fuc, neu5ac, sul/phos
  targets<-c('Hex:[0-9];|Hex:[0-9][0-9];', 'HexNAc:[0-9]}|HexNAc:[0-9];', 'Fuc:[0-9];', 'Neu5Ac:[0-9]}|Neu5Ac:[0-9];', '\\@sulfate:[0-9];')
  pepn<-gsub(paste0(glysplit,'.*'),'',glypep)
  pepn<-gsub('N\\(.*\\)','n',pepn)
  pepn<-gsub('\\(.*\\)','*',pepn)
  glyn<-gsub(paste0('.*',glysplit),'',glypep)
  for (j in 1:length(targets)){
    ii<-gregexpr(targets[j],glyn)
    if (ii[[1]][1]!=-1){
      glycomp<-substr(glyn,ii[[1]][1],ii[[1]][1]+attr(ii[[1]],'match.length'))
      ij<-gregexpr(':',glycomp)[[1]][1]
      ik<-gregexpr(';|\\}',glycomp)[[1]][1]
      glyu<-substr(glycomp,ij+1,ik-1)
      list[j]<-as.numeric(glyu)
    }
  }
  return(paste0(c(pepn,'{',paste0(c(list),collapse=';'),'}'),collapse=''))
}
