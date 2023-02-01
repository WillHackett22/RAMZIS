#' PeptideCollapser Subset peptides down to a list of sequons
#'
#' @param filename dataset to be subset
#' @param kmin minimum number of glycopeptide observations Default=2
#' @param rel Boolean indicating use of relativization. Default=TRUE
#' @param sequonvector Sequon Vector
#'
#' @return Returns dataset broken down by sequon (It does not separate by protein)
#' @export
#'
#' @examples #
PeptideCollapser<-function(filename,kmin=2,rel=TRUE,sequonvector){
  datfile1<-SimDataClean(filename,kmin,rel)
  #find sequon sections
  #find glycans associated with a sequon
  #merge same glycans for a given sequon
  #return datafile with sequon plus glycan
  Rnames<-row.names(datfile1)
  Cleaned<-gsub("\\s*\\{[^\\}]+\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames))
  Outdat<-data.frame(matrix(NA,nrow=1,ncol=dim(datfile1)[2]))
  idxterm<-0
  NewRnameList<-c()
  for (j in 1:length(sequonvector)){
    Selection<-grep(sequonvector[j],Cleaned)
    Glycans<-gsub(".*\\{|\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames[Selection]))
    UniGly<-unique(Glycans)
    for (l in 1:length(UniGly)){
      RowGlyPep<-paste0(sequonvector[j],'{',UniGly[l],'}')
      GlySel<-grep(UniGly[l],Glycans)
      Outdat[idxterm+l,]<-colSums(datfile1[Selection,][GlySel,],na.rm=T)
      NewRnameList<-c(NewRnameList,RowGlyPep)
    }
    idxterm<-idxterm+length(UniGly)

  }
  row.names(Outdat)<-NewRnameList
  return(Outdat)
}
