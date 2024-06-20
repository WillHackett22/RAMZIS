#' PeptideCollapser turns one protein's glycreread output into a break down by glycosite.
#'
#' @param filename dataset to be subset
#' @param Outfilename Defualt=NULL. If given string, writes each file with format Outfilename_sequonvectori.csv
#' @param verbose Default=T. Sends output to console as list of dataframes
#' @param sequonvector Sequon Vector to collapse peptides on. If Null, it goes by exact value. Be sure your sequons are distinct
#' @param gsep Default="Curly" The delineators that separate the glycan from the peptide. Alternate Option="Square"
#'
#' @return Returns dataset broken down by sequon
#' @export
#'
#' @examples #
PeptideCollapser<-function(filename,Outfilename=NULL,verbose=T,sequonvector=NULL,gseps="Curly"){
  df<-read.csv(filename,row.names = 1)
  #find sequon sections
  #find glycans associated with a sequon
  #merge same glycans for a given sequon
  #return datafile with sequon plus glycan
  Rnames<-row.names(df)
  if (gseps=="Curly"){
    Cleaned<-gsub("\\s*\\{[^\\}]+\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames))
  } else if (gseps=="Square"){
    Cleaned<-gsub("\\s*\\[[^\\]+\\]","",gsub("\\s*\\([^\\)]+\\)","",Rnames))
  }

  if (is.null(sequonvector)){
    sequonvector<-unique(Cleaned)
    exact=T
  } else {
    exact=F
  }
  outlist<-list()
  for (j in 1:length(sequonvector)){
    if (exact==T){
      Selection<-which(Cleaned==sequonvector[j])
    } else {
      Selection<-grep(sequonvector[j],Cleaned)
    }
    if (length(Selection)>0){
      if (gseps=="Curly"){
        Glycans<-gsub(".*\\{|\\}","",gsub("\\s*\\([^\\)]+\\)","",Rnames[Selection]))
      } else if (gseps=="Square"){
        Glycans<-gsub(".*\\[|\\]","",gsub("\\s*\\([^\\)]+\\)","",Rnames[Selection]))
      }
      UniGly<-unique(Glycans)
      out<-data.frame(matrix(NA,nrow=length(UniGly),ncol=ncol(df)))
      colnames(out)<-colnames(df)
      if (gseps=="Curly"){
        sepg=c("{","}")
      } else if (gseps=="Square"){
        sepg=c("[","]")
      }
      row.names(out)<-paste0(sequonvector[j],sepg[1],UniGly,sepg[2])
      for (l in 1:length(UniGly)){
        GlySel<-grep(UniGly[l],Glycans)
        out[paste0(sequonvector[j],sepg[1],UniGly[l],sepg[2]),]<-colSums(df[Selection,][GlySel,],na.rm=T)
      }
      outlist[[j]]<-out
      if (!is.null(Outfilename)){
        utils::write.csv(out,paste0(Outfilename,"_",sequonvector[j],'.csv'))
      }
    }
  }
  if (verbose){
    return(outlist)
  }
}
