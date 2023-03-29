#' BoxplotAbundances produces a paired ggplot of glycopeptide abundances
#'
#' @param filename1 The First .csv matrix RAMZIS_formatted file to be entered
#' @param filename2 The second .csv matrix RAMZIS_formatted file to be entered
#' @param normvector The normalization vectors of the files in list format. Default=list("None","None")
#' @param group1 Name of the first sample group. Corresponds to first file. Use string
#' @param group2 Name of the second sample group. Corresponds to second file. Use string
#' @param PlotTitle Title of the plot
#' @param axislabelsize size of axis text
#' @param zerofill If glycopeptides are absent from one of the samples, they will be filled in with zero values. Default=FALSE
#'
#' @return ggplot object
#' @export
#'
#' @examples #
BoxplotAbundances<-function(filename1,filename2,normvector=list("None","None"),group1,group2,PlotTitle,axislabelsize=6,zerofill=FALSE){
  dfNI1<-SimDataCleanJoint(filename1,filename2,normvector = normvector)
  df1<-dfNI1$DF1
  df2<-dfNI1$DF2
  ug<-unique(c(row.names(df1),row.names(df2)))
  if (zerofill==TRUE){
    ggd2<-data.frame(matrix(NA,nrow=length(ug)*dim(df1)[2]+length(ug)*dim(df2)[2],ncol=4))
    colnames(ggd2)<-c('Group','Sample','Identification','RelativeLogAbundance')
    ggd2$Group<-c(rep(group1,length(ug)*dim(df1)[2]),rep(group2,length(ug)*dim(df2)[2],ncol=4))
    for (j in 1:length(ug)){
      gptgt<-ug[j]
      gpabv<-GlycanSimplifier(gptgt)
      basej<-((j-1)*dim(df1)[2])
      jdx<-(basej+1):(basej+dim(df1)[2])
      ggd2$Sample[jdx]<-colnames(df1)
      ggd2$Identification[jdx]<-gpabv
      if (!(gptgt %in% row.names(df1))){
        ggd2$RelativeLogAbundance[jdx]<-0
      } else {
        ggd2$RelativeLogAbundance[jdx]<-unlist(df1[which(row.names(df1)==gptgt),])
      }
    }
    for (j in 1:length(ug)){
      gptgt<-ug[j]
      gpabv<-GlycanSimplifier(gptgt)
      basej<-(((j-1)*dim(df2)[2])+length(ug)*dim(df1)[2])
      jdx<-(basej+1):(basej+dim(df2)[2])
      ggd2$Sample[jdx]<-colnames(df2)
      ggd2$Identification[jdx]<-gpabv
      if (!(gptgt %in% row.names(df2))){
        ggd2$RelativeLogAbundance[jdx]<-0
      } else {
        ggd2$RelativeLogAbundance[jdx]<-unlist(df2[which(row.names(df2)==gptgt),])
      }
    }
  } else {
    ggd2<-data.frame(matrix(NA,nrow=dim(df1)[1]*dim(df1)[2]+dim(df2)[1]*dim(df2)[2],ncol=4))
    colnames(ggd2)<-c('Group','Sample','Identification','RelativeLogAbundance')
    ggd2$Group<-c(rep(group1,dim(df1)[1]*dim(df1)[2]),rep(group2,dim(df2)[1]*dim(df2)[2],ncol=4))
    for (j in 1:dim(df1)[1]){
      for (l in 1:dim(df1)[2]){
        idx<-(((j-1)*dim(df1)[2])+l)
        ggd2$Sample[idx]<-colnames(df1)[l]
        ggd2$Identification[idx]<-GlycanSimplifier(row.names(df1)[j])
        ggd2$RelativeLogAbundance[idx]<-df1[j,l]
      }
    }
    for (j in 1:dim(df2)[1]){
      for (l in 1:dim(df2)[2]){
        idx<-(((j-1)*dim(df2)[2])+l+dim(df1)[1]*dim(df1)[2])
        ggd2$Sample[idx]<-colnames(df2)[l]
        ggd2$Identification[idx]<-GlycanSimplifier(row.names(df2)[j])
        ggd2$RelativeLogAbundance[idx]<-df2[j,l]
      }
    }
  }
  numlines<-1.5*(length(ug)-2)
  p<-ggplot(data=ggd2,mapping=aes(x=Identification,y=RelativeLogAbundance,fill=Group))+
    geom_boxplot(lwd=0.1,outlier.size = 0.1)+theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line = element_line(colour = "black"),axis.text.x = element_text(angle=55,vjust=1,hjust=1),
          text=element_text(size=axislabelsize),plot.title = element_text(angle=0,size=10),legend.key.width = unit(0.25,'cm'))+
    geom_vline(xintercept=seq(1.5,numlines),linewidth=0.1)+ggtitle(PlotTitle)+xlab('Peptide{Hex;HexNAc;Fuc;Neu5AC;Sulfation}')
  return(p)
}
