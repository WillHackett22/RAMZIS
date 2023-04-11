#' SimPlotGGBase plots similarity distributions. It is primarily for use through other functions
#'
#' @param PlotTitle The title of the plot
#' @param dfSim TIDY Formatted Similarity Distributions
#' @param RefName Reference Distribution Name. ("Null","Internal1","Internal2","All")
#' @param OverlapData Overlap data from OverlapCalculator
#' @param ObservedSimObj List produced by ObservedSimilarityStats() Default=NULL
#' @param Confidence Internal Confidence Score. Do not touch if doing TestNull Plot. DEFAULT=NULL
#' @param legend.bool Boolean indicating whether or not to produce the legend. Default=T
#' @param themepick Default=theme_RAMZIS. Theme object from ggplot2 for formatting.
#' @param zopt Indicates use of Z-score for simulation validity over percentile. Default=T
#' @param xbounds Default=c(0,1.1)
#'
#' @return A similarity plot using ggplot2
#' @export
#' @import ggplot2
#'
#' @examples #
SimPlotGGBase<-function(PlotTitle,dfSim,RefName,OverlapData,ObservedSimObj=NULL,Confidence=NULL,legend.bool=T,themepick=theme_RAMZIS,zopt=T,xbounds=c(0,1.1),ybounds="Auto",...){
  FPPer<-round(OverlapData$FP*100,1)
  FNPer<-round(OverlapData$FN*100,1)
  if (!is.null(ObservedSimObj)){
    if (zopt){
      z<-round(ObservedSimObj$ZScore,2)
      obsqua<-paste0(', Z=',z)
    } else {
      per<-round(ObservedSimObj$Percentile,2)
      obsqua<-paste0(', %=',per)
    }
    obsloc<-round(ObservedSimObj$Actual,2)
    ObsLabel<-paste0(obsloc,obsqua)
  } else {
    ObsLabel<-"Observed"
  }
  if (RefName=="Null"){
    dfSim$Dist<-factor(dfSim$Dist,levels=c("Test","Null","Alpha","Beta"))
    group.labels<-c("Test","Null",paste0("Alpha=",FPPer,"%"),paste0("Beta=",FNPer,"%"))
    group.colors<-c(Test="Red",Null="Blue",Alpha="Black",Beta="Grey")
  } else if (RefName=="All"){
    dfSim$Dist<-factor(dfSim$Dist,levels=c("Test","Null","Internal1","Internal2","Alpha","Beta"))
    group.colors<-c(Test="Red",Null="Blue",Internal1="Green",Internal2="Yellow",Alpha="Black",Beta="Grey")
    group.labels<-c("Test","Null","Internal1","Internal2",paste0("Alpha=",FPPer,"%"),paste0("Beta=",FNPer,"%"))
  } else if (RefName=="Internal1"){
    dfSim$Dist<-factor(dfSim$Dist,levels=c("Test","Internal1","Alpha","Beta"))
    group.labels<-c("Test","Internal1",paste0("Alpha=",FPPer,"%"),paste0("Beta=",FNPer,"%"))
    group.colors<-c(Test="Red",Internal1="Green",Alpha="Black",Beta="Grey")
    if (!is.null(Confidence)){
      group.labels[2]<-paste0("Internal=",Confidence)
    } else{
      group.labels[2]<-paste0("Internal")
    }
  } else if (RefName=="Internal2"){
    dfSim$Dist<-factor(dfSim$Dist,levels=c("Test","Internal2","Alpha","Beta"))
    group.labels<-c("Test","Internal2",paste0("Alpha=",FPPer,"%"),paste0("Beta=",FNPer,"%"))
    group.colors<-c(Test="Red",Internal2="Green",Alpha="Black",Beta="Grey")
    if (!is.null(Confidence)){
      group.labels[2]<-paste0("Internal=",Confidence)
    } else{
      group.labels[2]<-paste0("Internal")
    }
  }
  maxh<-max(dfSim$Dens)*1.05
  if (any(ybounds=="Auto")){
    ybounds<-c(0.045*maxh,maxh)
  }
  p<-ggplot2::ggplot(dfSim,ggplot2::aes(x=Sim,y=Dens,fill=Dist))+
    ggplot2::geom_area(position='identity')+
    ggplot2::scale_fill_manual(values=group.colors,labels=group.labels)+
    ggplot2::geom_line(size=0.05)+
    themepick+ggplot2::coord_cartesian(ylim=ybounds,xlim=xbounds)+
    ggplot2::scale_x_continuous(name='Similarity',breaks=seq(0,1,by=.2))+
    ggplot2::scale_y_continuous(name='Density',breaks=seq(0,maxh,by=round(maxh/5,1)))+
    ggplot2::labs(title=PlotTitle,fill="Distribution",colour="Observed")+
    ggplot2::geom_vline(ggplot2::aes(xintercept=ObservedSimObj$Actual,color="Observed"),linewidth=0.75)+
    ggplot2::scale_color_manual(values=c(Observed="Black"),labels=ObsLabel)
  if (!legend.bool){
    p<-p+ggplot2::theme(legend.position = "none")
  }
  return(p)
}
## theme_RAMZIS is the ggplot2 theme used by simplotggbase
theme_RAMZIS<-ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank(),axis.line=ggplot2::element_line(colour='black',linewidth=1),
                    axis.text = ggplot2::element_text(size=6),axis.title = ggplot2::element_text(size=7),
                    plot.title = ggplot2::element_text(size=8),legend.title=ggplot2::element_text(size=7),
                    legend.text = ggplot2::element_text(size=6))
