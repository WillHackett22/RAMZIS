#' SimPlotGGBase plots similarity distributions. It is primarily for use through other functions
#'
#' @param PlotTitle The title of the plot
#' @param dfSim PlaceHolder
#' @param RefName PlaceHolder
#' @param OverlapData PlaceHolder
#' @param Confidence PlaceHolder
#' @param legend.bool PlaceHolder
#'
#' @return PlaceHolder
#' @export
#' @import ggplot2
#'
#' @examples #
SimPlotGGBase<-function(PlotTitle,dfSim,RefName,OverlapData,Confidence=NULL,legend.bool=T){
  FPPer<-round(OverlapData$FP*100,1)
  FNPer<-round(OverlapData$FN*100,1)
  group.labels<-c("Test","Null",paste0("Alpha=",FPPer,"%"),paste0("Beta=",FNPer,"%"))
  if (RefName=="Null"){
    group.colors<-c(Test="Red",Null="Blue",Alpha="Black",Beta="Grey")
  } else if (RefName=="All"){
    group.colors<-c(Test="Red",Null="Blue",Internal1="Green",Internal2="Yellow",Alpha="Black",Beta="Grey")
    group.labels<-c("Test","Null","Internal1","Internal2",paste0("Alpha=",FPPer,"%"),paste0("Beta=",FNPer,"%"))
  } else if (RefName=="Internal1"){
    group.colors<-c(Test="Red",Internal1="Green",Alpha="Black",Beta="Grey")
    if (!is.null(Confidence)){
      group.labels[2]<-paste0("Internal=",Confidence)
    } else{
      group.labels[2]<-paste0("Internal")
    }
  } else if (RefName=="Internal2"){
    group.colors<-c(Test="Red",Internal2="Green",Alpha="Black",Beta="Grey")
    if (!is.null(Confidence)){
      group.labels[2]<-paste0("Internal=",Confidence)
    } else{
      group.labels[2]<-paste0("Internal")
    }
  }
  maxh<-max(dfSim$Dens)*1.05
  p<-ggplot2::ggplot(dfSim,aes(x=Sim,y=Dens,fill=Dist))+geom_area(position='identity')+
    scale_fill_manual(values=group.colors,labels=group.labels)+geom_line(size=0.05)+
    theme_RAMZIS+coord_cartesian(ylim=c(0.009999,maxh),xlim=c(0,1.1))+
    scale_x_continuous(name='Similarity',breaks=seq(0,1,by=.2))+
    scale_y_continuous(name='Density',breaks=seq(0,maxh,by=round(maxh/5,1)))+
    labs(title=PlotTitle,fill="Distribution",colour="Similarity")+
    geom_vline(aes(xintercept=SimilarityObj$Actual,color="Observed"),size=1)+
    scale_color_manual(values=c(Observed="Black"))
  if (!legend.bool){
    p<-p+theme(legend.position = "none")
  }
  return(p)
}
## theme_RAMZIS is the ggplot2 theme used by simplotggbase
theme_RAMZIS<-ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank(),axis.line=ggplot2::element_line(colour='black',size=1),
                    axis.text = ggplot2::element_text(size=6),axis.title = ggplot2::element_text(size=8),
                    plot.title = ggplot2::element_text(size=9),legend.title=ggplot2::element_text(size=8),
                    legend.text = ggplot2::element_text(size=6))
