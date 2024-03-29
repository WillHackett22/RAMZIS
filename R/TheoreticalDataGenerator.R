#' TheoreticalDataGenerator makes a theoretical dataset of relative abundances
#'
#' @param n number of samples
#' @param g number of identifications
#' @param p average presence
#' @param samplevar controls variance; lower is less variant if given a vector length equal to g the variance is controlled by identification; if given a vector length equal to n the variance is controlled by sample. default=0
#' @param alp controls mean of the beta distribution used for sampling, if given a vector length g each identification has its own average
#' @param betcontrols mean of the beta distribution used for sampling, if given a vector length g each identification has its own average
#' @param maxim Boolean maximizing each column's identification. Default=FALSE
#' @param ratio_sig_N number of significant digits to round alpha and beta to after finding the means. requires integers. default=2
#' @param rat_mult controls the variance of the beta distribution production. default=10
#' @param fakeGPLeader name of identifications produced
#'
#' @return Theoretical dataset
#' @export
#' @import MASS
#'
#' @examples
#' #data<-TheoreticalDataGenerator(5,20,0.7,5.5,alp=seq(100,2000,by=100),bet=seq(2000,100,by=-100))
#' #data$Dataset is a dataset with 5 samples and 20 identifications (glycopeptides)
#' #the means of data will trend from 0 to 1
TheoreticalDataGenerator<-function(n,g,p,samplevar=1,alp=2,bet=2,maxim=FALSE,ratio_sig_N=2,rat_mult=10,fakeGPLeader="Identification"){
  #n samples
  #g identifications
  #p average presence
  ## p determines the presence distribution that appears across the sample
  ## binomial distribution rate of observation is proportional to signal intensity
  ## eg P(1,n,p_adj)=#observations. p_adj=p*exp(x-mean(x))
  #samplevar term controlling the variance in the data
  #alp and bet are terms for general data distribution
  datagen<-data.frame(matrix(0,nrow=g,ncol=n))
  rowsampname<-paste0(fakeGPLeader,1:g)
  row.names(datagen)<-rowsampname
  colnames(datagen)<-paste0('Sample',1:n)

  #initialize datavalues
  datavals<-rep(0,g)
  #scalar or vector for initialization
  if (length(alp)==1 & length(bet)==1){
    datavals<-rbeta(g,alp,bet)
  } else if (length(alp)==g & length(bet)==g){
    for (j in 1:g){
      datavals[j]<-rbeta(1,alp[j],bet[j])
    }
  }

  #USE MASS:::.rat(mean)$rat to get alphas and betas
  #use rounding to generate betas for each value so they all are within 0-1
  datmean<-mean(datavals)
  for (j in 1:g){
    if (length(samplevar)==n & length(samplevar)!=g){
      for (l in 1:n){
        #make sample variant data
        FractionObj<-MASS:::.rat(signif(datavals[j],ratio_sig_N))$rat
        FractionObj<-FractionObj*(rat_mult^(samplevar[l]))
        bet0<-FractionObj[2]-FractionObj[1]
        alp0<-FractionObj[1]
        datagen[j,l]<-rbeta(1,alp0,bet0)
      }
    } else if (length(samplevar)==g) {
      #make glycopeptide variant data
      FractionObj<-MASS:::.rat(signif(datavals[j],ratio_sig_N))$rat
      FractionObj<-FractionObj*(rat_mult^(samplevar[j]))
      bet0<-FractionObj[2]-FractionObj[1]
      alp0<-FractionObj[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    } else {
      #make uniform variant data
      FractionObj<-MASS:::.rat(signif(datavals[j],ratio_sig_N))$rat
      FractionObj<-FractionObj*(rat_mult^(samplevar[1]))
      bet0<-FractionObj[2]-FractionObj[1]
      alp0<-FractionObj[1]
      datagen[j,]<-rbeta(n,alp0,bet0)
    }
    if (p!=1){
      padj<-p*exp(mean(unlist(datagen[j,])-datmean))
      padj[padj>=1]<-0.99999999
      padj[padj<=0]<-0.00000001
    } else {
      padj<-p
    }

    keepnum<-rbinom(1,(n-1),padj)+1
    datagen[j,-match(sample(unlist(datagen[j,]),keepnum),datagen[j,])]<-0
  }
  if (maxim==TRUE){
    datagen[which(rowMeans(datagen)==max(rowMeans(datagen))),]<-rep(1,n)
  }
  FinalOut<-list("Dataset"=datagen,"parameters"=c(n,g,p),"variance_parameters"=list("VarExp"=samplevar,"varMult"=rat_mult,"SigFig"=ratio_sig_N),"alp"=alp,"bet"=bet,"DistributionCenters"=datavals)
  return(FinalOut)
}
