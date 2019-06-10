# #RamzisComparisonTest
# 
# #UnitTest k<-1
# k<-1
# df1<-data.frame(matrix(k,nrow=20,ncol=5))
# row.names(df1)<-paste0('Ident',1:20)
# 
# wmatrix<-matrix(0,nrow = 325,ncol=21)
# bmatrixab<-matrix(0,nrow=325,ncol=20)
# bmatrixba<-matrix(0,nrow=325,ncol=20)
# whold1<-WithinSimBootstrap2(df1)
# wmatrix[,1]<-whold1$Summary$Tanimoto
# 
# for (j in 1:20){
#   df2<-data.frame(matrix(k,nrow=20,ncol=5))
#   row.names(df2)<-paste0('Ident',(j+1):(j+20))
#   whold2<-WithinSimBootstrap2(df2)
#   wmatrix[,(1+j)]<-whold2$Summary$Tanimoto
#   bholdab<-BetweenSimBootstrap2(df1,df2,whold1$Boot)
#   bmatrixab[,j]<-bholdab$Summary$Tanimoto
#   bholdba<-BetweenSimBootstrap2(df2,df1,whold2$Boot)
#   bmatrixba[,j]<-bholdba$Summary$Tanimoto
# }
# 
# 
# #UnitTest k<-.5
# k<-.5
# df1<-data.frame(matrix(k,nrow=20,ncol=5))
# row.names(df1)<-paste0('Ident',1:20)
# 
# wmatrix2<-matrix(0,nrow = 325,ncol=21)
# bmatrixab2<-matrix(0,nrow=325,ncol=20)
# bmatrixba2<-matrix(0,nrow=325,ncol=20)
# whold1<-WithinSimBootstrap2(df1)
# wmatrix2[,1]<-whold1$Summary$Tanimoto
# 
# for (j in 1:20){
#   df2<-data.frame(matrix(k,nrow=20,ncol=5))
#   row.names(df2)<-paste0('Ident',(j+1):(j+20))
#   whold2<-WithinSimBootstrap2(df2)
#   wmatrix2[,(1+j)]<-whold2$Summary$Tanimoto
#   bholdab<-BetweenSimBootstrap2(df1,df2,whold1$Boot)
#   bmatrixab2[,j]<-bholdab$Summary$Tanimoto
#   bholdba<-BetweenSimBootstrap2(df2,df1,whold2$Boot)
#   bmatrixba2[,j]<-bholdba$Summary$Tanimoto
# }
# #Near Linear  Effect of shared 
# 
# 
# wmatrix3<-matrix(0,nrow = 325,ncol=21)
# bmatrixab3<-matrix(0,nrow=325,ncol=20)
# bmatrixba3<-matrix(0,nrow=325,ncol=20)
# wmatrix3[,1]<-whold1$Summary$Tanimoto
# 
# for (j in 1:19){
#   df2<-data.frame(matrix(k,nrow=(20-j),ncol=5))
#   row.names(df2)<-paste0('Ident',1:(20-j))
#   whold3<-WithinSimBootstrap2(df2)
#   wmatrix3[,(1+j)]<-whold2$Summary$Tanimoto
#   bholdab<-BetweenSimBootstrap2(df1,df2,whold1$Boot)
#   bmatrixab3[,j]<-bholdab$Summary$Tanimoto
#   bholdba<-BetweenSimBootstrap2(df2,df1,whold2$Boot)
#   bmatrixba3[,j]<-bholdba$Summary$Tanimoto
# }




Adet1<-seq(2000,40000,by=2000)
Bdet1<-seq(40000,2000,by=-2000)

n<-5
g<-20
p<-.65
w<-5

TD1<-TheoreticalDataGenerator(n,g,p,w,alp=Adet1,bet=Bdet1)
TD2<-TheoreticalDataGenerator(n,g,p,w,alp=Bdet1,bet=Adet1)

wh1<-WithinSimBootstrap2(TD1$Dataset,kmin=1)
wh2<-WithinSimBootstrap2(TD2$Dataset,kmin=1)
b12<-BetweenSimBootstrap2(TD1$Dataset,TD2$Dataset,wh1$Boot,kmin=1)
b21<-BetweenSimBootstrap2(TD2$Dataset,TD1$Dataset,wh2$Boot,kmin=1)


boxplot(wh1$RankInfoFinal[,order(colMeans(wh1$RankInfoW))],las=2)
boxplot(wh2$RankInfoFinal[,order(colMeans(wh2$RankInfoW))],las=2)

boxplot(b21$RankInfoFinal[,order(colMeans(b21$RankInfoW))],las=2)
boxplot(b12$RankInfoFinal[,order(colMeans(b12$RankInfoW))],las=2)

#RANK TROUBLESHOOTING
b21rel<-b21$RankInfoFinal/wh2$RankInfoFinal
b12rel<-b12$RankInfoFinal/wh1$RankInfoFinal
bebop<-bcomp<-data.frame(matrix(0,nrow=65,ncol=20))
bcomp<-data.frame(matrix(0,nrow=20,ncol=2))
colnames(bcomp)<-c('Mean','SD')
thold<-rep(0,20)
for (j in 1:20){
  bcomp$Mean[j]<-mean(c(unlist(b12rel[j]),unlist(b21rel[j])),na.rm=T)
  bcomp$SD[j]<-sd(c(unlist(b12rel[j]),unlist(b21rel[j])),na.rm=T)
  bebop[,j]<-abs(log(unlist(b12rel[j]))-log(unlist(b21rel[j])))
  thold[j]<-t.test(TD1$Dataset[j,],TD2$Dataset[j,])$p.value
}



#Checking against p value
RankInfoAggregate<-rep(0,200)
QualityInfoAggregate<-rep(0,200)
PValueAggregate<-rep(0,200)
Relative<-matrix(0,nrow=65,ncol=20)
RankInfoAggregateL<-rep(0,200)
QualityInfoAggregateL<-rep(0,200)
PValueAggregate<-rep(0,200)
RelativeL<-matrix(0,nrow=65,ncol=20)

for (l in 1:10){
  TD1<-TheoreticalDataGenerator(n,g,p,w)
  TD2<-TheoreticalDataGenerator(n,g,p,w)
  CompHold<-AdvancedSimPlot2(TD1$Dataset,TD2$Dataset,kmin=2,rel=TRUE,simtype='Tanimoto',PlotTitle = 'Example','DF1','DF2')
  while (dim(CompHold$Internal1$RankInfo)[2]<20 | dim(CompHold$Internal2$RankInfo)[2]<20){
    TD1<-TheoreticalDataGenerator(n,g,p,w)
    TD2<-TheoreticalDataGenerator(n,g,p,w)
    CompHold<-AdvancedSimPlot2(TD1$Dataset,TD2$Dataset,kmin=2,rel=TRUE,simtype='Tanimoto',PlotTitle = 'Example','DF1','DF2')
  }
  CompRank<-RankAggregate(CompHold$Internal1,CompHold$Internal2,CompHold$External12,CompHold$External21)
  CompRankL<-RankAggregate(CompHold$Internal1,CompHold$Internal2,CompHold$External12,CompHold$External21,RankMethod = 'Log')
  RankInfoAggregate[((l-1)*20+1):(l*20)]<-unlist(sqrt(colSums(CompRank$RankData^2)))
  RankInfoAggregateL[((l-1)*20+1):(l*20)]<-unlist(sqrt(colSums(CompRankL$RankData^2)))
  QualityInfoAggregate[((l-1)*20+1):(l*20)]<-unlist(CompRank$QualityCheck[8,])
  QualityInfoAggregateL[((l-1)*20+1):(l*20)]<-unlist(CompRankL$QualityCheck[8,])
  PTemp<-rep(0,20)
  for (j in 13:20){
    if (sum(TD1$Dataset[j,]==1)>1 & sum(TD2$Dataset[j,]==1)>1){
      PTemp[j]<-1
    } else {
      PTemp[j]<-t.test(TD1$Dataset[j,],TD2$Dataset[j,])$p.value
    }
  }
  PValueAggregate[((l-1)*20+1):(l*20)]<-PTemp
  Relative[,(2*(l-1)+1):(2*l)]<-matrix(unlist(CompRank$RankTotal),nrow=65,ncol=2)
}

plot(log(RankInfoAggregate),PValueAggregate,pch=16,col=rep(2,20)+2*QualityInfoAggregate)
abline(v=0,col=4,lty=6,lwd=1.5)
abline(h=0.05,col=2,lty=6,lwd=1.5)
plot(log(RankInfoAggregate)[which(QualityInfoAggregate==1)],log(PValueAggregate[which(QualityInfoAggregate==1)]),main='Ranking Information v P-Value',xlab='Ln(SQRT(rM^2+rS^2))',ylab = 'ln(P-value)',pch=16,col=rep(2,length(which(QualityInfoAggregate==1)))+2*(log(RankInfoAggregate)[which(QualityInfoAggregate==1)]<0))
abline(v=0,col=4,lty=6,lwd=1.5)
abline(h=log(0.05/seq(1,40)),col=16,lty=3,lwd=.5)
abline(h=log(0.05/seq(41,150)),col=18,lty=3,lwd=.5)
text(-2,-15,labels=c("26.7% Recommended"),col=4)
text(-2,log(0.05),labels=c("alpha=0.05"),col=16)



bebop<-data.frame(matrix(0,nrow=4225,ncol=20))
pos<-matrix(0,nrow=65,ncol=65)
pos[1,]<-1:65
for (j in 2:65){
  pos[j,]<-c(j:65,1:(j-1))
}
for (j in 1:65){
  bebop[((j-1)*65+1):(j*65),]<-abs(log(unlist(b12rel[j,]))-log(unlist(b21rel[pos[j,],])))
}

b21relex<-data.frame(matrix(0,nrow=1625,ncol=20))
b12relex<-data.frame(matrix(0,nrow=1625,ncol=20))


pos<-matrix(1:25,nrow=5,ncol=5)
for (l in 1:5){
  for (k in 1:5){
    b12relex[((pos[k,l]-1)*65+1):(pos[k,l]*65),]<-b12$RankInfoW[(seq(0,320,by=5)+k),]/wh1$RankInfoW[(seq(0,320,by=5)+l),]
  }
}
for (j in 1:20){
  b12relex[is.infinite(b12relex[,j]),j]<-NA
  b12relex[which(b12relex[,j]=='NaN'),j]<-NA
}

plot(log(b12rel[order(b12rel)]),las=2)
hist(b12rel[order(b12rel)],las=2,prob=T,breaks=c(seq(0,4,by=.1),seq(5,30,by=5)))
bopboo<-colMeans(b12rel,na.rm=T)*colMeans(b21rel,na.rm=T)/(colMeans(b12rel,na.rm=T)^2+colMeans(b21rel,na.rm=T)^2-colMeans(b21rel,na.rm=T)*colMeans(b12rel,na.rm=T))


b12mat<-matrix(0,nrow=dim(b12$RankInfoFinal)[1]/n,ncol=dim(b12$RankInfoFinal)[2])
colnames(b12mat)<-colnames(b12$RankInfoFinal)
wh1mat<-matrix(0,nrow=dim(wh1$RankInfoFinal)[1]/n,ncol=dim(wh1$RankInfoFinal)[2])
colnames(wh1mat)<-colnames(wh1$RankInfoFinal)

for (j in 1:65){
  b12mat[j,]<-colMeans(b12$RankInfoW[((j-1)*5+1):(j*5),])
  wh1mat[j,]<-colMeans(wh1$RankInfoW[((j-1)*5+1):(j*5),])
}
boxplot(b12mat,las=2)
boxplot(b12$RankInfoW,las=2)
boxplot(wh1mat,las=2)
boxplot(wh1$RankInfoW,las=2)
boxplot(wh1$Summary$Tanimoto,b12$Summary$Tanimoto)
boxplot(wh1$Hold$Tanimoto,b12$Hold$Tanimoto)
boxplot(wh2$Summary$Tanimoto,b21$Summary$Tanimoto)
boxplot(wh2$Hold$Tanimoto,b21$Hold$Tanimoto)



boxplot(b21$RankInfoFinal,las=2)
boxplot(b21$RankInfoW,las=2)





bRefs<-unique(c(colnames(b12$RankInfoFinal),colnames(wh1$RankInfoFinal)))
b12relm<-data.frame(matrix(0,nrow=dim(b12$RankInfoFinal)[1],ncol=length(bRefs)))
colnames(b12relm)<-bRefs
for (j in 1:length(bRefs)){
  b12relm[,bRefs[j]]<-b12$RankInfoFinal[,bRefs[j]]/wh1$RankInfoFinal[,bRefs[j]]
}
for (j in 1:dim(b12relm)[2]){
  b12relm[is.infinite(b12relm[,j]),j]<-NA
  b12relm[is.na(b12relm[,j]),j]<-NA
  
}

boxplot(b12relm,las=2)

b21relm<-b21$RankInfoW/wh2$RankInfoW
for (j in 1:dim(b12relm)[2]){
  b12relm[is.infinite(b12relm[,j]),j]<-NA
  b12relm[is.na(b12relm[,j]),j]<-NA
  
}

boxplot(b12relm,las=2)



#uniformly vary the alpha criterion to change the mean











