
callgp<-function(filename){
  
  
  #inizializzazione
gpregeOptions = list(indexRange=(1:2), explore=FALSE, exhaustPlotRes=30, exhaustPlotLevels=10,
                     exhaustPlotMaxWidth=100, iters=100, labels=rep(FALSE,2), display=FALSE)

AA=read.table(filename,header = FALSE)
MA=as.matrix(AA[,c(2:dim(AA)[2])])
vt=unique(as.numeric(MA[1,]))
ntimes=length(vt)
stdv=NULL
varv=NULL
for (i in c(1:ntimes)){
  ind=which(MA[1,]==vt[i])
  stdv[i]=sd(MA[2,ind])
  varv[i]=var(MA[2,ind])
}
#sigmaest=mean(stdv)
GlobalVar=var(MA[2,])
SigNoise=mean(varv)/GlobalVar
if (SigNoise>1)SigNoise=1
##########

#SigNoise=1-var(MA[2,])
sigmaest=1-SigNoise
#mod='08'
#ntimes=50

gpregeOptions$inithypers <- matrix( c(
      1/1000,	0,	1
      ,1/ntimes,	sigmaest, SigNoise
    ), ncol=3, byrow=TRUE)
#gpregeOptions$inithypers <- matrix( c(
 # 1/1000,  0,	1
  #,1/ntimes,	0.8,0.2
#), ncol=3, byrow=TRUE)


dat=read.table(filename, header = FALSE, sep = "\t") 
dvet=t(data.matrix(dat[1,2:211]))
dd=(data.matrix(dat[2,2:211]))
rownames(dd)='VI'
colnames(dd)=dvet
datadum=rbind(dd,dd)
gpregeOutput <- gprege(data=datadum, inputs=dvet, gpregeOptions=gpregeOptions)
bf=gpregeOutput$rankingScores[1]
return(bf)
}