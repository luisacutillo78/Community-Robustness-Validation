library('gprege')
source('~/GPranking/callgpC.R')

#####Facebook 348 Null Model Rewire

case='Facebook_fastgreedy'
outpat='~/GPranking/OutGp/'

path='~/'
files=list.files(path,pattern='BATS')
bf=NULL
for (i in c(1:length(files))){
  filename=paste(path,files[i],sep='')
  bf[i]=callgp(filename)
}


