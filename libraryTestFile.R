##Package Test File
library(tkSim)
source('~/Documents/GitHub/tkSim/R/llFactory.R')
source('~/Documents/GitHub/tkSim/R/logSumExp.R')
source('~/Documents/GitHub/tkSim/R/nllFactory.R')
source('~/Documents/GitHub/tkSim/R/getAbund.R')

library(tkSim)
bkm=basicKineticModel(times=0:30,synthRate=1:10,degRate = rep(0.5,10))
bkm=simulateData(bkm)
bkm=simulateReads(bkm,expectedLibSize=10^6,replicates=3,spikeInSizes = 200,dispersionModel=function(x){rep(10^3,length(x))},dispByGene=F) #CLUGE
bkm = estimateDispersions(bkm) #update disperson model based on @data
bkm=inferParameters(bkm)
plotParameterFit(bkm,geneIdx=1:3)


bkm=calculatePosteriors(bkm,alphaRange=c(.25,2))
plotPosteriors(bkm,3)
plotPosteriors(object=bkm,geneIdx=3,alphaRange=c(3000,4500),betaRange=c(.45,.55),relative=F,recalculate=T,paramSpaceSize=10^5,dispByGene=T)
# print(test$par[1]/test$par[2])

##test using program with experimental data
bkm2 = basicKineticModel(data=bkm@data,spikeIns=bkm@spikeIns,expMetadata = bkm@expMetadata)
bkm2 = estimateDispersions(bkm2)
bkm2=inferParameters(bkm2)
plotParameterFit(bkm2,geneIdx=1:3)
bkm2=calculatePosteriors(bkm2,alphaRange=c(.25,2))
plotPosteriors(bkm2,3)
plotPosteriors(object=bkm2,geneIdx=3,alphaRange=c(3000,4500),betaRange=c(.45,.55),relative=F,recalculate=T,paramSpaceSize=10^5,dispByGene=T)

##test for inference errors due to low sequencing depth causing most data to be zero
srate=2^(0:10)
drate=round(2^(seq(-0.01,-5,-0.5)),3)
rate.comb=expand.grid(s=srate,d=drate)
## Build models and simulate data
bkm=basicKineticModel(times=c(2,4,8,16,240),synthRate=rate.comb$s,degRate = rate.comb$d)
bkm=simulateData(bkm)
## Simulate reads at low coverage
bkm=simulateReads(bkm,expectedLibSize=nrow(rate.comb)*10,replicates=5,spikeInSizes = 200,dispersionModel=function(x){rep(10^5,length(x))}) #CLUGE
bkm=inferParameters(bkm) #should throw errors

