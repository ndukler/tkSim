##Package Test File
library(tkSim)
bkm=basicKineticModel(times=1:30,synthRate=1:10,degRate = rep(0.5,10))
bkm=simulateData(bkm)
bkm=simulateReads(bkm,expectedLibSize=10^6,replicates=3,spikeInSizes = 200,errorModel=function(x){rep(10^6,length(x))}) #CLUGE
bkm=inferParameters(bkm)
plotParameterFit(bkm,geneIdx=1:3)
# print(test$par[1]/test$par[2])
