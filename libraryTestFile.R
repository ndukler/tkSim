##Package Test File
library(tkSim)
bkm=basicKineticModel(times=0:30,synthRate=1:10,degRate = rep(0.5,10))
bkm=simulateData(bkm)
bkm=simulateReads(bkm,3,1,errorModel=function(x){rep(2,length(x))}) #CLUGE
head(bkm@data)
