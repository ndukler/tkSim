##Package Test File
library(tkSim)
bkm=basicKineticModel(times=1:30,synthRate=1:10,degRate = rep(0.5,10))
bkm=simulateData(bkm)
bkm=simulateReads(bkm,expectedLibSize=10^6,replicates=1,errorModel=function(x){rep(10^6,length(x))}) #CLUGE
bkm@data[10,] = rnbinom(30,mu=10,size=3)
plot(bkm@data[10,])
test=inferParameters(bkm)
print(test)
print(test$par[1]/test$par[2])
