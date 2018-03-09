x=rnbinom(10000,mu=10,size=3)

nll.nbinom <- function(x,mu.val,size.val=3){
  ## log(P(x|mu,size))
  log.prob=dnbinom(x,mu = mu.val,size = size.val,log = TRUE)
  ## Compute the j
  return(-sum(log.prob))
}

mu.vec=seq(0.01,100,00.01)
ll.vec=numeric(length(mu.vec))
for(i in 1:length(mu.vec)){
  ## Compute P(X|theta) for many values of theta
  ll.vec[i]=nll.nbinom(x,mu.vec[i])
}

plot(ll.vec,type="l")
abline(v=which.min(ll.vec))
mu.vec[which.min(ll.vec)]

##
plot(exp(-ll.vec)/sum(exp(-ll.vec)))


nll2.nbinom <- function(par=c(2,5),data){
  ## log(P(x|mu,size))
  log.prob=dnbinom(data,mu = par[1],size = par[2],log = TRUE)
  ## Compute the j
  return(-sum(log.prob))
}


optim(par = c(1,8),fn = nll2.nbinom,data=x,method = "L-BFGS-B",lower = c(10^-5,10^-5))
