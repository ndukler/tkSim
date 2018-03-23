#return theoretical values given an alpha, beta, vector of times, and an initial value
#collapses to a/b * (1-e^(-bt)) when X(0) is == 0
getAbund = function(alpha,beta,time,initVal)
{
  return(exp(-time * beta) * (initVal - alpha / beta) + alpha / beta)
}
