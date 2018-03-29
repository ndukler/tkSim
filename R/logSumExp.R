#trick to allow multiplying many small numbers in real space using addition in log space
logSumExp = function(x)
{
  a = max(x)
  return(a+log(sum(exp(x-a))))
}
