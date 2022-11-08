
CC.scale.obj.fun <- function(param,q.old,mu.old,param.old,perc.q,perc.mu.gamma) {
  
  shape.old <- param.old[1]
  rate.old  <- param.old[2]
  
  #perturb the gamma parameters in a way that gaurentees the perc.mu change
  shape.new <- shape.old*perc.mu.gamma*param
  rate.new <- rate.old*param
  q.new <- qgamma(.99999999999999,shape.new,rate.new)
  #calcualte error between percent change and target percent change
  e.q <- (q.new/q.old - perc.q)  
  return(e.q^2)
  
}
