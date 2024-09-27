
CC.scale.obj.fun <- function(param,q.old,mu.old,param.old,perc.q,perc.mu.gamma,q.max.prcp) {
  
  shape.old <- param.old[1]
  rate.old  <- param.old[2]
  
  #perturb the gamma parameters in a way that guarantees the perc.mu change
  shape.new <- shape.old*perc.mu.gamma*param
  rate.new <- rate.old*param
  
  q.new <- qgamma(q.max.prcp,shape=shape.new,rate=rate.new)
  
  #calculate error between percent change and target percent change
  e.q <- (q.new/q.old - perc.q)  
  return(e.q^2)
  
}
