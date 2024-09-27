CC.scale <- function(emission.fits.site,months,dates.sim,n.sites,perc.mu,perc.q,thshd.prcp,qq.month) {
  
  #this function takes fitted emission distributions to each site by month, and returns those parameters as 
  #an array of parameters with dimensions (n.sim, n.sites,n.parameters) 
  #it also enables a Clausius-Clapyeron scaling for the mean and the qqth quantile, and returns a similar
  #array of parameters for emission distributions with that scaling
  
  #Arguments:
  #emission.fits.site = a list of length 2 (gamma and gpd fits), with each list containing an array of dimension (n.parameters x n.months x n.sites) that contains emission distribution parameters for each site and and month
  #months = a vector of calendar months that define the season
  #dates.sim = a vector of dates associated with the simulation period
  #n.sites = the number of individual sites
  #perc.mu = the percentage by which to scale the mean
  #perc.q = the percentage by which to scale the upper tail for CC scaling
  #thshd.prcp = the threshold used to separate gamma from GPD
  #qq.month = the frequency prcp under threshold by month and site
  
  #month sequence for the simulation period
  months.sim <- as.numeric(format(dates.sim,"%m"))
  n.sim <- length(dates.sim)  #length of the simulation period
  
  #original fits  
  emission.old1 <- emission.fits.site[[1]]
  
  # monthly-only #
  #find the mean and upper quantile of the different distributions, and concatenate them with the parameters into a single array
  
  #q.max.prcp = the quantile for corresponding max prcp (month and site) under old gamma fits
  q.max.prcp <- emission.fits.site[[3]]
  
  q.old <- sapply(1:n.sites,function(x,z,m) {
    qgamma(q.max.prcp[m,x],shape=z[1,m,x],rate=z[2,m,x])
  }, z=emission.old1,m=months
  )
  mu.old <- apply(emission.old1,c(2,3),function(x) {x[1]/x[2]})
  
  #GPD means
  gpd_mu <- thshd.prcp + apply(emission.fits.site[[2]],2,function(x) {x[1]/(1-x[2])})
  gpd_mu.mat <- matrix(rep(gpd_mu,12),nrow=12,byrow=T)
  perc.mu.gamma <- (perc.mu*qq.month*mu.old - (perc.q-perc.mu)*(1-qq.month)*gpd_mu.mat)/(qq.month*mu.old)
  
  q.mu.old <- abind(q.old,mu.old,emission.old1,perc.mu.gamma,q.max.prcp,along=1)
  
  #adjust distribution for new mean, upper quantile
  emission.new1 <- apply(q.mu.old,c(2,3), function(x) {
    start.par <- 1
    lowerb <- 0.0001
    upperb <- 10
    opt <- optim(par=start.par,CC.scale.obj.fun,
                 q.old=x[1],mu.old=x[2],
                 param.old=x[3:4],
                 perc.q=perc.q,perc.mu.gamma=x[5],
                 q.max.prcp=x[6],
                 method="L-BFGS-B",lower=lowerb,upper=upperb)
    shape.new <- x[3]*x[5]*opt$par
    rate.new <- x[4]*opt$par
    return(c(shape.new,rate.new))
  })
  return(list(emission.old1,emission.new1))
  
}