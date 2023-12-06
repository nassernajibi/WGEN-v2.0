
perturb.climate <- function(prcp.site.sim,
                            tmin.site.sim,
                            tmax.site.sim,
                            emission.fits.site,months,dates.sim,n.sites,
                            qq,perc.mu,perc.q,Sbasin,cur.jitter,cur.tc.min,cur.tc.max,num.iter,thshd.prcp,qq.month) {
  
  #this function takes the simulated prcp, tmax, and tmin from the daily weather generator and perturbs the data to 
  #1) impose thermodynamic climate changes
  #2) impose copula-based jitters
  
  #Arguemnts:
  #prcp.site.sim = a list of length num.iter, with each list element containing a matrix of simulated daily prcp (rows: days; cols: sites)
  #tmax.site.sim = a list of length num.iter, with each list element containing a matrix of simulated daily tmax (rows: days; cols: sites)
  #tmin.site.sim = a list of length num.iter, with each list element containing a matrix of simulated daily tmin (rows: days; cols: sites)
  #emission.fits.site = a list of length 2 (gamma and gpd fits), with each list containing an array of dimension (n.parameters x n.months x n.sites) that contains emission distribution parameters for each site and and month
  #months = a vector of the calendar months included in each year
  #dates.sim = a vector of dates over the simulation period
  #n.sites = the number of sites across the basin
  #qq = quantile to anchor the CC-scaling
  #perc.mu = percent change in the mean for each month of non-zero prcp for CC-scaling
  #perc.q = percent change in the qqth quantile for each month of non-zero prcp for CC-scaling
  #Sbasin = the spearman correlation matrix (n.site x nsite) between the precipitation sites
  #cur.jitter = TRUE/FALSE indicating whether we randomly perturb the non-exceedance probabilities at each site conditional on the basin average value?
  #cur.tc.min = a step change to apply to all tmin temperatures
  #cur.tc.max = a step change to apply to all tmax temperatures
  #num.iter = the number of iterations for the simulation
  #thshd.prcp = the threshold used to separate gamma from GPD
  #qq.month = the frequency prcp under threshold by month and site
  
  #create time series of gamma parameters for both original and CC-scaled gamma distribution
  emission.fit.old.new1 <- CC.scale(emission.fits.site=emission.fits.site,months=months,dates.sim=dates.sim,
                                      n.sites=n.sites,perc.mu=perc.mu,perc.q=perc.q,thshd.prcp=thshd.prcp,qq.month=qq.month)
  emission.old1 <- emission.fit.old.new1[[1]]
  emission.new1 <- emission.fit.old.new1[[2]]
  emission.old2 <- emission.fits.site[[2]]
  #loop through each of the iterations and impose the appropriate scaling change
  for (j in 1:num.iter) {  
    my.systime <- Sys.time()
    #perturb the simulated time series of precipitation at each site
    prcp.site.sim[[j]] <- quantile.mapping(prcp.site=prcp.site.sim[[j]],
                                           Sbasin=Sbasin,thshd.prcp=thshd.prcp,perc.q=perc.q,
                                           emission.old1=emission.old1,
                                           emission.old2=emission.old2,
                                           emission.new1=emission.new1,
                                           n.sites=n.sites,months=months,dates.sim=dates.sim,cur.jitter=cur.jitter)
    
    #Add in temperature trends
    tmin.site.sim[[j]] <- tmin.site.sim[[j]] + cur.tc.min
    tmax.site.sim[[j]] <- tmax.site.sim[[j]] + cur.tc.max
    #dont let tmin exceed tmax - replace tmax with tmin if this situation occurs
    tmax.site.sim[[j]][tmin.site.sim[[j]] > tmax.site.sim[[j]]] <- tmin.site.sim[[j]][tmin.site.sim[[j]] > tmax.site.sim[[j]]]
    
    print(paste("Quant map",j,":", Sys.time()-my.systime))
  }
  return(list(prcp.site.sim,tmin.site.sim,tmax.site.sim))
}