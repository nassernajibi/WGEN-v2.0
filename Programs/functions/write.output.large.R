
write.output.large <- function(dir.to.sim.files,
                               prcp.site.sim,
                               tmin.site.sim,
                               tmax.site.sim,
                               mc.sim,resampled.date.sim,dates.sim,file.suffix)
{
  
  #this function writes out the "large" amount of the simulated data as .RData files
  
  #Arguments:
  #prcp.site.sim = a list of length iter, with each list element containing a matrix of simulated daily prcp (rows: days; cols: sites)
  #...
  #mc.sim = a list of length iter, with each list element containing a Markov chain vector
  #resampled.date.sim = a list of length iter, with each list element containing a vector of resampled dates
  #dates.sim = a vector of dates over the simulation period
  #file.suffix = a string to append to the end of each file name to track perturbations in each set of simulations
  
  save(prcp.site.sim,file= paste0(dir.to.sim.files,"/prcp.site.sim",file.suffix,".RData"))          
  save(tmin.site.sim,file= paste0(dir.to.sim.files,"/tmin.site.sim",file.suffix,".RData"))          
  save(tmax.site.sim,file= paste0(dir.to.sim.files,"/tmax.site.sim",file.suffix,".RData"))          
  save(mc.sim,file = paste0(dir.to.sim.files,"/mc.sim",file.suffix,".RData"))
  save(resampled.date.sim,file=  paste0(dir.to.sim.files,"/resampled.date.sim.sim",file.suffix,".RData"))  
  save(dates.sim,file = paste0(dir.to.sim.files,"/dates.sim",file.suffix,".RData"))
  
}