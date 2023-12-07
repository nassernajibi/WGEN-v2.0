execute.simulations <- function(){
  
  #########Precipitation/Temperature Inputs #########
  
  #location of obs weather data (RData format): weather data (e.g., precip and temp) as matrices (time x lat|lon: t-by-number of grids); dates vector for time; basin average precip (see the example meteohydro file)
  load(path.to.processed.data.meteohydro) #load in weather data
  n.sites <- dim(prcp.site)[2] # Number of gridded points for precipitation
  
  identical.dates.idx <- dates.weather%in%dates.synoptics
  dates.weather <- dates.weather[identical.dates.idx]
  months.weather <- as.numeric(format(dates.weather,'%m'))
  prcp.site <- prcp.site[identical.dates.idx,]
  tmax.site <- tmax.site[identical.dates.idx,]
  tmin.site <- tmin.site[identical.dates.idx,]
  prcp.basin <- prcp.basin[identical.dates.idx]
  
  identical.dates.idx <- dates.synoptics%in%dates.weather
  weather.state.assignments <- weather.state.assignments[identical.dates.idx]
  
  # sanitary check if all days for that specific month and grid is zero === causing problems in gamma fit
  # this will be zero after implementing the quantile mapping, so there is no effect whatsoever
  {
    indicators.month.sites <- array(NA,c(length(months),n.sites))
    indicators.month.unit <- array(NA,length(months))
    for (mo in months){
      for (r in 1:n.sites){
        indicators.month.sites[mo,r] <- sum(prcp.site[months.weather==mo,r]==0)/sum(months.weather==mo)
      }
      indicators.month.unit[mo] <- sum(months.weather==mo)
    }
    required.daily.vals <- 2
    low.sample.size <- min(1-required.daily.vals/indicators.month.unit)
    if (sum(indicators.month.sites>=low.sample.size)>0){
      mo.site <- which(indicators.month.sites>=low.sample.size,arr.ind = T)
      for (mi in 1:dim(mo.site)[1]){
        set.seed(100)
        escaping.value <- runif(sum(months.weather==mo.site[mi,1]),min=0,max=0.01) # giving random small values in mm
        prcp.site[months.weather==mo.site[mi,1],mo.site[mi,2]] <- escaping.value
      }
    }
    
    # check if all days for that specific month and grid is negative === causing problems in gamma fittings
    indicators.sites <- which(prcp.site<0,arr.ind=T)
    if (dim(indicators.sites)[1]>0){
      set.seed(1000)
      escaping.value <- runif(dim(indicators.sites)[1],min=0,max=0.001) # giving random values between 0 and 0.001 mm
      prcp.site[indicators.sites] <- round(escaping.value,4)
    }
  }

  #extreme threshold (Gamma-GPD cutoff) for each site
  thshd.prcp <- apply(prcp.site,2,function(x) {quantile(x[x!=0],qq,na.rm=T)})
  #The spearman correlation between the precipitation sites, used in the copula-based jitters
  Sbasin <- cor(prcp.site,method="spearman")
  
  
  ######### Gamma-GPD Distribution Parameters #########
  
  #fit emission distributions to each site by month. sites along the columns, parameters for month down the rows
  emission.fits.site <- fit.emission(prcp.site=prcp.site,
                                     months=months,
                                     months.weather=months.weather,
                                     n.sites=n.sites,
                                     thshd.prcp=thshd.prcp)
  
  
  #how often is prcp under threshold by month and site
  qq.month <- sapply(1:n.sites,function(i,x=prcp.site,m,mm) {
    sapply(m,function(m) {
      length(which(as.numeric(x[x[,i]!=0 & mm==m & !is.na(x[,i]),i])<=thshd.prcp[i]))/length(as.numeric(x[x[,i]!=0 & mm==m & !is.na(x[,i]),i]))
    })},
    m=months, mm=months.weather)
  
  #################################################
  
  ##########################Simulate model with perturbations#######################################

  #run the daily weather generate num.iter times using the num.iter Markov chains
  mc.sim <- resampled.date.sim <- resampled.date.loc.sim <- prcp.site.sim <- tmin.site.sim <- tmax.site.sim <-  list()
  start_time <- Sys.time()
  for (k in 1:num.iter) {
    my.itertime <- Sys.time()
    my.sim <- wgen.simulator(weather.state.assignments=weather.state.assignments,mc.sim=markov.chain.sim[[k]],
                             prcp.basin=prcp.basin,dates.weather=dates.weather,
                             first.month=first.month,last.month=last.month,dates.sim=dates.sim,
                             months=months,window.size=window.size,min.thresh=pr.trace)
    print(paste(k,":", Sys.time()-my.itertime))
    
    #each of these is a list of length iter
    mc.sim[[k]] <- my.sim[[1]]
    resampled.date.sim[[k]] <- my.sim[[2]]
    resampled.date.loc.sim[[k]] <- my.sim[[3]]
    prcp.site.sim[[k]] <- prcp.site[resampled.date.loc.sim[[k]],]
    tmin.site.sim[[k]] <- tmin.site[resampled.date.loc.sim[[k]],]
    tmax.site.sim[[k]] <- tmax.site[resampled.date.loc.sim[[k]],]
  }
  end_time <- Sys.time(); run.time.sim <- end_time - start_time
  print(paste("SIMULATION started at:",start_time,", ended at:",end_time)); print(round(run.time.sim,2))
  #remove for memory
  rm(my.sim)
  
  
  #once the simulations are created, we now apply post-process climate changes (and jitters)
  for (change in 1:nrow(change.list)) {
    start_time <- Sys.time()
    cur.tc.max <- change.list$tc.max[change]
    cur.tc.min <- change.list$tc.min[change]
    cur.pccc <- change.list$pccc[change]
    cur.pmuc <- change.list$pmuc[change]
    
    cur.tc <- mean(cur.tc.max,cur.tc.min)
    
    #precipitation scaling (temperature change dependent)
    perc.q <- (1 + cur.pccc)^cur.tc    #scaling in the upper tail for each month of non-zero prcp
    perc.mu <- (1 + cur.pmuc)          #scaling in the mean for each month of non-zero prcp
    
    #set the jitter
    cur.jitter <- to.jitter
    
    #perturb the climate from the simulations above (the longest procedure in this function is saving the output files)
    set.seed(1)   #this ensures the copula-based jitterrs are always performed in the same way for each climate change
    perturbed.sim <-  perturb.climate(prcp.site.sim=prcp.site.sim,
                                      tmin.site.sim=tmin.site.sim,
                                      tmax.site.sim=tmax.site.sim,
                                      emission.fits.site=emission.fits.site,
                                      months=months,dates.sim=dates.sim,n.sites=n.sites,
                                      qq=qq,perc.mu=perc.mu,perc.q=perc.q,Sbasin=Sbasin,cur.jitter=cur.jitter,
                                      cur.tc.min=cur.tc.min,cur.tc.max=cur.tc.max,
                                      num.iter=num.iter,thshd.prcp=thshd.prcp,qq.month=qq.month)
    set.seed(NULL)
    prcp.site.sim.perturbed <- perturbed.sim[[1]]
    tmin.site.sim.perturbed <- perturbed.sim[[2]]
    tmax.site.sim.perturbed <- perturbed.sim[[3]]
    
    #remove for memory
    rm(perturbed.sim)
    end_time <- Sys.time(); run.time.qmap <- end_time - start_time
    print(paste("QMAPPING started at:",start_time,", ended at:",end_time)); print(round(run.time.qmap,2))
    
    #how to name each file name to track perturbations in each set of simulations
    file.suffix <- paste0(".tmax.",cur.tc.max,".tmin.",cur.tc.min,"_p.CC.scale.",cur.pccc,"_p.mu.scale.",cur.pmuc,"_num.year.",number.years.long,"_with.",num.iter)
    
    print(paste("|---start saving---|"))
    write.output.large(dir.to.sim.files,
                       prcp.site.sim=prcp.site.sim.perturbed,
                       tmin.site.sim=tmin.site.sim.perturbed,
                       tmax.site.sim=tmax.site.sim.perturbed,
                       mc.sim=mc.sim,resampled.date.sim=resampled.date.sim,
                       dates.sim=dates.sim,file.suffix=file.suffix)
    print(paste("|---finished saving---|"))
    
  }
  
  #remove for memory
  rm(resampled.date.sim,resampled.date.loc.sim,markov.chain.sim,mc.sim)
  ##################################################################################################
  print(paste0("--- done.  state= ", num.states," --- ensemble member:",num.iter))
  gc()
  print(paste0("------------------------------------------------------"))
  print(paste0("-->> simulated files were saved at= ", dir.to.sim.files))
  

}
# The End
#####################################################################################