
create.delimited.outputs <- function(scenario = selected_scenario){
  
  # this function creates output files #
  #scenario = the row in ClimateChangeScenarios.csv for which to plot results
  
  ##// weather data and synoptic dates ---##
  #location of obs weather data
  load(path.to.processed.data.meteohydro) #load in weather data
  dates.user.specific <- seq(as.Date(start.date.weather),as.Date(end.date.weather),by="day")
  identical.dates.idx <- dates.weather%in%dates.user.specific
  
  prcp.site.data <- prcp.site[identical.dates.idx,]
  tmin.site.data <- tmin.site[identical.dates.idx,]
  tmax.site.data <- tmax.site[identical.dates.idx,]
  
  dates.weather <- dates.user.specific
  years.weather <- as.numeric(format(dates.weather,'%Y'))
  months.weather <- as.numeric(format(dates.weather,'%m'))
  days.weather <- as.numeric(format(dates.weather,'%d'))
  
  list.file.names <- gsub("data", "sim", list.file.names)

  n.sites <- dim(prcp.site)[2] # Number of gridded points for precipitation
  cur.jitter <- to.jitter      #whether jitters were activated or not
  
  # Sim. file
  cur.tc <- change.list$tc[scenario]
  cur.pccc <- change.list$pccc[scenario]
  cur.pmuc <- change.list$pmuc[scenario]
  
  ##// sim files ---##
  {
    simulated.file.run.model.saved <- paste0(".temp.",cur.tc,"_p.CC.scale.",cur.pccc,"_p.mu.scale.",cur.pmuc,"_hist.state.",use.non_param.WRs,"_jitter.",cur.jitter,"_s",num.states,"_with_",num.iter)
    
    prcp.site.sim_sfx <- "prcp.site.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    prcp.site.data.sim <- prcp.site.sim[[scenario]] # because there was only one iteration
    
    prcp.site.sim_sfx <- "tmin.site.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    tmin.site.data.sim <- tmin.site.sim[[scenario]]
    
    prcp.site.sim_sfx <- "tmax.site.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    tmax.site.data.sim <- tmax.site.sim[[scenario]]
    
    tmean.site.data.sim <- (tmax.site.data.sim+tmin.site.data.sim)/2
    
    prcp.site.sim_sfx <- "dates.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    #dates.sim
    
    months.sim <- as.numeric(format(dates.sim,"%m"))
    years.sim <- as.numeric(format(dates.sim,"%Y"))
    days.sim <- as.numeric(format(dates.sim,"%d"))
    
  }
  
  {
    # generating synthetic dates/years #
    a <- rle(years.sim)
    a.freq <- a$lengths
    a.val <- a$values
    yr=1
    yr.vec <- c()
    for (k in 1:length(a.freq)){
      yr.vec <- append(yr.vec,rep(yr,a.freq[k]))
      yr <- yr+1
    }
    long.years.vec <- yr.vec
    long.months.vec <- as.numeric(months.sim)
    total.num.yrs <- max(long.years.vec)
    
  }
  
  
  ##/ create a single long trace for each site
  
  long.dates.sim <- cbind(long.years.vec,long.months.vec,days.sim)
  long.dates.obs <- cbind(years.weather,months.weather,days.weather)
  
  start_time <- Sys.time()
  
  dir.create(file.path(dir.to.output.files), showWarnings = FALSE)
  
  for (k in 1:n.sites){
    
    # sim output
    my.filename <- paste0(list.file.names[k])
    dates.vec <- long.dates.sim
    prcp.atemp <- cbind(prcp.site.data.sim[,k],
                        tmax.site.data.sim[,k],
                        tmin.site.data.sim[,k])
    mat.to.save <- as.matrix(cbind(dates.vec,prcp.atemp)) # YEAR MONTH DAY PRCP MAX-TEMP MIN-TEMP
    write.table(mat.to.save,file=paste0(dir.to.output.files,my.filename),
                row.names = FALSE, col.names = FALSE, sep = " , ",quote = FALSE)
    
    print(paste('site',k,'out of',n.sites,
                '--',round(k/n.sites*100,2),"%"))
  }
  
  end_time <- Sys.time(); run.time.saving.sim <- end_time - start_time
  print(round(run.time.saving.sim,2))
  
  
}
