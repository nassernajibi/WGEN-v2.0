
create.delimited.outputs <- function(scenario = selected_scenario){
  
  # this function creates output files #
  #scenario = the row in ClimateChangeScenarios.csv for which to plot results
  
  ##// weather data and synoptic dates ---##
  #location of obs weather data
  load(path.to.processed.data.meteohydro) #load in weather data

  # Sim. file
  cur.tc.max <- change.list$tc.max[scenario]
  cur.tc.min <- change.list$tc.min[scenario]
  cur.pccc <- change.list$pccc[scenario]
  cur.pmuc <- change.list$pmuc[scenario]
  
  ##// sim files ---##
  {
    simulated.file.run.model.saved <- paste0(".tmax.",cur.tc.max,".tmin.",cur.tc.min,"_p.CC.scale.",cur.pccc,"_p.mu.scale.",cur.pmuc,"_num.year.",number.years.long,"_with.",num.iter)
    
    load(paste0(dir.to.sim.files,"/","prcp.site.sim",simulated.file.run.model.saved,".RData"))
    prcp.site.data.sim <- prcp.site.sim[[1]] # because there was only one iteration
    
    load(paste0(dir.to.sim.files,"/","tmin.site.sim",simulated.file.run.model.saved,".RData"))
    tmin.site.data.sim <- tmin.site.sim[[1]]
    
    load(paste0(dir.to.sim.files,"/","tmax.site.sim",simulated.file.run.model.saved,".RData"))
    tmax.site.data.sim <- tmax.site.sim[[1]]
    
    #dates.sim
    load(paste0(dir.to.sim.files,"/","dates.sim",simulated.file.run.model.saved,".RData"))
    
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
