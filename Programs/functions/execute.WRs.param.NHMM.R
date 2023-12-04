
execute.WRs.param.NHMM <- function(){
  
  
  dir.create(file.path(dir.to.sim.WRs.files), showWarnings = FALSE)
  n.seasons <- length(seasons)
  

  #########################################################################
  
  #dates for the historical WRs
  dates.weather <- seq(as.Date(start.date.weather),as.Date(end.date.weather),by="day") # also known as dates.shared amongst
  
  # below should make sense in terms of leap years (starting with leap year (1948); ending to a year before leap year in history)
  dates.shared.par <- seq(as.Date(start.date.par),as.Date(end.date.par),by="day")
  
  dates.synoptic <- seq(as.Date(start.date.synoptic),as.Date(end.date.synoptic), by="days")
  months.synoptic <- as.numeric(format(dates.synoptic,'%m'))
  
  identical.dates.idx <- dates.synoptic%in%dates.shared.par # indexes of dates to replicate different segments considering leap years
  
  #create dates for the simulated WRs
  my.num.sim = ceiling(num.years.sim.WRs/length(unique(format(dates.synoptic[identical.dates.idx],'%Y')))) # number of chunks of historical periods
  long.dates.sim <- rep(dates.synoptic[identical.dates.idx],times=my.num.sim)
  target.long.dates <- list(long.dates.sim[as.numeric(format(long.dates.sim[identical.dates.idx],"%m"))%in%seasons[[1]]],
                            long.dates.sim[as.numeric(format(long.dates.sim[identical.dates.idx],"%m"))%in%seasons[[2]]])
  
  
  #// Load processed 500mb-GPH hgt region data and dates
  hgt.synoptic.region <- readRDS(file=path.to.processed.GPHAs)
  
  
  #####Load in covariates for the cold season
  #Covariates should be a matrix with the first column as dates, and the second column as 
  # normalized pPC1 (scaled and centered)
  SPI.dates.PCA.extracted <- readRDS(file=path.to.processed.SPI.PCs)
  n.cov <- nrow(SPI.dates.PCA.extracted)
  k.cov <- ncol(SPI.dates.PCA.extracted)
  dates.covariates.cold <- SPI.dates.PCA.extracted$WY_dates.added
  covariates.cold <- SPI.dates.PCA.extracted[,2:k.cov]
  names(covariates.cold) <- paste("X",1:(k.cov-1),sep="")
  ####
  
  
  ####Load in covariates for the warm season (if any)###
  covariates.warm <- NULL
  dates.covariates.warm <- NULL
  
  covariates.all <- list(covariates.cold,covariates.warm)
  dates.covariates.all <- list(dates.covariates.cold,dates.covariates.warm)
  # indexes of dates to replicate different segments considering leap years:
  identical.cov.dates.idx <- dates.covariates.cold%in%dates.shared.par
  
  ############################################################################
  
  
  ############Load in geopotential heights and separate by season##############
  
  covariates.final <- list()
  dates.final <- list()
  hgt.final <- list()
  for (s in 1:n.seasons) {
    covariates.cur <- covariates.all[[s]]
    dates.covariates.cur <- dates.covariates.all[[s]]
    
    if (is.null(covariates.cur)) {
      dates.final[[s]] <- dates.synoptic[months.synoptic%in%seasons[[s]]]
      hgt.final[[s]] <- hgt.synoptic.region[months.synoptic%in%seasons[[s]],]
      covariates.final[[s]] <- NA
    } else {
      months.covariates.cur <- as.numeric(format(dates.covariates.cur,"%m"))
      dates.final[[s]] <- dates.covariates.cur[months.covariates.cur%in%seasons[[s]] & dates.covariates.cur%in%dates.synoptic]
      hgt.final[[s]] <- hgt.synoptic.region[dates.synoptic%in%dates.final[[s]],]
      covariates.final[[s]] <- covariates.cur[dates.covariates.cur%in%dates.final[[s]],]
    }
  }
  
  #####################################################################################
  
  
  ##################fit NHMMs by season with appropriate covariates#####################
  
  fit.mod.NHMM <- list()    #final NHMMs fit by season
  
  for (s in 1:n.seasons) {  
    #organize hgt and covariate for current season
    hgt.synoptic.region.pca <- prcomp(hgt.final[[s]],center=T,scale=T)
    n.eofs <- num_eofs.season[s]
    synoptic.pcs <- hgt.synoptic.region.pca$x[,1:n.eofs]
    n.states <- num_WRs.season[s]
    
    #here we define the covariate matrix and the formula (fo) of the multinomial regression in the NHMM
    #define day of year for all dates, and then define formula for the NHMM
    tsteps <- as.numeric(format(dates.final[[s]], "%j")) #seq_len(nrow(my.synoptic.pcs))
    omegaT <- (2*pi*tsteps)/365
    if (is.na(all(covariates.final[[s]]))) {
      my.covar.trans <- data.frame(intercept=1, CosT = cos(omegaT), SinT = sin(omegaT))
      fo <- "~ -1 + intercept + CosT + SinT"     #if we want to drop intercept: "~ -1 + CosT + SinT"
      fo <- as.formula(fo)  
    } else{
      my.covar.trans <- data.frame(intercept=1, CosT = cos(omegaT), SinT = sin(omegaT),covariates.final[[s]])
      num.cov <- ncol(covariates.final[[s]])
      fo <- "~ -1 + intercept + CosT + SinT"
      for (h in 1:num.cov) {
        fo <- paste(fo," + X",h,sep="")
      }
      fo <- as.formula(fo)     
    }
    
    
    # fit NHMM. this can take a while (tries multiple random starts until convergence)
    fit.mod.NHMM[[s]] <- fit.NHMM(my.nstates=n.states,
                                  my.synoptic.pcs=synoptic.pcs,
                                  my.covar.trans=my.covar.trans,
                                  fo=fo,
                                  n.eofs=n.eofs)
  }
  
  saveRDS(fit.mod.NHMM,file = paste0(dir.to.sim.WRs.files,"/fit.mod.NHMM.rds"))
  #fit.mod.NHMM <- readRDS(paste0(dir.to.sim.WRs.files,"/fit.mod.NHMM.rds"))
  #####################################################################################
  
  
  
  ##################Provide arbitrary dates for any length##############################

  s=1 # for cold season only in which we have paleo PCs
  u <- covariates.final[[s]]
  u.limited <- u[identical.cov.dates.idx,]
  uRepped <- u.limited[rep(seq_len(nrow(u.limited)), my.num.sim),]
  covariates.final.long <- list(uRepped,'NA')
  
  
  ##################Simulate from the NHMMs to get simulated WRs########################
  
  sim.mod.NHMM <- list()    #simulations from NHMMs across iterations
  WR.historical.s <- list()
  WR.simulation.s <- list()
  dates.sim.s <- list()
  
  for (s in 1:n.seasons) {
    
    n.states <- num_WRs.season[s]
    
    #######here define the covariate matrix to use for the simulation#########
    
    #in this simple case, we use the same covariates as used for fitting . but they could be of different length or new values.
    dates.sim.s[[s]] <- target.long.dates[[s]]
    covariates.sim.s <- covariates.final.long[[s]]
    
    tsteps <- as.numeric(format(dates.sim.s[[s]], "%j")) 
    omegaT <- (2*pi*tsteps)/365
    if (is.na(all(covariates.sim.s))) {
      my.covar.trans.sim <- data.frame(intercept=1, CosT = cos(omegaT), SinT = sin(omegaT))
    } else{
      my.covar.trans.sim <- data.frame(intercept=1, CosT = cos(omegaT), SinT = sin(omegaT),
                                       covariates.sim.s)
    }
    
    
    #create the time varying transition probabilities based on the covariates, stored in an mc object. there will be n.sim different tpm's
    mcObject.final <- create.mc.object(my.NHMM=fit.mod.NHMM[[s]][[1]],my.nstates=n.states,my.covar.trans.sim=my.covar.trans.sim)
    #define historical WRs
    WR.historical.s[[s]] <- fit.mod.NHMM[[s]][[2]]
    #use that mc object to simulate n.sim WRs
    WR.simulation <- array(NA,c(length(mcObject.final@markovchains),num.iter))
    for (k in 1:num.iter) {
      WR.simulation[,k] <- simulate.NHMM(mcObject.final=mcObject.final)
    }
    WR.simulation.s[[s]] <- WR.simulation
    
  }
  
  
  #here we glue together the different seasons
  #for each season, put in the final simulations into the right places for each iteration.
  #do the same thing for the historical
  #re-relabel the next season WRs for the final sequence of WRs
  #dates.sim <- as.Date(sort(unlist(dates.sim.s)),origin="1970-01-01")
  dates.sim <- long.dates.sim
  WR.simulation <- array(NA,c(length(dates.sim),num.iter))
  dates.historical <- as.Date(sort(unlist(dates.final)),origin="1970-01-01")
  WR.historical <- array(NA,length(dates.historical))
  for (s in 1:n.seasons) {
    if (s>1) {to.add <- sum(num_WRs.season[1:(s-1)])} else {to.add <- 0}     #relabel the WRs, based on WRs from previous seasons
    for (k in 1:num.iter) {
      WR.simulation[dates.sim%in%dates.sim.s[[s]],k] <- WR.simulation.s[[s]][,k]+to.add
    }
    WR.historical[dates.historical%in%dates.final[[s]]] <- WR.historical.s[[s]]+to.add
  }
  
  
  
  #final output to be saved and loaded into config.simulations:
  #dates.historical
  #WR.historical
  #dates.sim
  #WR.simulation
  #filename_path
  filename_path <- paste0(dir.to.sim.WRs.files,"/final.NHMM.param.output.rds")
  final.NHMM.output <- list(dates.historical,WR.historical,dates.sim,WR.simulation,filename_path)
  names(final.NHMM.output) <- c('dates.historical','WR.historical','dates.sim','WR.simulation','filename_path')
  saveRDS(final.NHMM.output,file = filename_path)
  
  gc()
  print(paste0("-->> Saved at: ", filename_path))
  
  return(final.NHMM.output)
}

#####################################################################################