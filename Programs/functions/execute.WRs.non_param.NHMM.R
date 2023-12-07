
execute.WRs.non_param.NHMM <- function(){
  
  
  dir.create(file.path(dir.to.sim.WRs.files), showWarnings = FALSE)
  n.seasons <- length(seasons)
  
  
  ###### output path and filenames ####
  filename_path <- paste0(dir.to.sim.WRs.files,'/final.NHMM.non_param.output.user.developed.rds')
  
  #########################################################################
  dates.synoptic <- seq(as.Date(start.date.synoptic),as.Date(end.date.synoptic), by="days")
  months.synoptic <- as.numeric(format(dates.synoptic,'%m'))
  
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
  
  saveRDS(fit.mod.NHMM,file = paste0(dir.to.sim.WRs.files,"/fit.mod.NHMM.user.developed.rds"))
  
  # fit.mod.NHMM <- readRDS(paste0(dir.to.sim.WRs.files,"/fit.mod.NHMM.user.developed.rds")) # load in this to save running time for now
  # fit.mod.NHMM <- readRDS(paste0(dir.to.sim.WRs.files,"/fit.mod.NHMM.rds")) # load in this to save running time for now
  
  
  #here we glue together the different seasons
  #for each season, put in the final simulations into the right places for the historical
  #re-relabel the next season WRs for the final sequence of WRs
  #dates.sim <- as.Date(sort(unlist(dates.sim.s)),origin="1970-01-01")
  WR.historical.s <- list()
  dates.historical <- as.Date(sort(unlist(dates.final)),origin="1970-01-01")
  WR.historical <- array(NA,length(dates.historical))
  for (s in 1:n.seasons) {
    WR.historical.s[[s]] <- fit.mod.NHMM[[s]][[2]]
    if (s>1) {to.add <- sum(num_WRs.season[1:(s-1)])} else {to.add <- 0}     #relabel the WRs, based on WRs from previous seasons
    WR.historical[dates.historical%in%dates.final[[s]]] <- WR.historical.s[[s]]+to.add
  }
  #
  
  
  
  weather.state.assignments <- WR.historical
  # given the length of non-parametric simulating segments, below should make sense in terms of leap years (every 4 years)
  # e.g., segment for 4-yr; 1948 ('1948') is a leap year; 2020 ('2019+1') is a leap year
  # dates.WRs.specific
  weather.nonpar.state.assignments <- WR.historical[dates.synoptic%in%dates.WRs.specific]
  
  dates.weather.nonpar <- dates.synoptic[dates.synoptic%in%dates.WRs.specific]
  months.weather.nonpar <- as.numeric(format(dates.weather.nonpar,'%m'))
  years.weather.nonpar <- as.numeric(format(dates.weather.nonpar,'%Y'))
  
  my.num.sim = ceiling(num.years.sim.WRs/length(unique(format(dates.weather.nonpar,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
  ##################################################
  
  
  ##################################################
  
  years.vec <- unique(years.weather.nonpar)
  nyears <- length(years.vec)
  n.segment <- 4 # length/number of segments for this case
  my.n.segment <- base::split(1:nyears,ceiling(seq(1,nyears)/n.segment))
  num.segment <- length(my.n.segment)
  size.trace <- nyears*my.num.sim # how many years of simulated data at the end
  total.seqments <- ceiling(size.trace/n.segment)
  
  
  ###get avg frequencies of WRs in each historical block ###
  num_periods <- length(my.n.segment)
  num_WRs <- length(unique(weather.nonpar.state.assignments))
  WR_prob_chunk <- array(NA,c(num_WRs,num_periods))
  for (p in 1:num_periods) {
    chunk.years <- years.vec[my.n.segment[[p]]]
    idx.chunks <- which(years.weather.nonpar%in%chunk.years)
    WR_chunk <- weather.nonpar.state.assignments[idx.chunks]
    WR_prob_chunk[,p] <- table(WR_chunk)/length(WR_chunk)
  }
  
  #####################################
  
  
  
  
  #####################################
  
  output <- WR_prob_estimator(WR_prob_chunk,WR_prob_change,lp.threshold,piecewise_limit,
                              weather.nonpar.state.assignments)
  optim_status <- output[[1]]
  sampling_probabilities <- output[[2]]
  deviation_from_WR_target <- output[[3]]
  deviation_from_prob_unif <- output[[4]]
  deviation_from_prob_unif_piecewise <- output[[5]]
  objval <- output[[6]]
  
  # print(sampling_probabilities)
  
  
  ###Verify that the actual WR probabilities from sampling will match the target
  #calculate target probabilities for WRs
  WR_prob_avg <- table(weather.nonpar.state.assignments)/length(weather.nonpar.state.assignments)
  target_probabilities <- WR_prob_avg*(1+WR_prob_change)
  target_probabilities <- target_probabilities/sum(target_probabilities)   #need to rescale to sum to 1
  
  #expected probabilities from sampling with sampling probabilities for periods
  expected_probabilities <- WR_prob_chunk%*%sampling_probabilities
  
  ###############################################################
  
  
  
  
  ###############################################################
  my.prob <- sampling_probabilities
  set.seed(n.segment)
  my.replication <- sample(seq(1,num.segment),total.seqments,
                           replace = TRUE,
                           prob=my.prob)
  
  dates.sim <- c()
  markov.chain.sim <- c()
  for (p in my.replication){
    
    chunk.years <- years.vec[my.n.segment[[p]]]
    idx.chunks <- which(years.weather.nonpar%in%chunk.years)
    
    dates.sim <- append(dates.sim,
                        dates.weather.nonpar[idx.chunks])
    markov.chain.sim <- append(markov.chain.sim,
                               weather.nonpar.state.assignments[idx.chunks])
  }
  markov.chain.sim <- rep(list(markov.chain.sim),times=num.iter)
  n.length <- length(dates.sim)
  
  
  #final output to be saved and loaded into config.simulations:
  #dates.historical
  #WR.historical
  #dates.sim
  #WR.simulation
  #filename_path
  
  ################################################################################
  final.NHMM.output <- list(dates.historical,
                            weather.state.assignments,
                            dates.sim,
                            markov.chain.sim,
                            filename_path)
  names(final.NHMM.output) <- c('dates.historical','WR.historical',
                                'dates.sim','WR.simulation',
                                'filename_path')
  
  
  saveRDS(final.NHMM.output,filename_path)
  
  gc()
  print(paste0("-->> Saved at: ", filename_path))
  
  return(final.NHMM.output)
  
}

#####################################################################################