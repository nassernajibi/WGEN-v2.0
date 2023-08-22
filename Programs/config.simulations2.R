config.simulations <- function(){
  
  ######  ------------------------
  ######  ------------------------
  ######----> General Settings <-------------------------------------------------######
  ######  ------------------------
  ######  ------------------------
  {
    basin.cnt <- 'Tuolumne.River' # for a set of 12 randomly selected Livneh grids in the Tuolumne River Basin
    # basin.cnt <- 'SFB' # for San Francisco Bay: HUC4-1805 (not yet; will run if we get some extra time)
    
    
    ##adjust main directory and directory for simulation files
    mainDir <- "D:/Projects/Tuolumne_River_Basin/GitHub_WGENv2.0"
    setwd(mainDir)
    
    dir.to.sim.files <- "./Data/simulated.data.files/WGEN.out"
    dir.create(file.path(mainDir, dir.to.sim.files), showWarnings = FALSE)
    
    ##location of obs weather data (RData format): weather data (e.g., precip and temp) as matrices (time x lat|lon: t-by-number of grids); dates vector for time; basin average precip (see the example meteohydro file)
    path.to.processed.data.meteohydro <- paste0("./Data/processed.data.files/processed.meteohydro/processed.meteohydro.",basin.cnt,".RData")
    months <- seq(1,12) # Jun-Dec calendar year
    
    
    ##threshold for mixed Gamma-GPD population separation
    qq <- .99  
    
    
    ##bootstrapping choices##
    window.size <- rep(3,length(months))   #the size of the window (in days) from which runs can be bootstrapped around the current day of simulation, by month: Jan -- Dec
    pr.trace <- 0.25     # {0.25 mm, 0.01 in} trace prcp threshold. 0.25 mm (for Livneh dataset); or 0.01 inches: lower threshold below which day is considered dry
    
    ##-------------Define perturbations-------------##
    ##climate changes and jitter to apply:
    change.list <- data.frame("tc"=  c(0), # {e.g., }0, 1, 2, ...} (changes in temperature)
                              "jitter"= c(TRUE),
                              "pccc"= c( 0), # {e.g., 0, 0.07, 0.14, ...} (changes for precipitation extreme quantile -- CC)
                              "pmuc"= c( 0)# {e.g., 0, -.125, .125, ...} (changes in precipitation mean)
    )
    ##----------------------------------------------##
    
    ##load in supporting functions
    files.sources = list.files("./Programs/functions",full.names = TRUE)
    my.functions <- sapply(files.sources, source)
    
  }
  
  ######  -----------------------
  ######  -----------------------
  ######----> Dates and Years <------------------------------------------------------------######
  ######  -----------------------
  ######  -----------------------
  {
    ##length of final simulated weather (in calendar years)##
    number.years.long <- 1000 # {e.g., 500, 1000, 2000, 3000, 5000 years,...} [note: current NHMM output (parametric) is for 1036 years; current non-parametric is for 3050 years]
    num.iter <- 1 # A single long trace (e.g., thousand years) is sufficient although we create like 5 ensembles in the simulated WRs
    
    start.date.weather="1950-01-01"; end.date.weather="2013-12-31" # Toulamne River
    # start.date.weather="1948-01-01"; end.date.weather="2018-12-31" # SFB
    
    start.date.synoptic="1948-01-01"; end.date.synoptic="2021-12-31" # from processed GPHA file
    
    start.date.par="1948-01-01"; end.date.par="2019-12-31" # proper leap year orders (starting with leap year of 1948, ending a year before (i.e., 2019) the leap year of 2020)
    start.date.nonpar="1948-01-01"; end.date.nonpar="2019-12-31" # proper leap year orders (starting with leap year of 1948, ending a year before (i.e., 2019) the leap year of 2020)
  }
  
  ######  -----------------------------------------------------------------
  ######  -----------------------------------------------------------------
  ######----> Hyperparameters of the WRs Identification and Simulation <---------######
  ######  -----------------------------------------------------------------
  ######  -----------------------------------------------------------------
  {
    num.years.sim.WRs <- 1000 # e.g., 500, 1000, 2000, 3000, 5000 years, etc [note: current NHMM output (parametric) is for 1036 years]
    
    dir.to.sim.WRs.files <- "./Data/simulated.data.files/WRs.out" # dir.to.sim.WRs.files
    num.iter.WRs <- 1   #number of iterations to simulate sequence of WRs
    path.to.processed.GPHAs <- './Data/processed.data.files/processed.hgt/hgt.500.Pacific.NorthAmer.synoptic.region_19480101_20211231.rds'
    # Covariates should be a matrix with the first column as dates, and the second column as 
    #      ... normalized pPC1 (scaled and centered)
    path.to.processed.SPI.PCs <- './Data/processed.data.files/processed.NHMM.data/paleo.norm.4.cold.PCs.dates_extracted.rds'
    
    #######################define seasons and covariates for NHMM models of WRs#############################
    cold.months <- c(11,12,1,2,3,4) # Nov-Apr
    warm.months <- c(5,6,7,8,9,10)  #May-Oct
    num.PCs <- 10
    seasons <- list(cold.months,warm.months)
    num_eofs.season <- rep(num.PCs,length(seasons))  #number of PCs to use for geopotential heights per season
    num_WRs.season <- c(7,3)    #number of WRs to fit per season
    
    
    ##Choose below whether through parametric or non-parametric way to create the simulated WRs ##
    use.non_param.WRs <- TRUE #{TRUE, FALSE}: TRUE for non-parametric, FALSE for parametric simulated WRs
    
    dynamic.scenario  <- 0 # {0, 1, 2}: 0: no dynamic change; 1: dynamic scenario #1 (30% increase in WR3); or 2: dynamic scenario #2 (linear trend)
    
    if (use.non_param.WRs){      #----------- 1+2 dynamic scenarios ----------#
      if (dynamic.scenario==0){
        ##===> Attempt #0 (thermodynamic only; no change to freq of WRs) ===##
        # #specify target change (as a percent) for WR probabilities
        WR_prob_change <- c(0,0,0,0,0,0,0,0,0,0) # between 0 and 1
        # #how close (in % points) do the WR frequencies (probabilities) need to be to the target
        lp.threshold <- 0.00001
        # #how much change do we allow in a sub-period sampling probability before incurring a larger penalty in the optimization
        piecewise_limit <- .02
        
      }else if(dynamic.scenario==1){
        ##===> Attempt #1 (dynamic scenario #1) ===##
        # #specify target change (as a percent) for WR probabilities (if, increasing WR3 in future)
        WR_prob_change <- c(0,0,.3,0,0,0,0,0,0,0) # between 0 and 1
        # #how close (in % points) do the WR frequencies (probabilities) need to be to the target
        lp.threshold <- 0.007
        # #how much change do we allow in a sub-period sampling probability before incurring a larger penalty in the optimization
        piecewise_limit <- .02
        
      }else if(dynamic.scenario==2){
        ##===> Attempt #2 (dynamic scenario #2) ===##
        # specify target change (as a percent) for WR probabilities (if, continuing their current trends in future)
        WR_prob_change <- c(-0.09969436,  0.27467048,  0.33848792,
                            -0.28431861, -0.23549986,  0.03889970,
                            -0.05628958, 0.38059153, -0.16636739, -0.17995965) # between 0 and 1
        # how close (in % points) do the WR frequencies (probabilities) need to be to the target
        lp.threshold <- 0.008
        # how much change do we allow in a sub-period sampling probability before incurring a larger penalty in the optimization
        piecewise_limit <- .02
      }
    }
    
  }
  
  ######  -----------------------------------
  ######  -----------------------------------
  ######----> Create Figures <---------######
  ######  -----------------------------------
  ######  -----------------------------------
  {
    # labels for x and y-axes
    label1 <- paste0('Obs [',format(as.Date(start.date.weather),'%Y'),
                     '-',format(as.Date(end.date.weather),'%Y'),']')
    label2 <- 'Sim (WGEN: baseline)'
  }
  
  ######  ----------------------------------------
  ######  ----------------------------------------
  ######----> Create Output files <---------######
  ######  ----------------------------------------
  ######  ----------------------------------------
  {
    # directory to store output files
    dir.to.output.files <- './Data/output.data.files/'
  }
  
  # returning the entire values inserted here to 'run.stochastic.weather.generator'
  values = as.list(environment())
  return(values)
}