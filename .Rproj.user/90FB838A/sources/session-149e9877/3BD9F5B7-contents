config.simulations <- function(){
  
  ######  --------------------------------------------------------------------------
  ######  --------------------------------------------------------------------------
  ######  ---------------> Simulation Length <------------------------------- ######
  ######  --------------------------------------------------------------------------
  ######  --------------------------------------------------------------------------
  {
    ##length of final simulated weather (in calendar years)##
    simulation.length <- suppressWarnings(read.table('SimulationLength.csv',header=TRUE,sep=","))
    number.years.long <- simulation.length$number_of_years_per_ensemble_member # {e.g., 500, 1000, 2000, 3000, 5000 years,...} [note: current NHMM output (parametric) is for 1036 years; current non-parametric is for 3050 years]
    num.iter <- simulation.length$number_of_ensemble_members # A single long trace (e.g., thousand years) is sufficient, although more can be developed if desired
  }
  
  
  ######  --------------------------------------------------------------------------
  ######  --------------------------------------------------------------------------
  ######----> Thermodynamic Climate Change Scenarios<-------------------------######
  ######  --------------------------------------------------------------------------
  ######  --------------------------------------------------------------------------
  {
    ##-------------Define perturbations-------------##
    ##climate changes and jitter to apply:
    climate.change.scenarios <- suppressWarnings(read.table('ClimateChangeScenarios.csv',header=TRUE,sep=","))
    
    change.list <- data.frame("tc.max"=  climate.change.scenarios$max_temperature_change_degC, # {e.g., 0, 1, 2, ...} (changes in temperature)
                              "tc.min"=  climate.change.scenarios$min_temperature_change_degC, # {e.g., 0, 1, 2, ...} (changes in temperature)
                              "pmuc"= climate.change.scenarios$mean_precipitation_change_percent/100, # {e.g., 0, -.125, .125, ...} (changes in precipitation mean)
                              "pccc"= climate.change.scenarios$extreme_precipitation_scaling_rate_percent/100 # {e.g., 0, 0.07, 0.14, ...} (changes for precipitation extreme quantile -- CC)
    )
    ##----------------------------------------------##
  }
  
  
  ######  --------------------------------------------------------------------------
  ######  --------------------------------------------------------------------------
  ######  ----------------> Dates for Weather Data <---------------------------######
  ######  --------------------------------------------------------------------------
  ######  --------------------------------------------------------------------------
  {
    lst.import.datafile <- tryCatch(suppressMessages(readRDS(file='./Data/processed.data.files/processed.meteorology/lst.import.datafile.rds')),
                                    error=function(e) {
                                      message('You have not yet run process.meteorology.R')
                                      print(e)
                                    })
    start.date.weather <- lst.import.datafile$seq.of.dates[1]
    end.date.weather <- lst.import.datafile$seq.of.dates[length(lst.import.datafile$seq.of.dates)]
  }
  
  ######  ----------------------------------------
  ######  ----------------------------------------
  ######----> Directories <---------######
  ######  ----------------------------------------
  ######  ----------------------------------------
  {
    dir.to.sim.files <- "./Data/simulated.data.files/WGEN.out"
    dir.create(file.path(dir.to.sim.files), showWarnings = FALSE)
    
    # directory to store output files
    dir.to.output.files <- './Data/output.data.files/'
    
    ##location of obs weather data (RData format): weather data (e.g., precip and temp) as matrices (time x lat|lon: t-by-number of grids); dates vector for time; basin average precip (see the example meteohydro file)
    path.to.processed.data.meteohydro <- "./Data/processed.data.files/processed.meteorology/processed.meteorology.RData"
  }
  
  
  ######  -----------------------------------------------------
  ######  -----------------------------------------------------
  ######  -------> WGEN Hyperparameters <--------------- ######
  ######  -----------------------------------------------------
  ######  -----------------------------------------------------
  {
    first.month <- 1
    last.month <- 12
    months <- seq(first.month,last.month) # Jan-Dec calendar year
    ##threshold for mixed Gamma-GPD population separation
    qq <- .99  
    
    #keep the jittering on
    to.jitter <- TRUE
    
    ##bootstrapping choices##
    window.size <- rep(3,length(months))   #the size of the window (in days) from which runs can be bootstrapped around the current day of simulation, by month: Jan -- Dec
    pr.trace <- 0.25     # {0.25 mm, 0.01 in} trace prcp threshold. 0.25 mm (for Livneh dataset); or 0.01 inches: lower threshold below which day is considered dry
    
    ##load in supporting functions
    files.sources = list.files("./Programs/functions",full.names = TRUE)
    my.functions <- sapply(files.sources, source)
  }
  
  ######  -----------------------------------------------------------------------------
  ######  -----------------------------------------------------------------------------
  ######----> Hyperparameters of the WRs Identification and Simulation <---------######
  ######  -----------------------------------------------------------------------------
  ######  -----------------------------------------------------------------------------
  {
    ##Choose below whether to use provided WRs (TRUE), or run WRs identification from scratch (FALSE)
    use.provided.WRs <- TRUE #{TRUE, FALSE}: TRUE for the WRs already provided for the Pacific/North American sector in 1948-2021
    
    start.date.synoptic="1948-01-01"; end.date.synoptic="2021-12-31" # from processed GPHA file
    
    start.date.WRs="1948-01-01"; end.date.WRs="2019-12-31" # proper leap year orders (starting with leap year of 1948, ending a year before (i.e., 2019) the leap year of 2020)
    
    dates.WRs.specific <- seq(as.Date(start.date.WRs),as.Date(end.date.WRs),by="day")
    
    num.years.sim.WRs <- number.years.long # e.g., 500, 1000, 2000, 3000, 5000 years, etc [note: current NHMM output (parametric) is for 1036 years]
    
    dir.to.sim.WRs.files <- "./Data/simulated.data.files/WRs.out" # dir.to.sim.WRs.files
    num.iter.WRs <- num.iter   #number of iterations to simulate sequence of WRs
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
    dynamic.scenario  <- 0 # {0, 1}: 0: no dynamic change; 1: yes dynamic change 
    
    if (dynamic.scenario==0){
      ##===> Attempt #0 (thermodynamic only; no change to freq of WRs) ===##
      # #specify target change (as a percent) for WR probabilities
      WR_prob_change <- c(0,0,0,0,0,0,0,0,0,0) # between 0 and 1
      # #how close (in % points) do the WR frequencies (probabilities) need to be to the target
      lp.threshold <- 0.00001
      # #how much change do we allow in a sub-period sampling probability before incurring a larger penalty in the optimization
      piecewise_limit <- .02
      
      #   --------- NOTE: some of these hyper-parameters may need tuning depending on the dynamic climate change selected
      #   --------- Attempt with caution!!!!
    }else if(dynamic.scenario==1){
      ##===> Attempt #1 (dynamic scenario #1) ===##
      # #specify target change (as a percent) for WR probabilities (if, increasing WR3 in future)
      WR_prob_change <- c(0,0,.3,0,0,0,0,0,0,0) # between 0 and 1
      # #how close (in % points) do the WR frequencies (probabilities) need to be to the target
      lp.threshold <- 0.007
      # #how much change do we allow in a sub-period sampling probability before incurring a larger penalty in the optimization
      piecewise_limit <- .02
      
      ##===> Other option explored in final report===##
      # specify target change (as a percent) for WR probabilities (if, continuing their current trends in future)
      #WR_prob_change <- c(-0.09969436,  0.27467048,  0.33848792,
      #                    -0.28431861, -0.23549986,  0.03889970,
      #                    -0.05628958, 0.38059153, -0.16636739, -0.17995965) # between 0 and 1
      # how close (in % points) do the WR frequencies (probabilities) need to be to the target
      #lp.threshold <- 0.008
      # how much change do we allow in a sub-period sampling probability before incurring a larger penalty in the optimization
      #piecewise_limit <- .02
      
    }
    
    
    # returning the entire values inserted here to 'run.stochastic.weather.generator'
    values = as.list(environment())
    return(values)
  }
}