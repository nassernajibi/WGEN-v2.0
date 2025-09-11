
# rm(list=ls())

# Load required libraries -------------------------------------------------

### Use pacman to automatically install missing libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)

p_load(MASS)        # Gamma fit
p_load(lpSolve)     # lp optimization
p_load(mvtnorm)     # MVN
p_load(moments)     # computation
p_load(abind)       # computation
p_load(lubridate)   # dates
p_load(depmixS4)    # HMMs/NHMMs fit
p_load(markovchain) # HMMs/NHMMs fit
p_load(rebmix)      # split/WRs
p_load(evmix)       # GPD fit
p_load(eva)         # GPD fit
p_load(POT)         # (plot) event-based computations
p_load(extRemes)    # (plot)
p_load(ismev)       # (plot)
p_load(fExtremes)   # (plot)
p_load(parallel)    # (plot)
p_load(zoo)         # (plot)
p_load(proxy)       # (plot)
p_load(scales)      # (plot)
p_load(readxl)      # output

start_time0 <- Sys.time()

# Load config file --------------------------------------------------------

source("./Programs/config.simulations.R") # config file

lst <- config.simulations() # call in configuration inputs
for (i in 1:length(lst)) {assign(names(lst[i]), lst[[i]]) }; rm(lst)



# Weather Regimes Module --------------------------------------------------

#use provided WRs
if (use.provided.WRs){
  final.NHMM.output <- readRDS('./Data/simulated.data.files/WRs.out/final.NHMM.non_param.output.rds')
  weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
  num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
  dates.sim <- final.NHMM.output$dates.sim
  markov.chain.sim <- final.NHMM.output$WR.simulation
  dates.synoptics <- final.NHMM.output$dates.historical
#simulate your own WRs
} else{
  final.NHMM.output <- execute.WRs.non_param.NHMM()
  weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
  num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
  dates.sim <- final.NHMM.output$dates.sim
  markov.chain.sim <- final.NHMM.output$WR.simulation
  dates.synoptics <- final.NHMM.output$dates.historical
}
rm(final.NHMM.output) # for memory



# Weather Generation Module -----------------------------------------------

start_time1 <- Sys.time()
execute.simulations(parallel = FALSE, number_of_cores = NULL)
# execute.simulations(parallel = TRUE, number_of_cores = 2)
Sys.time() - start_time1
invisible(gc())



# EXTRA -------------------------------------------------------------------

### Below are auxiliary functions to do a list of tasks
# - create sample figures for selected scenario
# - generate individual output files in tab or text delimited formats

#this is the scenario (i.e., the row in ClimateChangeScenarios.csv) for which to make plots and write out the data as .csv files
selected_scenario = 1

## Figures -----------------------------------------------------------------
# arguments are labels for x and y-axes
start_time2 <- Sys.time()
create.figures.baselines.stacked(scenario = selected_scenario)
Sys.time() - start_time2

## Output csv files --------------------------------------------------------

# YYYY, MM, DD, P(mm), Tmax(C), Tmin(C) in .csv individual lat/lon file #
# for simulated data #
start_time3 <- Sys.time()
create.delimited.outputs(scenario = selected_scenario)
Sys.time() - start_time3

end_time0 <- Sys.time()
time_taken0 <- end_time0 - start_time0

print(paste("Starting time: ", start_time0))
print(paste("Ending time: ", end_time0))
print(paste("Total runtime: "))
time_taken0
