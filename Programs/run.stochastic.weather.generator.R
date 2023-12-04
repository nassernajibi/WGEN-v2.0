
rm(list=ls())

library(MASS) # Gamma fit
library(evmix) # GPD fit
library(eva) # GPD fit
library(depmixS4) # HMMs/NHMMs fit
library(markovchain) # HMMs/NHMMs fit
library(rebmix) # split/WRs
library(lpSolve) # lp optimization
library(mvtnorm) # MVN
library(lubridate) # dates
library(tictoc) # run time
library(moments) # computation
library(abind) # computation
#--------------------------------
library(zoo) # (plot)
library(fExtremes) # (plot)
library(scales) # (plot)
library(parallel) # (plot)
library(proxy) # (plot)
library(POT) # (plot) event-based computations
library(extRemes) # (plot)
library(ismev) # (plot)
library(readxl) # output

source("./Programs/config.simulations.R") # config file

lst <- config.simulations() # call in configuration inputs
for (i in 1:length(lst)) {assign(names(lst[i]), lst[[i]]) }; rm(lst)


#*************************************
#--- Weather Regimes Module ---#
if (use.provided.WRs){
  final.NHMM.output <- readRDS('./Data/simulated.data.files/WRs.out/final.NHMM.non_param.output.rds')
  weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
  num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
  dates.sim <- final.NHMM.output$dates.sim
  markov.chain.sim <- final.NHMM.output$WR.simulation
  dates.synoptics <- final.NHMM.output$dates.historical
} else{
  if (use.non_param.WRs){      #----------- 1+2 dynamic scenarios ----------#
    final.NHMM.output <- execute.WRs.non_param.NHMM()
    weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
    num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
    dates.sim <- final.NHMM.output$dates.sim
    markov.chain.sim <- final.NHMM.output$WR.simulation
    dates.synoptics <- final.NHMM.output$dates.historical
    
  } else {
    final.NHMM.output <- execute.WRs.param.NHMM()
    weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
    num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
    dates.sim <- final.NHMM.output$dates.sim
    markov.chain.sim <- as.list(data.frame(final.NHMM.output$WR.simulation[,1:num.iter]))
    dates.synoptics <- final.NHMM.output$dates.historical
  }
}
rm(final.NHMM.output) # for memory


#*************************************
#--- Weather Generation Module ---#
execute.simulations()
# done. #


# EXTRA #
### Below are auxiliary functions to do a list of tasks
#*************************************
# - create sample figures for selected scenario
# - generate individual output files in tab or text delimited formats

#this is the scenario (i.e., the row in ClimateChangeScenarios.csv) for which to make plots and write out the data as .csv files
selected_scenario = 1

#--- figures ---#
# arguments are labels for x and y-axes
start_time <- Sys.time()
create.figures.baselines.stacked(scenario = selected_scenario)
Sys.time() - start_time


#--- outputs ---#
# YYYY, MM, DD, P(mm), Tmax(C), Tmin(C) in .csv individual lat/lon file #
# for simulated data #
create.delimited.outputs(scenario = selected_scenario)

