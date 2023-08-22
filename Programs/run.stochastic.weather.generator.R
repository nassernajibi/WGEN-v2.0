
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
if (use.non_param.WRs){      #----------- 1+2 dynamic scenarios ----------#
  final.NHMM.output <- execute.WRs.non_param.NHMM(mainDir,
                                                  dir.to.sim.WRs.files,
                                                  num.iter.WRs,
                                                  seasons,
                                                  num_eofs.season,
                                                  num_WRs.season,
                                                  num.years.sim.WRs,
                                                  start.date.synoptic,
                                                  end.date.synoptic,
                                                  path.to.processed.SPI.PCs,
                                                  path.to.processed.GPHAs,
                                                  WR_prob_change,
                                                  lp.threshold,
                                                  piecewise_limit,
                                                  start.date.weather,
                                                  end.date.weather,
                                                  start.date.nonpar,
                                                  end.date.nonpar,
                                                  dynamic.scenario)
  weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
  num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
  dates.sim <- final.NHMM.output$dates.sim
  markov.chain.sim <- final.NHMM.output$WR.simulation
  
} else {
  final.NHMM.output <- execute.WRs.param.NHMM(mainDir,
                                              dir.to.sim.WRs.files,
                                              num.iter.WRs,
                                              seasons,
                                              num_eofs.season,
                                              num_WRs.season,
                                              num.years.sim.WRs,
                                              start.date.synoptic,
                                              end.date.synoptic,
                                              path.to.processed.SPI.PCs,
                                              path.to.processed.GPHAs,
                                              start.date.weather,
                                              end.date.weather,
                                              start.date.par,
                                              end.date.par)
  weather.state.assignments <- final.NHMM.output$WR.historical # this is the historical WRs 
  num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
  dates.sim <- final.NHMM.output$dates.sim
  markov.chain.sim <- as.list(data.frame(final.NHMM.output$WR.simulation[,1:num.iter]))
}
rm(final.NHMM.output) # for memory


#*************************************
#--- Weather Generation Module ---#
execute.simulations(mainDir,
                    dir.to.sim.files,
                    num.iter,
                    basin.cnt,
                    number.years.long,
                    change.list,
                    path.to.processed.data.meteohydro,
                    qq,
                    window.size,
                    pr.trace,
                    weather.state.assignments,
                    num.states,
                    dates.sim,
                    markov.chain.sim,
                    start.date.weather,
                    end.date.weather)
# done. #


# EXTRA #
### Below are auxiliary functions to do a list of tasks
#*************************************
# - create sample figures for baseline
# - generate individual output files in tab or text delimited formats
# - ...add yours!

#--- figures ---#
create.figures.baselines.stacked(mainDir,
                                 dir.to.sim.files,
                                 path.to.processed.data.meteohydro,
                                 num.iter,
                                 basin.cnt,
                                 change.list,
                                 use.non_param.WRs,
                                 num.states,
                                 qq,
                                 start.date.weather,
                                 end.date.weather,
                                 label1,
                                 label2)

#--- outputs ---#
# YYYY, MM, DD, P(mm), Tmax(C), Tmin(C) in .csv individual lat/lon file #
# for both simulated and observed data #
create.delimited.outputs(mainDir,
                         dir.to.sim.files,
                         dir.to.output.files,
                         path.to.processed.data.meteohydro,
                         num.iter,
                         basin.cnt,
                         change.list,
                         use.non_param.WRs,
                         num.states,
                         start.date.weather,
                         end.date.weather)

