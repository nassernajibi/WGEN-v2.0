
rm(list=ls())

library(MASS)
library(fExtremes)
library(evmix)
library(zoo)
library(abind)
library(parallel)
library(mvtnorm)
library(tictoc)
library(lubridate)

#adjust main directory and directory for simulation files
mainDir <- "D:/Projects/Tuolumne_River_Basin/GitHub_WGENv2.0"
setwd(mainDir)

dir.to.sim.files <- "./Data/simulated.data.files/WGEN.out"
dir.create(file.path(mainDir, dir.to.sim.files), showWarnings = FALSE)

num.iter <- 1 # A single long trace (e.g., thousand years) is sufficient although we create like 5 ensembles in the simulated WRs
basin.cnt <- 'myPilot' # for a set of 12 randomly selected Livneh grids in the Tuolumne River Basin 
number.years.long <- 3050 # e.g., 500, 1000, 2000, 3000, 5000 years, etc [note: current NHMM output (parametric) is for 1036 years; current non-parametric is for 3050 years]

##############################Define perturbations#######################################
## Climate changes and jitter to apply:
change.list <- data.frame("tc"=  c(0), # e.g., 1, 2, ...
                          "jitter"=  c(TRUE),
                          "pccc"=c( 0), # e.g., 0.07, 0.14, ...
                          "pmuc"=c( 0)# e.g., -.125
)

# For simulating the WRs (i.e., 'config.WRs.non_param.R'; 'config.WRs.param.NHMM.R'), do you use non-parametric or parametric method
use.non_param.WRs <- TRUE #TRUE for non-parametric, FALSE for parametric simulated WRs
####### Choose A (FALSE) or B (TRUE) below for the simulated WRs #################
#-- A) load in NHMM with WRs (parametric)
tmp.list <- readRDS(paste0("./Data/simulated.data.files/WRs.out/final.NHMM.output.rds"))
weather.state.assignments <- tmp.list$WR.historical # this is for BOTH the parametric and non-parametric approach 
sim.weather.state.assignments <- tmp.list$WR.simulation[,1:num.iter] # this is only for the parametric approach of simulating WRs
num.states <- length(unique(as.vector(weather.state.assignments)))    #number of WRs in the model
rm(tmp.list) # for memory

#-- B) load in non-parametric simulation with WRs
np.list.sim.weather.state.assignments <- readRDS(paste0("./Data/simulated.data.files/WRs.out/final.non_param.WRs.output.rds")) #this is only for the non-parametric approach of simulating WRs
num.states <- length(unique(np.list.sim.weather.state.assignments$markov.chain.sim[[num.iter]]))    #number of WRs in the model

# load in supporting functions
files.sources = list.files("./Programs/functions",full.names = TRUE)
my.functions <- sapply(files.sources, source)
##################################

#dates for the historical WRs
start_date_synoptic="1948-01-01"; end_date_synoptic="2021-12-31"
dates.synoptic <- seq(as.Date(start_date_synoptic),as.Date(end_date_synoptic), by="days")
#create dates for the simulated WRs
my.num.sim = ceiling(number.years.long/length(unique(format(dates.synoptic,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
long.dates.sim <- rep(dates.synoptic,times=my.num.sim)

##################################


#########Precipitation characteristics#########

#location of obs weather data (RData format): weather data (e.g., precip and temp) as matrices (time x lat|lon: t-by-number of grids); dates vector for time; basin average precip (see the example meteohydro file)
processed.data.meteohydro <- paste0("./Data/processed.data.files/processed.meteohydro/processed.meteohydro.myPilot.RData")
load(processed.data.meteohydro) #load in weather data


qq <- .99              # percentile threshold to separate Gamma and GPD distributions
thshd.prcp <- apply(prcp.site,2,function(x) {quantile(x[x!=0],qq,na.rm=T)})


#Bootstrapping choices###
window.size <- rep(3,length(months))   #the size of the window (in days) from which runs can be bootstrapped around the current day of simulation, by month: Jan -- Dec

trace <- 0.25     #trace prcp threshold. 0.25 mm (for Livneh dataset); or 0.01 inches


#The spearman correlation between basin and site precipitation, used in the copula-based jitters
prcp.basin.site <- cbind(prcp.basin,prcp.site)
S <- cor(prcp.basin.site,method="spearman")
n.sites <- dim(prcp.site)[2] # Number of gridded points for precipitation


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

#decide whether or not to use historic or simulated WRs
if(use.non_param.WRs) {
  
  dates.sim <- np.list.sim.weather.state.assignments$dates.sim
  markov.chain.sim <- np.list.sim.weather.state.assignments$markov.chain.sim
  identical.dates.idx <- dates.synoptic%in%dates.weather
  weather.state.assignments <- weather.state.assignments[identical.dates.idx]
  
} else {
  dates.sim <- long.dates.sim
  markov.chain.sim <- as.list(data.frame(sim.weather.state.assignments))
}

#run the daily weather generate num.iter times using the num.iter Markov chains
mc.sim <- resampled.date.sim <- resampled.date.loc.sim <- prcp.site.sim <- tmin.site.sim <- tmax.site.sim <-  list()
start_time <- Sys.time()
for (k in 1:num.iter) {
  my.itertime <- Sys.time()
  my.sim <- wgen.simulator(weather.state.assignments=weather.state.assignments,mc.sim=markov.chain.sim[[k]],
                           prcp.basin=prcp.basin,dates.weather=dates.weather,
                           first.month=first.month,last.month=last.month,dates.sim=dates.sim,
                           months=months,window.size=window.size,min.thresh=trace)
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
  cur.tc <- change.list$tc[change]
  cur.jitter <- change.list$jitter[change]
  cur.pccc <- change.list$pccc[change]
  cur.pmuc <- change.list$pmuc[change]
  
  #precipitation scaling (temperature change dependent)
  perc.q <- (1 + cur.pccc)^cur.tc    #scaling in the upper tail for each month of non-zero prcp
  perc.mu <- (1 + cur.pmuc)          #scaling in the mean for each month of non-zero prcp
  
  #perturb the climate from the simulations above (the longest procedure in this function is saving the output files)
  set.seed(1)   #this ensures the copula-based jitterrs are always performed in the same way for each climate change
  perturbed.sim <-  perturb.climate(prcp.site.sim=prcp.site.sim,
                                    tmin.site.sim=tmin.site.sim,
                                    tmax.site.sim=tmax.site.sim,
                                    emission.fits.site=emission.fits.site,
                                    months=months,dates.sim=dates.sim,n.sites=n.sites,
                                    qq=qq,perc.mu=perc.mu,perc.q=perc.q,S=S,cur.jitter=cur.jitter,cur.tc=cur.tc,
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
  file.suffix <- paste0(".temp.",cur.tc,"_p.CC.scale.",cur.pccc,"_p.mu.scale.",cur.pmuc,"_hist.state.",use.non_param.WRs,"_jitter.",cur.jitter,"_s",num.states,"_with_",num.iter,".",basin.cnt)
  
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

# The End
#####################################################################################