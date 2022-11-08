
rm(list=ls())

library(depmixS4) # HMMs fit
library(rapportools)
library(markovchain)
library(rebmix)
library(moments)
library(MASS)
library(abind)
library(forecast)
library(biwavelet)
library(parallel)

#adjust main directory and directory for simulation files
mainDir <- "D:/Projects/Tuolumne_River_Basin/GitHub_WGENv2.0"
setwd(mainDir)

dir.to.sim.files <- "./Data/simulated.data.files/WRs.out"
dir.create(file.path(mainDir, dir.to.sim.files), showWarnings = FALSE)

# load in supporting functions
files.sources = list.files("./Programs/functions",full.names = TRUE)
my.functions <- sapply(files.sources, source)

num.iter <- 5   #number of iterations to simulate

#######################define seasons and covariates for NHMM models of WRs#############################
cold.months <- c(11,12,1,2,3,4) # Nov-Apr
warm.months <- c(5,6,7,8,9,10)  #May-Oct
seasons <- list(cold.months,warm.months)
n.seasons <- length(seasons)
num_eofs.season <- rep(10,n.seasons)  #number of PCs to use for geopotential heights per season
num_WRs.season <- c(7,3)    #number of WRs to fit per season

number.years.long <- 1000 # e.g., 1000 years; 2000 years, etc



#####Load in covariates for the cold season
#Covariates should be a matrix with the first column as dates, and the second column as 
# normalized pPC1 (scaled and centered)
SPI.dates.PCA.extracted <- readRDS(file='./Data/processed.data.files/processed.NHMM.data/paleo.norm.4.cold.PCs.dates_extracted.rds')
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
#// Load processed 500mb-GPH hgt region data and dates
hgt.synoptic.region <- readRDS(file='./Data/processed.data.files/processed.hgt/hgt.500.Pacific.NorthAmer.synoptic.region_19480101_20211231.rds')
start_date_synoptic="1948-01-01"; end_date_synoptic="2021-12-31"
dates.synoptic <- seq(as.Date(start_date_synoptic),as.Date(end_date_synoptic), by="days")
months.synoptic <- as.numeric(format(dates.synoptic,'%m'))

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
  
  ##intcpts <- intcpt.months.final[[s]]
  
  # my.formula <- list()
  # my.family <- list()
  # my.formula[[1]] <- as.formula("PC1~1")
  # my.family[[1]] <- gaussian()
  # for (eof in 2:n.eofs) {
  #   my.formula[[eof]] <- as.formula(paste("PC",eof,"~1",sep=""))
  #   my.family[[eof]] <- gaussian()
  # }
  # 
  # #fit initial HMM to get initial estimates of parameters
  # modHMMs <- depmix(my.formula,
  #                   nstates = n.states,
  #                   family=my.family,
  #                   ntimes =  nrow(synoptic.pcs),
  #                   data = data.frame(synoptic.pcs))
  # fit.mod.HMM <- fit(modHMMs)
  

  
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

saveRDS(fit.mod.NHMM,file = paste0(dir.to.sim.files,"/fit.mod.NHMM.rds"))
#fit.mod.NHMM <- readRDS(paste0(dir.to.sim.files,"/fit.mod.NHMM.rds"))
#####################################################################################



##################Provide arbitrary dates for any length##############################
##// following lines are for any length (incomplete)
# number.years.long <- 1036
# end_year <- as.numeric(format(as.Date(start_date_synoptic),"%Y"))+number.years.long
# vec_dates <- seq(as.Date(start_date_synoptic),as.Date(paste0(end_year,'-12-31')),"days")
# 
# target.long.dates <- list(vec_dates[as.numeric(format(vec_dates,"%m"))%in%seasons[[1]]],
#                           vec_dates[as.numeric(format(vec_dates,"%m"))%in%seasons[[2]]])
##//

#create dates for the simulated WRs
my.num.sim = ceiling(number.years.long/length(unique(format(dates.synoptic,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
long.dates.sim <- rep(dates.synoptic,times=my.num.sim)
target.long.dates <- list(long.dates.sim[as.numeric(format(long.dates.sim,"%m"))%in%seasons[[1]]],
                          long.dates.sim[as.numeric(format(long.dates.sim,"%m"))%in%seasons[[2]]])
s=1 # for cold season only in which we have paleo PCs
u <- covariates.final[[s]]
uRepped <- u[rep(seq_len(nrow(u)), my.num.sim),]
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
  #these two variables would be read in from separate files that were created elsewhere
  #/ for identical length as observation;
  #dates.sim.s[[s]] <- dates.final[[s]]
  #covariates.sim.s <- covariates.final[[s]]
  #//for longer length than observation;
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

final.NHMM.output <- list(dates.historical,WR.historical,dates.sim,WR.simulation)
names(final.NHMM.output) <- c('dates.historical','WR.historical','dates.sim','WR.simulation')
saveRDS(final.NHMM.output,file = paste0(dir.to.sim.files,"/final.NHMM.output.rds"))


gc()


#####################################################################################