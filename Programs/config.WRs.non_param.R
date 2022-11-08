
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

number.years.long <- 3050 # e.g., 500, 1000, 2000, 3000, 5000 years, etc [note: current NHMM output (parametric) is for 1036 years; current non-parametric is for 3050 years]


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

#create dates for the simulated WRs
my.num.sim = ceiling(number.years.long/length(unique(format(dates.synoptic,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
long.dates.sim <- rep(dates.synoptic,times=my.num.sim)
target.long.dates <- list(long.dates.sim[as.numeric(format(long.dates.sim,"%m"))%in%seasons[[1]]],
                          long.dates.sim[as.numeric(format(long.dates.sim,"%m"))%in%seasons[[2]]])

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


#here we glue together the different seasons
#for each season, put in the final simulations into the right places for the historical
#re-relabel the next season WRs for the final sequence of WRs
#dates.sim <- as.Date(sort(unlist(dates.sim.s)),origin="1970-01-01")
WR.historical.s <- list()
dates.sim <- long.dates.sim
dates.historical <- as.Date(sort(unlist(dates.final)),origin="1970-01-01")
WR.historical <- array(NA,length(dates.historical))
for (s in 1:n.seasons) {
  WR.historical.s[[s]] <- fit.mod.NHMM[[s]][[2]]
  if (s>1) {to.add <- sum(num_WRs.season[1:(s-1)])} else {to.add <- 0}     #relabel the WRs, based on WRs from previous seasons
  WR.historical[dates.historical%in%dates.final[[s]]] <- WR.historical.s[[s]]+to.add
}
#

#dates for the historical WRs
start_date <- "1950-01-01"; end_date <- "2013-12-31"
dates.weather <- seq(as.Date(start_date),as.Date(end_date),by="day")
identical.dates.idx <- dates.synoptic%in%dates.weather

weather.state.assignments <- WR.historical[identical.dates.idx]


#create dates for the simulated WRs
my.num.sim = ceiling(number.years.long/length(unique(format(dates.weather,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
long.dates.sim <- rep(dates.weather,times=my.num.sim)

years.weather <- format(dates.weather,'%Y')

years.vec <- unique(years.weather)
nyears <- length(years.vec)
n.segment <- 5 # length of segments for this case
my.n.segment <- base::split(1:nyears,ceiling(seq(1,nyears)/n.segment))
num.segment <- length(my.n.segment)
size.trace <- nyears*my.num.sim
total.seqments <- ceiling(size.trace/n.segment)
my.prob <- rep(1,num.segment) # equal prob
#my.prob <- c(rep(0.25,9),rep(0.5,6))
#my.prob <- c(rep(0.25,9),rep(0.75,6))
#my.prob <- c(rep(0.1,9),rep(0.9,6))
set.seed(n.segment)
my.replication <- sample(seq(1,num.segment),total.seqments,
                         replace = TRUE,
                         prob=my.prob)
dates.sim <- c()
markov.chain.sim <- c()
for (p in my.replication){
  
  chunk.years <- years.vec[my.n.segment[[p]]]
  idx.chunks <- which(years.weather%in%chunk.years)
  
  dates.sim <- append(dates.sim,
                      dates.weather[idx.chunks])
  markov.chain.sim <- append(markov.chain.sim,
                             weather.state.assignments[idx.chunks])
}
markov.chain.sim <- rep(list(markov.chain.sim),times=num.iter)


#final output to be saved and loaded into config.simulations:
#dates.sim
#markov.chain.sim
final.WRs.output <- list(dates.sim,markov.chain.sim)
names(final.WRs.output) <- c('dates.sim','markov.chain.sim')
saveRDS(final.WRs.output,file = paste0(dir.to.sim.files,"/final.non_param.WRs.output.rds"))


gc()


#####################################################################################