rm(list=ls())
library(stringr) # for str_split

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# This script provides the 'meteohydro' input required for the WGEN run #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

#------------------------------------------------------------------
##/ Step 1: Importing raw meteorological time series \##
# creating a vector of dates for this dataset
start_date <- "1915-01-01"; end_date <- "2018-12-31"
seq.of.dates <- seq(as.Date(start_date),as.Date(end_date),by = "day")

dir.to.all.files <- "./Data/processed.data.files/SanFranciscoBay_HUC4_1805_data/"

list.grids <- list.files(dir.to.all.files)
lst.grids.data <- list()

ascii.array.precip.temp <- array(NA,c(length(seq.of.dates),length(list.grids),3))
ascii.lonlat <- array(NA,c(length(list.grids),2))
for (f in 1:length(list.grids)){
  dum <- unlist(str_split(list.grids[f],pattern='_'))
  ascii.lonlat[f,] <- cbind(dum[3],dum[2])
  
  my.file <- read.table(paste0(dir.to.all.files,list.grids[f]))
  ascii.array.precip.temp[,f,1] <- my.file[,4] # precip -- verify this
  ascii.array.precip.temp[,f,2] <- my.file[,5] # tmax -- verify this
  ascii.array.precip.temp[,f,3] <- my.file[,6] # tmin -- verify this
  prct.print <- round(f/length(list.grids)*100,2)
  print(paste('done -> gridded location',f,'out of',
              length(list.grids),', %',prct.print))
}



lst.import.datafile <- list('ascii.array.precip.temp'=ascii.array.precip.temp,
                            'ascii.lonlat'=ascii.lonlat,
                            'seq.of.dates'=seq.of.dates)

saveRDS(lst.import.datafile,
        file='./Data/processed.data.files/processed.meteohydro/lst.import.datafile.rds')



#------------------------------------------------------------------
##/ Step 2: Organizing for WGEN runs \##
rm(list=ls())

lst.import.datafile <- readRDS(file='./Data/processed.data.files/processed.meteohydro/lst.import.datafile.rds')

ascii.array.precip.temp <- lst.import.datafile$ascii.array.precip.temp
ascii.lonlat <- lst.import.datafile$ascii.lonlat
rm(lst.import.datafile)

start_date <- "1915-01-01"; end_date <- "2018-12-31"
dates.weather <- seq(as.Date(start_date),as.Date(end_date),by = "day")
years.weather <- as.numeric(format(format(dates.weather,'%Y')))
months.weather <- as.numeric(format(format(dates.weather,'%m')))
first.month <- 1; last.month <- 12; first.day <- 1; last.day <- 31
months <- seq(1,12)

wateryears.weather <- (years.weather+1)*(months.weather>=10) + (years.weather)*(months.weather<10)

n.sites <- dim(ascii.array.precip.temp)[2]
prcp.site <- ascii.array.precip.temp[,,1]
tmax.site <- ascii.array.precip.temp[,,2]
tmin.site <- ascii.array.precip.temp[,,3]
lon_lat <- ascii.lonlat

rm(ascii.array.precip.temp,ascii.lonlat)

prcp.basin <- apply(prcp.site,1,mean)
tmax.basin <- apply(tmax.site,1,mean)
tmin.basin <- apply(tmin.site,1,mean)

#create annual average for basin precipitation
prcp.basin.annual <- aggregate(prcp.basin,FUN=mean,by=list(years.weather),na.rm=T)
prcp.basin.annual[,2] <- scale(prcp.basin.annual[,2])[,1]
#create annual average for basin temperature mean.
tmin.basin.annual <- aggregate(tmin.basin,FUN=mean,by=list(years.weather),na.rm=T)
tmin.basin.annual[,2] <- scale(tmin.basin.annual[,2])[,1]
tmax.basin.annual <- aggregate(tmax.basin,FUN=mean,by=list(years.weather),na.rm=T)
tmax.basin.annual[,2] <- scale(tmax.basin.annual[,2])[,1]

file.name.to.save <- "processed.meteohydro"


numberdummy <- c('SFB')

save(list = ls(envir = environment(), all.names = TRUE), 
     file = paste0("./Data/processed.data.files/processed.meteohydro/",
                   file.name.to.save,".",numberdummy,".RData"),
     envir = environment())

gc()
# end
