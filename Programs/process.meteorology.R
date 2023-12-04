rm(list=ls())

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# This script provides the 'meteorology' input required for the WGEN run #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

#------------------------------------------------------------------
##/ Step 1: Importing raw meteorological time series \##

###----- raw data files MUST BE formatted as: yyyy mm dd P[mm] Tmax[C] Tmin[C]  *.csv file ----- ####

dir.to.all.raw.files <- "./Data/raw.data.files/"
list.locations <- list.files(dir.to.all.raw.files) # list of gridded location or station data
my.file <- read.table(paste0(dir.to.all.raw.files,list.locations[1]),sep= ",")
start_date <- paste(my.file[1,1],my.file[1,2],my.file[1,3],sep="-")
end_date <- paste(my.file[nrow(my.file),1],my.file[nrow(my.file),2],my.file[nrow(my.file),3],sep="-")

seq.of.dates <- seq(as.Date(start_date),as.Date(end_date),by = "day")
lst.grids.data <- list()

ascii.array.precip.temp <- array(NA,c(length(seq.of.dates),length(list.locations),3))
for (f in 1:length(list.locations)){

  my.file <- read.table(paste0(dir.to.all.raw.files,list.locations[f]),
                        sep= ",")
  ascii.array.precip.temp[,f,1] <- my.file[,4] # precip -- verify this
  ascii.array.precip.temp[,f,2] <- my.file[,5] # tmax -- verify this
  ascii.array.precip.temp[,f,3] <- my.file[,6] # tmin -- verify this
  prct.print <- round(f/length(list.locations)*100,2)
  print(paste('done -> location',f,'/',
              length(list.locations),', %',prct.print))
}


lst.import.datafile <- list('ascii.array.precip.temp'=ascii.array.precip.temp,
                            'seq.of.dates'=seq.of.dates,
                            'file.names'=list.locations)

saveRDS(lst.import.datafile,
        file='./Data/processed.data.files/processed.meteorology/lst.import.datafile.rds')



#------------------------------------------------------------------
##/ Step 2: Organizing input data for WGEN runs \##
rm(list=ls())

lst.import.datafile <- readRDS(file='./Data/processed.data.files/processed.meteorology/lst.import.datafile.rds')

ascii.array.precip.temp <- lst.import.datafile$ascii.array.precip.temp
dates.weather <- lst.import.datafile$seq.of.dates
list.file.names <- lst.import.datafile$file.names
rm(lst.import.datafile) # for memory

years.weather <- as.numeric(format(format(dates.weather,'%Y')))
months.weather <- as.numeric(format(format(dates.weather,'%m')))
wateryears.weather <- (years.weather+1)*(months.weather>=10) + (years.weather)*(months.weather<10)

n.sites <- dim(ascii.array.precip.temp)[2]
prcp.site <- ascii.array.precip.temp[,,1]
tmax.site <- ascii.array.precip.temp[,,2]
tmin.site <- ascii.array.precip.temp[,,3]

rm(ascii.array.precip.temp)

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

save(list = ls(envir = environment(), all.names = TRUE), 
     file = "./Data/processed.data.files/processed.meteorology/processed.meteorology.RData",
     envir = environment())

gc()

# end
