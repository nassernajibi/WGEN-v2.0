rm(list=ls())
# for 1000 or 2000, ... 5000, years of single trace generated weather variables #
library(ggridges)
library(ggplot2)
library(ggstance)
library(viridisLite)
library(viridis)
library(gridExtra)
library(moments)
library(reshape2)
library(matrixStats)
library(RColorBrewer)
library(wesanderson)
library(parallel)
library(tictoc)
library(abind)
library(hydroGOF) # for 'index of agreement'
library(readxl)
library(ggthemes)
library(lubridate)
library(scales)
# -------------------------------------------
## Diagnostic individual basin plots ##
## comparing sim. vs obs. for a list of selected diagnostics
mainDir <- "D:/Projects/Tuolumne_River_Basin/GitHub_WGENv2.0"
setwd(mainDir)

basin.cnt <- 'myPilot' # for a set of 12 randomly selected Livneh grids in Tuolumne River Basin 

use.non_param.WRs <- TRUE #TRUE for non-parametric, FALSE for parametric simulated WRs
num.states <- c(10) # K WRs
num.iter <- 1 # Number of ensemble members generated using WR-SWG


dir.to.sim.files <- "./Data/simulated.data.files/WGEN.out"

##// weather and synoptic dates ---##
#location of obs weather data
processed.data.meteohydro <- paste0("./Data/processed.data.files/processed.meteohydro/processed.meteohydro.myPilot.RData")
load(processed.data.meteohydro) #load in weather data

n.sites <- dim(prcp.site)[2] # Number of gridded points for precipitation

# for parametric
number.years.long <- 1000 # e.g., 1000 years; 2000 years, etc
start_date_synoptic="1948-01-01"; end_date_synoptic="2021-12-31"
dates.synoptic <- seq(as.Date(start_date_synoptic),as.Date(end_date_synoptic), by="days")
my.num.sim = ceiling(number.years.long/length(unique(format(dates.synoptic,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical

# for non-parametric
number.years.long <- 3050 # e.g., 1000 years; 2000 years, etc
start_date <- "1950-01-01"; end_date <- "2013-12-31"
dates.weather <- seq(as.Date(start_date),as.Date(end_date),by="day")
my.num.sim = ceiling(number.years.long/length(unique(format(dates.weather,'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
np.list.sim.weather.state.assignments <- readRDS(paste0("./Data/simulated.data.files/WRs.out/final.non_param.WRs.output.rds"))
dates.sim <- np.list.sim.weather.state.assignments$dates.sim

# Sim. file
change = 1
change.list <- data.frame("tc"=  c(0), # e.g., 1, 2, ...
                          "jitter"=  c(TRUE),
                          "pccc"=c( 0), # e.g., 0.07, 0.14, ...
                          "pmuc"=c( 0)# e.g., -.125
)
cur.tc <- change.list$tc[change]
cur.jitter <- change.list$jitter[change]
cur.pccc <- change.list$pccc[change]
cur.pmuc <- change.list$pmuc[change]

simulated.file.run.model.saved <- paste0(".temp.",cur.tc,"_p.CC.scale.",cur.pccc,"_p.mu.scale.",cur.pmuc,"_hist.state.",use.non_param.WRs,"_jitter.",cur.jitter,"_s",num.states,"_with_",num.iter,".",basin.cnt)
prcp.site.sim_sfx <- "prcp.site.sim"
load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
prcp.site.data.sim <- prcp.site.sim

prcp.site.sim_sfx <- "tmin.site.sim"
load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
tmin.site.data.sim <- tmin.site.sim

prcp.site.sim_sfx <- "tmax.site.sim"
load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
tmax.site.data.sim <- tmax.site.sim

prcp.site.sim_sfx <- "dates.sim"
load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
dates.sim.hist <- dates.sim

months.sim <- months.sim.hist <- format(dates.sim.hist,"%m")
years.long.dates.sim <- years.sim <- years.sim.hist <- as.numeric(format(dates.sim.hist,"%Y"))

# water years dates #
years.weather <- as.numeric(format(dates.weather,"%Y"))
months.weather <- as.numeric(format(dates.weather,"%m"))
years.synoptic <- format(dates.synoptic,"%Y")
num.yr <- length(unique(years.weather))
years.sim <- as.numeric(format(dates.sim,"%Y"))

#total.num.yrs <- my.num.sim*length(unique(years.sim))
total.num.yrs <- round(sum(as.numeric(table(years.sim)))/365.25)

dates.sim.hist.DOY <- as.numeric(format(dates.sim.hist,'%j'))
dates.sim.hist.day.in.month <- as.numeric(format(dates.sim.hist,'%d'))
dates.sim.hist.year <- unique(as.numeric(format(dates.sim.hist,'%Y')))
dates.sim.hist.month <- as.numeric(format(dates.sim.hist,'%m'))
DOY.in.long.trace <- dates.sim.hist.DOY
MONTH.in.long.trace <- dates.sim.hist.month
DAY.in.long.trace <- dates.sim.hist.day.in.month
#days.in.yr <- as.numeric(table(years.sim)/my.num.sim)
all.days.in.yr <- c(abs(diff(DOY.in.long.trace)-1)[abs(diff(DOY.in.long.trace)-1)>0],tail(dates.sim.hist.DOY,1))
YEAR.in.long.trace <- NULL
for(z in 1:total.num.yrs){
  YEAR.in.long.trace <- append(YEAR.in.long.trace,
                               rep(z,times=all.days.in.yr[z]))
}

long.dates.sim.hist <- cbind(YEAR.in.long.trace,MONTH.in.long.trace,DAY.in.long.trace)
long.dates.sim.vec <- long.dates.sim.hist[,1]*10000+long.dates.sim.hist[,2]*100+long.dates.sim.hist[,3]
long.years.vec <- YEAR.in.long.trace
long.months.vec <- MONTH.in.long.trace

wateryears.sim.new <- (long.years.vec+1)*(long.months.vec>=10) + (long.years.vec)*(long.months.vec<10)
wateryears.obs.new <- (years.weather+1)*(months.weather>=10) + (years.weather)*(months.weather<10)

site.data <- prcp.site
site.data.sim <- prcp.site.data.sim

wy.data.obs <- apply(site.data,2,function(x){
  aggregate(x,FUN=sum,by=list(wateryears.obs.new),na.rm=T)[,2]})
wy.data.obs <- apply(wy.data.obs,2,function(x){
  return(x[2:(length(x)-1)])})#drop first and last to only keep full water years

wy.data.sim <- apply(site.data.sim[[1]],2,function(x){
  aggregate(x,FUN=sum,by=list(wateryears.sim.new),na.rm=T)[,2]})
wy.data.sim <- apply(wy.data.sim,2,function(x){
  return(x[2:(length(x)-1)])})#drop first and last to only keep full water years


###///////////////////###
##-- Precipitation --##
###///////////////////###

# START #
#================================================================================
## Plot: severity of long-term droughts  #
#================================================================================
prcp.min.extreme.diagnostics <- function(prcp,my.range=NULL,mc.cores=1) {
  d <- length(prcp)
  n.sites <- dim(prcp[[1]])[2]
  #ma.options.rolling <- c(1095,1825,3650) # n-day moving averages
  ma.options.rolling <- seq(1,20)*365 # n-day moving averages
  #stat.names2 <- c('3.year.min','5.year.min','10.year.min')
  stat.names2 <- paste0(seq(1,20),'.year.min')
  #calculate at-site statistics: min
  n.stats <- length(stat.names2)
  if(is.null(my.range)) {my.range <- 1:dim(prcp[[1]])[1]}
  stat.list2 <- mclapply(X=1:d,mc.preschedule=TRUE,mc.cores=mc.cores,FUN=function(i){  
    my.stats2 <- array(NA,c(n.sites,n.stats))
    cur.prcp <- prcp[[i]][my.range,]
    if (sum(is.infinite(cur.prcp))>0){
      idx.dummy <- which(is.infinite(cur.prcp),arr.ind = T)
      cur.prcp[idx.dummy[,1],idx.dummy[,2]] <- NA} # checkpoints to not have "Inf" before going ahead
    for (s in 1:n.stats) {
      my.stats2[,s] <- apply(cur.prcp,2,function(x){min(rollapply(data=x, FUN = mean, width = ma.options.rolling[s]))})
    }
    colnames(my.stats2) <- stat.names2
    return(my.stats2)
  })
  my.stats2 <- unlist(stat.list2)
  dim(my.stats2) <- c(n.sites,n.stats,d)    
  dimnames(my.stats2)[[2]] <- stat.names2
  return(my.stats2)
}

years.sim <- years.sim.hist; months.sim <- months.sim.hist
n.dig <- 5

#individual site diagnostics
stat.names2 <- paste0(seq(1,20),'.year')
stat.names22 <- paste0(seq(1,20))
obs.diagnostics.min <- prcp.min.extreme.diagnostics(list(site.data))
obs.diagnostics.min.extreme <- as.matrix(obs.diagnostics.min[,,1])
obs.stats.min <- round(obs.diagnostics.min.extreme,digits = n.dig)

sim.diagnostics.min <- prcp.min.extreme.diagnostics(site.data.sim)
sim.min.median <- round(apply(abind(sim.diagnostics.min,along = 3),FUN=median,c(1,2),na.rm=T),5)
sim.min.q975 <- round(apply(abind(sim.diagnostics.min,along = 3),FUN=quantile,c(1,2),0.975, na.rm=T), digits = n.dig)
sim.min.q025 <- round(apply(abind(sim.diagnostics.min,along = 3),FUN=quantile,c(1,2),0.025,na.rm=T), digits = n.dig)

#n-year drought across all sites (cumulative precip) for more n-year inputs, i.e., 1 to 20-yr droughts
ww.mom <- 12; hh.mom <- 12
fig.stat.name1 <- paste0("Fig.n.year.drought.",basin.cnt,"_s",
                         num.states,"_with_",num.iter,"_ens_together.all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)

n.yr.obs.diagnostics.min.extreme.all <- colSums(obs.diagnostics.min.extreme)
n.yr.sim.diagnostics.min.all <- apply(sim.diagnostics.min,FUN=sum,c(2))

min.obs.sim <- min(c(n.yr.obs.diagnostics.min.extreme.all,n.yr.sim.diagnostics.min.all))
max.obs.sim <- max(c(n.yr.obs.diagnostics.min.extreme.all,n.yr.sim.diagnostics.min.all))
par(mfrow=c(1,1),mar=c(5.1, 6.1, 4.1, 2.1))
plot(n.yr.sim.diagnostics.min.all,
     n.yr.obs.diagnostics.min.extreme.all,
     main='Worst Droughts',font.main = 1,
     cex.axis=2,cex.lab=1.5,cex.main=2.25,frame=F,
     xlab=paste0('Sim [WGEN: ',total.num.yrs,'-yr]'),cex=2,xlim=c(min.obs.sim,
                                         max.obs.sim),
     ylim=c(min.obs.sim,
            max.obs.sim),
     ylab='Obs',col=seq(1,20),pch=seq(1,20))
abline(h=n.yr.obs.diagnostics.min.extreme.all,
       lty=2,col=alpha(seq(1,20),0.25))
abline(v=n.yr.sim.diagnostics.min.all,
       lty=2,col=alpha(seq(1,20),0.25))
abline(a=0,b=1,lty=2,
       col='red',lwd=1)
legend("bottomright", legend=stat.names22,pch=seq(1,20),cex=2,
       col=seq(1,20),title ='n-year drought',
       horiz = FALSE,inset=0.01,ncol=2)

dev.off()
#


#================================================================================
## Plot: frequency of and/or annual rolling mean drought  #
#================================================================================
ww.mom <- 20; hh.mom <- 12
fig.stat.name1 <- paste0("Fig.rolling.mean_precip_Obs.vs.Sim.",basin.cnt,"_s",
                        num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(4,5),mar=c(5.1, 6.1, 4.1, 2.1))
drought.rol.mean <- seq(1,20)
labs.metrics <- paste0(drought.rol.mean,'-yr')
annual.prcp <- apply(wy.data.obs,1,mean)
annual.prcp.sim <- apply(wy.data.sim,1,mean)
for (my.opt in drought.rol.mean){
  obs.HRU <- rollmean(annual.prcp,drought.rol.mean[my.opt])
  wgen.trace <- rollmean(annual.prcp.sim,drought.rol.mean[my.opt])
  x_min_max <- c(min(obs.HRU,wgen.trace),max(obs.HRU,wgen.trace))
  hgA <- hist(wgen.trace,plot = FALSE,freq = T,
              main=paste0(labs.metrics[my.opt]),font.main = 1,
              cex.axis=2,cex.lab=1.5,cex.main=2.25,
              xlab='Precipitation [WY]',xlim=x_min_max,
              ylab='frequency [#WYs]',col='gray60',lty="blank",breaks = 25)
  hgB <- hist(obs.HRU,plot = FALSE,freq = T,
              main=paste0(labs.metrics[my.opt]),breaks = 25)
  plot(hgA, col = alpha('gray',0.3),lty="blank",
       main=paste0(labs.metrics[my.opt]),font.main = 1,
       cex.axis=2,cex.lab=1.5,cex.main=2.25,
       xlab='Precipitation [WY]',xlim=x_min_max,
       ylab='#WYs') # Plot 1st histogram using a transparent color
  plot(hgB, col = alpha('red',0.35), add = TRUE,lty="blank") # Add 2nd histogram using different color
  abline(v=median(wgen.trace),
         col='black',lwd=2,lty=1)
  abline(v=median(obs.HRU),
         col='red',lwd=2,lty=1)
  abline(v=min(wgen.trace),
         col='black',lwd=1.5,lty=2)
  abline(v=min(obs.HRU),
         col='red',lwd=1.5,lty=2)
  abline(v=max(wgen.trace),
         col='black',lwd=1.5,lty=2)
  abline(v=max(obs.HRU),
         col='red',lwd=1.5,lty=2)
  if (my.opt==1){
    legend("topright",
           legend=c('Sim [min,50th,max]','Obs [min,50th,max]'),
           col=c('gray','red'),
           cex=1.5,
           title ='',fill=c('gray','red'),
           bty='n',border = c(NA,NA),
           horiz = FALSE,inset=0.01)
  }
}
dev.off()


#================================================================================
## Plot: GEV, NEP, and flood plots  #
#================================================================================
my.gev.rl <- function (a, mat, dat, ymax, site.num) 
{
  require(ismev)
  eps <- 1e-06
  a1 <- a
  a2 <- a
  a3 <- a
  a1[1] <- a[1] + eps
  a2[2] <- a[2] + eps
  a3[3] <- a[3] + eps
  f <- c(seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5, 
         0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999, 0.99999)
  q <- gevq(a, 1 - f)
  d <- t(gev.rl.gradient(a = a, p = 1 - f))
  v <- apply(d, 1, q.form, m = mat)
  
  #ymax <- max(dat, q+ 1.96 * sqrt(v))
  ymax <- max(dat, ymax, q+ 1.96 * sqrt(v))
  ymin <- min(dat, q- 1.96 * sqrt(v))
  
  plot(-1/log(f), q, log = "x", type = "n", xlim = c(0.1, 7000), panel.first=grid(equilogs=TRUE,lty=1,col="gray95",lwd=0.5),
       ylim = c(ymin, ymax), xlab = "Return Period [year]",font.main = 1, 
       ylab = "Precipitation [Return Level]",xaxt = 'n',
       cex=2,main=paste0('Site ',site.num),
       cex.axis=1.5,cex.lab=1.25,cex.main=2, frame=FALSE)
  # newx <- -1/log(f)
  # newy <- q
  # polygon(c(rev(newx), newx), c(rev(newy),newy), density=20, col = rgb(0,1,0,0.3), border = NA)
  
  lines(-1/log(f), q, col='red',lty=1, lwd=1.5)
  lines(-1/log(f), q + 1.96 * sqrt(v), col = alpha("red",.3), lty=2, lwd=0.75)
  lines(-1/log(f), q - 1.96 * sqrt(v), col = alpha("red",.3), lty=2, lwd=0.75)
  points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat),col='red', pch=4)
}

fct.prcp.sim.annmax <- function(prcp,my.years,my.range=NULL,mc.cores=1){
  uniq.years <- unique(my.years)
  d <- length(prcp)
  n.sites <- dim(prcp[[1]])[2]
  if(is.null(my.range)) {my.range <- 1:dim(prcp[[1]])[1]}
  stat.rp.rl.list <- mclapply(X=1:d,mc.preschedule=TRUE,mc.cores=mc.cores,FUN=function(i){  
    cur.prcp <- prcp[[i]][my.range,]
    if (sum(is.infinite(cur.prcp))>0){
      idx.dummy <- which(is.infinite(cur.prcp),arr.ind = T)
      cur.prcp[idx.dummy[,1],idx.dummy[,2]] <- NA} # checkpoints to not have "Inf" before going ahead
    
    annual.my.stats <- apply(cur.prcp,2,function(x){
      amax <- aggregate(x,FUN=max,by=list(my.years))[,2]
      return(amax)
    })
  })
  my.stats <- unlist(stat.rp.rl.list)
  dim(my.stats) <- c(length(uniq.years),n.sites,d)
  return(my.stats)
}

# obs annual max
psite.obs.amax <- sapply(1:n.sites,function(y){
  return(aggregate(prcp.site[,y],FUN=max,by=list(wateryears.obs.new))[,2])})
# sim annual max
prcp.sim.annmax0 <- fct.prcp.sim.annmax(prcp=prcp.site.sim,my.years=wateryears.sim.new)
prcp.sim.annmax <- apply(prcp.sim.annmax0,2,FUN=as.vector)
ymax <- max(c(psite.obs.amax,prcp.sim.annmax))

file.name.fig <- paste0("Fig_obs.GEV_sim.NEP.",basin.cnt,"_s",
                        num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
fname.fig <- paste0("./Figures/",file.name.fig)
ww.mom <- 12; hh.mom <- 10
png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(3,4))
for (site.num in 1:n.sites){
  require(ismev)
  my.gev.fit1 <- gev.fit(psite.obs.amax[,site.num],show = F)
  my.gev.rl(my.gev.fit1$mle, my.gev.fit1$cov, my.gev.fit1$data,ymax,site.num)
  dat <- prcp.sim.annmax[,site.num]
  points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat),
         col="black",cex=.3,pch=19)
  myTicks = axTicks(1)
  axis(1, at = myTicks, labels = formatC(myTicks, format = 'd'),cex.axis=1.15)
  if(site.num==1){
    legend("topright", legend=c("obs","sim (WGEN)"),
           col=c("red", "black"), lty=c(1,1),cex=1.5,pch=c(4,16),
           bty = "n",horiz = FALSE)
  }
}
dev.off()
#

#================================================================================
## Plot: monthly distributions #
#================================================================================
years.months.weather <- format(dates.weather,'%Y-%m')
unique.years.months.weather <- format(ym(unique(years.months.weather)),'%m')

dummy.long.years.vec <- long.years.vec+2000
years.months.sim <- paste0(dummy.long.years.vec,'-',long.months.vec)
unique.years.months.sim <- format(ym(unique(years.months.sim)),'%m')

psite.obs.metric <- sapply(1:n.sites,function(y){
  return(aggregate(prcp.site[,y],FUN=sum,by=list(years.months.weather))[,2])})
psite.obs.metric2 <- sapply(1:n.sites,function(y){
  return(aggregate(psite.obs.metric[,y],FUN=median,by=list(unique.years.months.weather))[,2])})
psite.obs.metric.q1th <- sapply(1:n.sites,function(y){
  return(aggregate(psite.obs.metric[,y],FUN=quantile,0.05,by=list(unique.years.months.weather))[,2])})
psite.obs.metric.q2th <- sapply(1:n.sites,function(y){
  return(aggregate(psite.obs.metric[,y],FUN=quantile,0.95,by=list(unique.years.months.weather))[,2])})
  
psite.sim.metric <- lapply(1:num.iter,function(y){
  return(aggregate(prcp.site.sim[[y]],FUN=sum,by=list(years.months.sim))[,-1])})
psite.sim.metric2 <- lapply(1:num.iter,function(y){
  return(aggregate(psite.sim.metric[[y]],FUN=median,by=list(unique.years.months.sim))[,-1])})
psite.sim.metric.q1th <- lapply(1:num.iter,function(y){
  return(aggregate(psite.sim.metric[[y]],FUN=quantile,0.05,by=list(unique.years.months.sim))[,-1])})
psite.sim.metric.q2th <- lapply(1:num.iter,function(y){
  return(aggregate(psite.sim.metric[[y]],FUN=quantile,0.95,by=list(unique.years.months.sim))[,-1])})

##reordering.dummy.aggre <- c(1,10,11,12,2,3,4,5,6,7,8,9) #stupidty of aggregate for sim
reordering.dummy.aggre <- c(1,5,6,7,8,9,10,11,12,2,3,4)
psite.sim.metric2 <- lapply(1:num.iter,function(y){
  return(psite.sim.metric2[[y]][reordering.dummy.aggre,])})
psite.sim.metric.q1th <- lapply(1:num.iter,function(y){
  return(psite.sim.metric.q1th[[y]][reordering.dummy.aggre,])})
psite.sim.metric.q2th <- lapply(1:num.iter,function(y){
  return(psite.sim.metric.q2th[[y]][reordering.dummy.aggre,])})


file.name.fig <- paste0("Fig_monthly.dist.obs_long.sim.",basin.cnt,"_s",
                        num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
fname.fig <- paste0("./Figures/",file.name.fig)
ww.mom <- 12; hh.mom <- 10
png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(3,4))
for (k in 1:num.iter){
  for (site.num in 1:n.sites){
    
    plot(months,psite.sim.metric2[[k]][,site.num],
         pch=3,col='black',cex=1.5,ylab='Precipitation [total]',frame=F,
         font.main = 1,sub='',
         cex.lab=1.5,cex.axis=1.5,cex.main=2,
         main=paste0('Site ',site.num),
         ylim=c(min(psite.sim.metric.q1th[[k]][,site.num]),max(psite.sim.metric.q2th[[k]][,site.num])))
    arrows(months,psite.obs.metric.q1th[,site.num]
           ,months,psite.obs.metric.q2th[,site.num],
           length = 0,lwd=1,col=alpha('red',0.65))
    
    arrows(months,psite.sim.metric.q1th[[k]][,site.num]
                      ,months,psite.sim.metric.q2th[[k]][,site.num],
          length = 0,lwd=6,col=alpha('gray',0.5))
    points(psite.obs.metric2[,site.num],col="red",cex=1,pch=19)
    if(site.num==1){
      legend("top", legend=c("Sim (WGEN: 5th,median,95th)","Obs (5th,median,95th)"),
             col=c("black", "red"), lty=c(1,0),cex=1.15,
             pch=c(3,19),
             horiz = FALSE)
    }
  }
}
dev.off()


#================================================================================
## Plot: water-year cumulative precipitations #
#================================================================================
years.months.weather <- format(dates.weather,'%Y-%m')
unique.years.months.weather <- format(ym(unique(years.months.weather)),'%m')
doy.obs <- format(dates.weather,'%j')

years.months.sim <- paste0(years.long.dates.sim,'-',months.sim.hist)
unique.years.months.sim <- format(ym(unique(years.months.sim)),'%m')
doy.sim <- format(dates.sim,'%j')

sample.WY.dates <- c(dates.weather[275:366],dates.weather[1:274]) # with leap year: 1948
sample.WY.months.days <- format(sample.WY.dates,'%m-%d')
sample.WY.months <- format(sample.WY.dates,'%m')
sample.WY.Months <- format(sample.WY.dates,'%b')
doy.WY.sample.idx <- as.numeric(format(c(dates.weather[275:366],dates.weather[1:274]),'%j'))
seq.days <- sort(doy.WY.sample.idx)

start_wy <- which(format(dates.weather,'%m-%d')=='10-01')
end_wy <- which(format(dates.weather,'%m-%d')=='09-30')
num.WY.available <- length(start_wy)-1
cumsum.WY.obs.sites <- array(NA,c(num.WY.available,n.sites,366))
for(y in 1:num.WY.available){
  extracted.wy.dates <- dates.weather[start_wy[y]:end_wy[y+1]]
  doy.WY.single.wy.idx <- as.numeric(format(extracted.wy.dates,'%j'))
  
  for (k in 1:n.sites){
    
    single.wy <- prcp.site[start_wy[y]:end_wy[y+1],k]
    cumsum.WY.obs.sites[y,k,(1:length(single.wy))] <- cumsum(single.wy)
  }
}
psite.obs.mean.metric2 <- apply(cumsum.WY.obs.sites,c(2,3),FUN=mean,na.rm=T)
psite.obs.metric2 <- apply(cumsum.WY.obs.sites,c(2,3),FUN=median,na.rm=T)
psite.obs.metric2.10th <- apply(cumsum.WY.obs.sites,c(2,3),FUN=quantile,0.1,na.rm=T)
psite.obs.metric2.90th <- apply(cumsum.WY.obs.sites,c(2,3),FUN=quantile,0.9,na.rm=T)

q=1 # single iteration number
start_wy <- which(format(dates.sim,'%m-%d')=='10-01')
end_wy <- which(format(dates.sim,'%m-%d')=='09-30')
num.WY.available <- length(start_wy)-1
cumsum.WY.sim.sites <- array(NA,c(num.WY.available,n.sites,366))
for(y in 1:num.WY.available){
  extracted.wy.dates <- dates.sim[start_wy[y]:end_wy[y+1]]
  doy.WY.single.wy.idx <- as.numeric(format(extracted.wy.dates,'%j'))
  
  for (k in 1:n.sites){
    
    single.wy <- prcp.site.sim[[q]][start_wy[y]:end_wy[y+1],k]
    cumsum.WY.sim.sites[y,k,(1:length(single.wy))] <- cumsum(single.wy)
  }
}

psite.sim.mean.metric2 <- apply(cumsum.WY.sim.sites,c(2,3),FUN=mean,na.rm=T)
psite.sim.metric2 <- apply(cumsum.WY.sim.sites,c(2,3),FUN=median,na.rm=T)
psite.sim.metric2.10th <- apply(cumsum.WY.sim.sites,c(2,3),FUN=quantile,0.1,na.rm=T)
psite.sim.metric2.90th <- apply(cumsum.WY.sim.sites,c(2,3),FUN=quantile,0.9,na.rm=T)

start_months <- ceiling_date(sample.WY.dates %m-% months(1), 'month')
idx.fst.d.month <- which(sample.WY.dates %in% start_months)
my.seq.days <- 1:365

# median
file.name.fig <- paste0("Fig_monthly.cumsum.50th.10th_90th.obs_long.sim.",basin.cnt,"_s",
                        num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
fname.fig <- paste0("./Figures/",file.name.fig)
ww.mom <- 18; hh.mom <- 10
png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(3,4))
for (site.num in 1:n.sites){
  
  plot(my.seq.days,psite.sim.metric2[site.num,my.seq.days],
       type='l',col='black',cex=1.5,ylab='Precipitation [cumulative]',frame=F,
       font.main = 1,sub='',xaxt = "n",xlab='WY',lwd=1.5,ylim=c(min(c(psite.sim.metric2.10th[site.num,],
                                                                      psite.obs.metric2.10th[site.num,])),
                                                                max(c(psite.sim.metric2.90th[site.num,],
                                                                      psite.obs.metric2.90th[site.num,]))),
       cex.lab=1.5,cex.axis=1.5,cex.main=2,
       main=paste0('Site ',site.num))
  
  axis(1,at = my.seq.days[idx.fst.d.month],col='gray50',cex.lab=1.5,cex.axis=1,
       labels = sample.WY.Months[idx.fst.d.month])
  
  abline(v=my.seq.days[idx.fst.d.month],col='gray90',lty=2)
  
  lines(psite.obs.metric2[site.num,my.seq.days],
        col="red",cex=1,lwd=1.5)
  
  lines(psite.sim.metric2.10th[site.num,my.seq.days],
        col="gray45",cex=1,lwd=1,lty=2)
  lines(psite.sim.metric2.90th[site.num,my.seq.days],
        col="gray45",cex=1,lwd=1,lty=2)
  
  lines(psite.obs.metric2.10th[site.num,my.seq.days],
        col=rgb(1,0,0,0.5),cex=1,lwd=1,lty=2)
  lines(psite.obs.metric2.90th[site.num,my.seq.days],
        col=rgb(1,0,0,0.5),cex=1,lwd=1,lty=2)
  

  if(site.num==1){
    legend("bottomright", legend=c("Sim (WGEN: 10,median,90th)",
                                   "Obs (10,median,90th)"),
           col=c("black", "red"), lty=c(1,1),lwd=c(1.5,1.5),cex=1.5,
           horiz = FALSE)
  }
}
dev.off()


# mean
file.name.fig <- paste0("Fig_monthly.cumsum.mean.10th_90th.obs_long.sim.",basin.cnt,"_s",
                        num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
fname.fig <- paste0("./Figures/",file.name.fig)
ww.mom <- 18; hh.mom <- 10
png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(3,4))
for (site.num in 1:n.sites){
  
  plot(my.seq.days,psite.sim.mean.metric2[site.num,my.seq.days],
       type='l',col='black',cex=1.5,ylab='Precipitation [cumulative]',frame=F,
       font.main = 1,sub='',xaxt = "n",xlab='WY',lwd=1.5,ylim=c(min(c(psite.sim.metric2.10th[site.num,],
                                                                      psite.obs.metric2.10th[site.num,])),
                                                                max(c(psite.sim.metric2.90th[site.num,],
                                                                      psite.obs.metric2.90th[site.num,]))),
       cex.lab=1.5,cex.axis=1.5,cex.main=2,
       main=paste0('Site ',site.num))
  
  axis(1,at = my.seq.days[idx.fst.d.month],col='gray50',cex.lab=1.5,cex.axis=1,
       labels = sample.WY.Months[idx.fst.d.month])
  
  abline(v=my.seq.days[idx.fst.d.month],col='gray90',lty=2)
  
  lines(psite.obs.mean.metric2[site.num,my.seq.days],
        col="red",cex=1,lwd=1.5,lty=1)
  
  lines(psite.sim.metric2.10th[site.num,my.seq.days],
        col="gray45",cex=1,lwd=1,lty=2)
  lines(psite.sim.metric2.90th[site.num,my.seq.days],
        col="gray45",cex=1,lwd=1,lty=2)
  
  lines(psite.obs.metric2.10th[site.num,my.seq.days],
        col=rgb(1,0,0,0.5),cex=1,lwd=1,lty=2)
  lines(psite.obs.metric2.90th[site.num,my.seq.days],
        col=rgb(1,0,0,0.5),cex=1,lwd=1,lty=2)
  
  
  if(site.num==1){
    legend("bottomright", legend=c("Sim (WGEN: 10,mean,90th)",
                                   "Obs (10,mean,90th)"),
           col=c("black", "red"), lty=c(1,1),lwd=c(1.5,1.5),cex=1.5,
           horiz = FALSE)
  }
}
dev.off()



#================================================================================
## Plot: short/long-term extreme precipitations #
#================================================================================
# less than a year
#n-days/years precip max
prcp.max.extreme.diagnostics <- function(prcp,my.years,my.range=NULL,mc.cores=1) {
  d <- length(prcp)
  n.sites <- dim(prcp[[1]])[2]
  ma.options.rolling <- c(7,10,30,90,180,seq(1,10)*365) # n-day moving averages
  stat.names2 <- c(paste0(c(7,10),'.day.max'),paste0(c(1,3,6),'.month.max'),paste0(seq(1,10),'.year.max'))
  #calculate at-site statistics: min
  n.stats <- length(stat.names2)
  if(is.null(my.range)) {my.range <- 1:dim(prcp[[1]])[1]}
  stat.list2 <- mclapply(X=1:d,mc.preschedule=TRUE,mc.cores=mc.cores,FUN=function(i){  
    my.stats2 <- array(NA,c(n.sites,n.stats))
    cur.prcp <- matrix(prcp[[i]][my.range,])
    if (sum(is.infinite(cur.prcp))>0){
      idx.dummy <- which(is.infinite(cur.prcp),arr.ind = T)
      cur.prcp[idx.dummy[,1],idx.dummy[,2]] <- NA} # checkpoints to not have "Inf" before going ahead
    for (s in 1:n.stats) {
      my.stats2[,s] <- apply(cur.prcp,2,function(x){max(rollapply(data=x, FUN = mean, width = ma.options.rolling[s]))})
    }
    colnames(my.stats2) <- stat.names2
    return(my.stats2)
  })
  my.stats2 <- unlist(stat.list2)
  dim(my.stats2) <- c(n.sites,n.stats,d)    
  dimnames(my.stats2)[[2]] <- stat.names2
  return(my.stats2)
}


stat.names2 <-c(paste0(c(7,10),'.day.max'),paste0(c(1,3,6),'.month.max'),paste0(seq(1,10),'.year.max'))

years.sim <- wateryears.sim.new
years.obs <- wateryears.obs.new
obs.array <- sim.array <- array(NA,c(length(stat.names2),n.sites))

for (v in 1:n.sites){
  site.data <- as.matrix(prcp.site[,v])
  site.data.sim <- list()
  for (k in 1:num.iter){
    site.data.sim[[k]] <- as.matrix(prcp.site.sim[[num.iter]][,v])
  }
  
  obs.diagnostics.min <- prcp.max.extreme.diagnostics(prcp=list(site.data),my.years=years.weather)
  obs.array[,v] <-  obs.diagnostics.min.extreme <- as.matrix(obs.diagnostics.min[,,1])
  
  sim.diagnostics.min <- prcp.max.extreme.diagnostics(prcp=site.data.sim,my.years=years.sim)
  sim.array[,v] <- sim.diagnostics.min[1,,]
}

ww.mom <- 18; hh.mom <- 18
fig.stat.name1 <- paste0("Fig.n.days_n.yrs.max.precip_with_",basin.cnt,"_s",
                         num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)

par(mfrow=c(3,3),mar=c(5,6,4,2)+0.1)
nd.selec <- c(1,2,3,4,5,6,7,8,10)
for (nd in nd.selec){
  
  my.x.y.min.max <- c(min(sim.array[nd,],obs.array[nd,]),max(sim.array[nd,],obs.array[nd,]))
  if(nd==1){
    plot(sim.array[nd,],obs.array[nd,],
         xlim=my.x.y.min.max,ylim=my.x.y.min.max,
         pch=seq(0,14),col=rainbow(9),frame=F,
         main=paste0(stat.names2[nd]),
         sub='',cex.sub=1.5,
         xlab='Sim [WGEN]',ylab='',
         cex.axis=2.5,cex.lab=2,cex.main=3.5, cex=4,
         font.main=1)
    mtext('Obs', side=2, line=3)
    abline(a=0,b=1, col='gray75',lwd=1.5, lty=2)
    
    legend('bottomright',
           legend=c(paste('Site',seq(1,n.sites))),
           pch=seq(0,8),cex=2.5,col=rainbow(9),
           inset=c(0,0), xpd=TRUE,ncol=1,
           horiz = FALSE,
           title='Prcp [WY]')
  }
  else{
    plot(sim.array[nd,],obs.array[nd,],
         xlim=my.x.y.min.max,ylim=my.x.y.min.max,
         pch=seq(0,8),col=rainbow(9),frame=F,
         main=paste0(stat.names2[nd]),
         xlab='Sim [WGEN]',ylab='Obs',
         cex.axis=2.5,cex.lab=2,cex.main=3.5, cex=4,
         font.main=1)
    
    abline(a=0,b=1, col='gray75',lwd=1.5, lty=2)
  }
  
}
dev.off()


#================================================================================
## Plot: whiskers 4-D-panel important metrics #
#================================================================================
site.data <- prcp.site
site.data.sim <- prcp.site.sim
lst.annual.obs <- lst.lst.iqr.annual.sim <- list()
lst.annual.obs[[1]] <- aggregate(site.data,FUN=mean,by=list(wateryears.obs.new))[,-1]
lst.lst.iqr.annual.sim[[1]] <- lapply(1:num.iter,function(y){
  return(aggregate(site.data.sim[[y]],FUN=mean,by=list(wateryears.sim.new))[,-1])})
lst.annual.obs[[2]] <- aggregate(site.data,FUN=skewness,by=list(wateryears.obs.new))[,-1]
lst.lst.iqr.annual.sim[[2]] <- lapply(1:num.iter,function(y){
  return(aggregate(site.data.sim[[y]],FUN=skewness,by=list(wateryears.sim.new))[,-1])})
lst.annual.obs[[3]] <- aggregate(site.data,FUN=sd,by=list(wateryears.obs.new))[,-1]
lst.lst.iqr.annual.sim[[3]] <- lapply(1:num.iter,function(y){
  return(aggregate(site.data.sim[[y]],FUN=sd,by=list(wateryears.sim.new))[,-1])})
lst.annual.obs[[4]] <- aggregate(site.data,FUN=quantile,0.99,by=list(wateryears.obs.new))[,-1]
lst.lst.iqr.annual.sim[[4]] <- lapply(1:num.iter,function(y){
  return(aggregate(site.data.sim[[y]],FUN=quantile,0.99,by=list(wateryears.sim.new))[,-1])})

ww.mom <- 12; hh.mom <- 12
fig.stat.name1 <- paste0("Fig.prcp_WY.whiskers.mean.iqr.std.99th.",basin.cnt,"_s",
                         num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(2,2),mar=c(5.1, 6.1, 4.1, 2.1))
labs.metrics <- c('Mean','Skewness','Stdev','99th')

for (k in 1:length(labs.metrics)){
  annual.obs <- lst.annual.obs[[k]]
  lst.iqr.annual.sim <- lst.lst.iqr.annual.sim[[k]]
  my.sim.y <- lst.iqr.annual.sim[[num.iter]]
  my.obs.y <- annual.obs
  
  my1.x0 <- apply(my.sim.y,2,FUN=median)
  my1.y0 <- apply(my.obs.y,2,FUN=quantile,0.25)
  my1.y1 <- apply(my.obs.y,2,FUN=quantile,0.75)
  
  my2.x0 <- my2.x1 <- apply(my.obs.y,2,FUN=median)
  my2.y0 <- apply(my.sim.y,2,FUN=quantile,0.25)
  my2.y1 <- apply(my.sim.y,2,FUN=quantile,0.75)
  x_y_lim <- c(min(c(my1.y0,my2.y0)),max(c(my1.y1,my2.y1)))
  
  plot(my1.x0, my2.x0, xlim=x_y_lim, ylim=x_y_lim,
       pch=1, cex=2, frame=F,col=NULL,
       xlab='',ylab='',
       cex.axis=1.5,
       main='',sub='(25th,median,75th)',col.sub=alpha('blue',0.5),
       cex.lab=1.25, cex.main=2, font.main=1,cex.sub=1.15)
  mtext(paste0('Sim [WGEN, ',total.num.yrs,' years]'),1,col='gray50',line=2.25)
  mtext('Obs',2,col='red',line=2.25)
  
  # abline(v=my1.x0,col=alpha('gray',0.8),lwd=0.25,lty=4)
  # abline(h=my2.x0,col=alpha('gray',0.8),lwd=0.25,lty=4)
  arrows(x0=my2.y0,
         y0=my2.x0, x1=my2.y1, y1=my2.x0, code=3, angle=90, length=0.025, col=alpha("gray40",0.75), lwd=0.5)
  arrows(x0=my1.x0,
         y0=my1.y0, x1=my1.x0, y1=my1.y1, code=3, angle=90, length=0.025, col=alpha("red",0.75), lwd=0.5)
  abline(0,1,lwd=1.5,lty=2)
  points(my1.x0, my2.x0, xlim=x_y_lim, ylim=x_y_lim, pch=19, 
         cex=1.5,col=alpha("gray40",0.65))
  legend("topleft", legend=paste0(labs.metrics[k]),pch=19,cex=1.5,
         col='gray60',title ='(At-Site)',bty = 'n',
         horiz = FALSE,inset=0.01)
}
dev.off()

#

#================================================================================
## Plot: Set of metrics as boxplot for PBIAS across sites#
#================================================================================
# for Precipitation #
data.obs <- prcp.site
data.sim <- prcp.site.sim[[1]]


long.dates.sim.hist <- cbind(YEAR.in.long.trace,MONTH.in.long.trace,DAY.in.long.trace)
long.dates.sim.vec <- long.dates.sim.hist[,1]*10000+long.dates.sim.hist[,2]*100+long.dates.sim.hist[,3]
long.years.vec <- YEAR.in.long.trace
long.months.vec <- MONTH.in.long.trace

wateryears.sim.new <- (long.years.vec+1)*(long.months.vec>=10) + (long.years.vec)*(long.months.vec<10)
wateryears.obs.new <- (years.weather+1)*(months.weather>=10) + (years.weather)*(months.weather<10)

wy.data.obs <- apply(data.obs,2,function(x){
  aggregate(x,FUN=sum,by=list(wateryears.obs.new),na.rm=T)[,2]})
wy.data.obs <- apply(wy.data.obs,2,function(x){
  return(x[2:(length(x)-1)])})#drop first and last to only keep full water years

wy.data.sim <- apply(data.sim,2,function(x){
  aggregate(x,FUN=sum,by=list(wateryears.sim.new),na.rm=T)[,2]})
wy.data.sim <- apply(wy.data.sim,2,function(x){
  return(x[2:(length(x)-1)])})#drop first and last to only keep full water years


##Metrics for PBIAS
nsites <- n.sites
metrics <- cbind(
  'SD'=sapply(1:nsites,function(x){return(100*(sd(wy.data.sim[,x])-sd(wy.data.obs[,x]))/sd(wy.data.obs[,x]))}),
  'MEAN'=sapply(1:nsites,function(x){return(100*(mean(wy.data.sim[,x])-mean(wy.data.obs[,x]))/mean(wy.data.obs[,x]))}),
  'Q99'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.99)-quantile(wy.data.obs[,x],0.99))/quantile(wy.data.obs[,x],0.99))}),
  'Q90'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.9)-quantile(wy.data.obs[,x],0.9))/quantile(wy.data.obs[,x],0.9))}),
  'MEDIAN'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.5)-quantile(wy.data.obs[,x],0.5))/quantile(wy.data.obs[,x],0.5))}),
  'Q05'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.05)-quantile(wy.data.obs[,x],0.05))/quantile(wy.data.obs[,x],0.05))}),
  'MIN'=sapply(1:nsites,function(x){return(100*(min(wy.data.sim[,x])-min(wy.data.obs[,x]))/min(wy.data.obs[,x]))}),
  'MIN2'=sapply(1:nsites,function(x){return(100*(min(rollmean(wy.data.sim[,x],2))-min(rollmean(wy.data.obs[,x],2)))/min(rollmean(wy.data.obs[,x],2)))}),
  'MIN3'=sapply(1:nsites,function(x){return(100*(min(rollmean(wy.data.sim[,x],3))-min(rollmean(wy.data.obs[,x],3)))/min(rollmean(wy.data.obs[,x],3)))}),
  'MIN5'=sapply(1:nsites,function(x){return(100*(min(rollmean(wy.data.sim[,x],5))-min(rollmean(wy.data.obs[,x],5)))/min(rollmean(wy.data.obs[,x],5)))}),
  'MAX'=sapply(1:nsites,function(x){return(100*(max(wy.data.sim[,x])-max(wy.data.obs[,x]))/max(wy.data.obs[,x]))}),
  'MAX2'=sapply(1:nsites,function(x){return(100*(max(rollmean(wy.data.sim[,x],2))-max(rollmean(wy.data.obs[,x],2)))/max(rollmean(wy.data.obs[,x],2)))}),
  'MAX3'=sapply(1:nsites,function(x){return(100*(max(rollmean(wy.data.sim[,x],3))-max(rollmean(wy.data.obs[,x],3)))/max(rollmean(wy.data.obs[,x],3)))}),
  'MAX5'=sapply(1:nsites,function(x){return(100*(max(rollmean(wy.data.sim[,x],5))-max(rollmean(wy.data.obs[,x],5)))/max(rollmean(wy.data.obs[,x],5)))})
)

fct.spell.stats <- function(prcp.file){
  prcp.thresh <- 0    #trace precipitation threshold (0 inches)
  apply(prcp.file,2,function(x) {
    y <- x[!is.na(x)]
    nn <- length(y)
    wet.dry <- y
    wet.dry[wet.dry<=prcp.thresh] <- 0
    wet.dry[wet.dry>prcp.thresh] <- 1
    first.of.run <- c(1,which(diff(wet.dry)!=0) + 1)
    last.of.run <- c(which(diff(wet.dry)!=0),nn)
    run.length <- last.of.run-first.of.run + 1
    run.type <- wet.dry[first.of.run]
    wet.dry.runs <- data.frame('type'=run.type,'first'=first.of.run,'last'=last.of.run,'length'=run.length)
    result <- c(
      mean(wet.dry.runs$length[wet.dry.runs$type==1]),
      max(wet.dry.runs$length[wet.dry.runs$type==1]),
      mean(wet.dry.runs$length[wet.dry.runs$type==0]),
      max(wet.dry.runs$length[wet.dry.runs$type==0])
    )
    return(result)
  })
}

spell.stats.metrics <- t(100*(fct.spell.stats(data.sim) - fct.spell.stats(data.obs))/fct.spell.stats(data.obs))
pBias.full.set.metrics <- cbind(metrics,spell.stats.metrics)
colnames(pBias.full.set.metrics) <- c(colnames(metrics),c('MEAN.W.SPELL','MAX.W.SPELL','MEAN.D.SPELL','MAX.D.SPELL'))
rownames(pBias.full.set.metrics) <- NULL

my.median.values <- paste0(round(apply(pBias.full.set.metrics,2,median),0),"%")
# Plot: PBIAS for Obs vs. Baseline
ww.mom <- 8; hh.mom <- 6
fig.stat.name1 <- paste0("Fig.boxplot_pbias_precip_Obs.vs.baseline.",basin.cnt,"_s",
                         num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)
boxplot(pBias.full.set.metrics,boxwex=.5,main=paste0('Baseline [WGEN: ',total.num.yrs,'-year] vs. Obs'),
        col=alpha('blue',0.5),outline=F,frame=F,cex.axis=1.25,cex.lab=1.25,cex.main=1.25,font.main=1,
        xlab='Precipitation [water year]',ylab='PBIAS[%]',xaxt = "n",ylim=c(-80,80),
        border=alpha('blue',0.5),whisklty = 1,medcol = alpha('red',0.5))
abline(v=1:dim(pBias.full.set.metrics)[2],col=alpha('gray',0.25),lty=2)
abline(h=0,col=alpha('red',0.15),lwd=4)
axis(side = 1, labels = FALSE)
text(x = 1:dim(pBias.full.set.metrics)[2],
     y = par("usr")[3],labels = colnames(pBias.full.set.metrics),
     xpd = NA,adj=c(rep(0,dim(metrics)[2]),rep(0.1,dim(spell.stats.metrics)[2])),
     ## Rotate the labels by 35 degrees.
     srt = 90,cex = 1, col='gray50')

text(x = 1:dim(pBias.full.set.metrics)[2],
     y = 87.5,labels = my.median.values,
     xpd = NA,adj=0,
     ## Rotate the labels by 35 degrees.
     srt = 45,cex = 0.75, col=alpha('red',0.5))
text(x = 0,
     y = 87.5,labels = 'median',
     xpd = NA,adj=0,
     ## Rotate the labels by 35 degrees.
     srt = 90,cex = 0.65, col=alpha('red',0.5))

legend('topleft',paste0('Sites: 1-',nsites),fill=alpha('blue',0.10),cex=1.15,bty = "n")
dev.off()



#================================================================================
## Plot: CDF of precipitation per site
#================================================================================
#// CDF
#get annual total prcp from observation
ww.mom <- 18; hh.mom <- 14
fig.stat.name1 <- paste0("Fig.EP_precip_Obs.vs.baseline.",basin.cnt,"_s",
                         num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(3,4),mar=c(5.1, 6.1, 4.1, 2.1))
for (v in 1:n.sites){
  annual.prcp <- aggregate(prcp.site[,v],FUN=sum,by=list(wateryears.obs.new))[,2]
  annual.prcp <- annual.prcp[2:(length(annual.prcp)-1)]   #drop first and last to only keep full water years
  annual.prcp.boot.sort <- array(NA,c(length(annual.prcp),1000))
  annual.prcp.sd <- array(NA,1000)
  for (i in 1:1000) {
    annual.prcp.boot.sort[,i] <- sort(sample(annual.prcp,size=length(annual.prcp),replace=TRUE))
    annual.prcp.sd[i] <- sd(sample(annual.prcp,size=length(annual.prcp),replace=TRUE))
  }
  boot.95 <- apply(annual.prcp.boot.sort,1,FUN=quantile,.95)
  boot.05 <- apply(annual.prcp.boot.sort,1,FUN=quantile,.05)
  #get annual total prcp from simulation
  prcp.basin.sim <- prcp.site.sim[[1]][,v]
  annual.prcp.sim <- aggregate(prcp.basin.sim,FUN=sum,by=list(wateryears.sim.new),na.rm=T)[,2]
  annual.prcp.sim <- annual.prcp.sim[2:(length(annual.prcp.sim)-1)] #drop first and last to only keep full water years
  
  
  y.hat <- annual.prcp
  EP.hat <- (1:length(y.hat))/(length(y.hat)+1)
  y <- annual.prcp.sim
  EP <- (1:length(y))/(length(y)+1)
  plot(EP.hat,sort(y.hat),log="y",type='l',frame=F,main=paste0('Site ',v),
       col="red",lwd=2.5,xlab='EP',ylab='Precipitation [WY Total]',cex.axis=1.75,cex.main=2.5,cex.lab=1.65, font.main=2)
  lines(EP,sort(y),lty=1, col='gray50',lwd=3)
  polygon(x=c(EP.hat,rev(EP.hat)),y=c(boot.05,rev(boot.95)),col=alpha('red',0.15),border=NA)
  legend('topleft',c('Obs [1000 boots]',paste0('Sim [WGEN:',total.num.yrs,'-yr]')),horiz = F,cex=1.5,
         fill=c('red','gray50'),bty = 'n')
}
dev.off()
#


#================================================================================
## Plot: SD of precipitation per site
#================================================================================
#get annual total prcp from observation
ww.mom <- 16; hh.mom <- 14
fig.stat.name1 <- paste0("Fig.SD_boxplot.precip_Obs.vs.baseline.",basin.cnt,"_s",
                         num.states,"_with_",num.iter,"_ens_all.sites.",n.sites,".png")
png(paste0("./Figures/",fig.stat.name1),width=ww.mom,height=hh.mom,units="in",res=300)
par(mfrow=c(3,4),mar=c(5.1, 6.1, 4.1, 2.1))
for (v in 1:n.sites){
  annual.prcp <- aggregate(prcp.site[,v],FUN=sum,by=list(wateryears.obs.new))[,2]
  annual.prcp <- annual.prcp[2:(length(annual.prcp)-1)]   #drop first and last to only keep full water years
  annual.prcp.sd <- array(NA,1000)
  for (i in 1:1000) {
    annual.prcp.sd[i] <- sd(sample(annual.prcp,size=length(annual.prcp),replace=TRUE))
  }
  #get annual total prcp from simulation
  prcp.basin.sim <- prcp.site.sim[[1]][,v]
  annual.prcp.sim <- aggregate(prcp.basin.sim,FUN=sum,by=list(wateryears.sim.new),na.rm=T)[,2]
  annual.prcp.sim <- annual.prcp.sim[2:(length(annual.prcp.sim)-1)] #drop first and last to only keep full water years
  
  
  boxplot(annual.prcp.sd,col=alpha('red',0.2),boxwex=0.5,
          frame=F,main=paste0('Site ',v),
          xlab='',ylab='Annual SD [WY Prcp Total]',cex.axis=2,cex.main=2.5,cex.lab=1.65, font.main=1)
  points(sd(annual.prcp),col="red",pch=16,cex=3)
  points(sd(annual.prcp.sim),col="gray50",pch=16,cex=3)
  legend('topleft',c('Obs [1000 boots]',paste0('Sim [WGEN:',total.num.yrs,'-yr]')),horiz = F,cex=1.25,
         col=c('red','gray50'),pch=16,bty = 'n')
}
dev.off()
#

# THE END #
