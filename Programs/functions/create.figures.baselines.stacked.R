
create.figures.baselines.stacked <- function(scenario=selected_scenario){
  # stacked version #
  
  #scenario = the row in ClimateChangeScenarios.csv for which to plot results
  
  # ------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------- #
  # A set of 8 figures stacked to show diagnostics for WGEN baseline (thousands-yr): precipitation and temperature #
  # ------------------------------------------------------------------------------- #
  # ------------------------------------------------------------------------------- #
  # including:
  #- Precipitation:
  # 1) extreme precipitation, wet-spells, flood stats
  # 2) dry-spells, drought stats
  # 3) seasonality, spatial correlation for extremes
  # 4) percent bias stats
  #- Temperature:
  # 1) spatial pattern
  # 2) heat waves, heat stress
  # 3) cold waves, cold stress
  
  #location of obs weather data
  load(path.to.processed.data.meteohydro) #load in weather data

  #preset labels for figures below
  identical.dates.idx <- dates.weather%in%dates.synoptics
  dates.weather <- dates.weather[identical.dates.idx]
  years.weather <- as.numeric(format(dates.weather,'%Y'))
  months.weather <- as.numeric(format(dates.weather,'%m'))
  
  label1 = paste0('Obs [',format(as.Date(range(dates.weather)[1]),'%Y'),
                  '-',format(as.Date(range(dates.weather)[2]),'%Y'),']')
  label2 = 'Sim (WGEN: baseline)'
  
  # use.non_param.WRs <- TRUE #TRUE for non-parametric, FALSE for parametric simulated WRs
  # num.states <- c(10) # K WRs
  # qq <- 0.99              # quantile to anchor the CC-scaling, and to split gamma and GPD
  # num.iter <- 1 # Number of ensemble members generated using WR-SWG
  # 
  # # directory to simulation files
  # dir.to.sim.files <- "./Data/simulated.data.files/WGEN.out/"
  
  ##// weather data and synoptic dates ---##
  prcp.site.data <- prcp.site[identical.dates.idx,]
  tmin.site.data <- tmin.site[identical.dates.idx,]
  tmax.site.data <- tmax.site[identical.dates.idx,]
  tmean.site.data <- (tmin.site.data+tmax.site.data)/2
  
  n.sites <- dim(prcp.site)[2] # Number of gridded points for precipitation
  
  #select the scenario for plotting
  cur.tc.max <- change.list$tc.max[scenario]
  cur.tc.min <- change.list$tc.min[scenario]
  cur.pccc <- change.list$pccc[scenario]
  cur.pmuc <- change.list$pmuc[scenario]

  cur.jitter <- to.jitter
  
  ##// sim files ---##
  {
    simulated.file.run.model.saved <- paste0(".tmax.",cur.tc.max,".tmin.",cur.tc.min,"_p.CC.scale.",cur.pccc,"_p.mu.scale.",cur.pmuc,"_num.year.",number.years.long,"_with.",num.iter)

    prcp.site.sim_sfx <- "prcp.site.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    prcp.site.data.sim <- prcp.site.sim[[1]] # becase there was only one iteration
    
    prcp.site.sim_sfx <- "tmin.site.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    tmin.site.data.sim <- tmin.site.sim[[1]]
    
    prcp.site.sim_sfx <- "tmax.site.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    tmax.site.data.sim <- tmax.site.sim[[1]]
    
    tmean.site.data.sim <- (tmax.site.data.sim+tmin.site.data.sim)/2
    
    prcp.site.sim_sfx <- "dates.sim"
    load(paste0(dir.to.sim.files,"/",prcp.site.sim_sfx,simulated.file.run.model.saved,".RData"))
    #dates.sim
    
    months.sim <- as.numeric(format(dates.sim,"%m"))
    years.sim <- as.numeric(format(dates.sim,"%Y"))
  }
  
  
  # water years dates/values #
  {
    # generating synthetic dates/years #
    a <- rle(years.sim)
    a.freq <- a$lengths
    a.val <- a$values
    yr=1
    yr.vec <- c()
    for (k in 1:length(a.freq)){
      yr.vec <- append(yr.vec,rep(yr,a.freq[k]))
      yr <- yr+1
    }
    long.years.vec <- yr.vec
    long.months.vec <- as.numeric(months.sim)
    total.num.yrs <- max(long.years.vec)
    
    # water years dates #
    wateryears.weather <- (years.weather+1)*(months.weather>=10) + (years.weather)*(months.weather<10)
    wateryears.obs <- wateryears.weather
    wateryears.sim <- (long.years.vec+1)*(long.months.vec>=10) + (long.years.vec)*(long.months.vec<10)
    
    # water years precipitation total #
    wy.data.obs <- apply(prcp.site.data,2,function(x){
      aggregate(x,FUN=sum,by=list(wateryears.obs),na.rm=T)[,2]})
    wy.data.obs <- apply(wy.data.obs,2,function(x){
      return(x[2:(length(x)-1)])})#drop first and last to only keep full water years
    
    wy.data.sim <- apply(prcp.site.data.sim,2,function(x){
      aggregate(x,FUN=sum,by=list(wateryears.sim),na.rm=T)[,2]})
    wy.data.sim <- apply(wy.data.sim,2,function(x){
      return(x[2:(length(x)-1)])})#drop first and last to only keep full water years
    
    # water years temperature average #
    wy.tmean.data.obs <- apply(tmean.site.data,2,function(x){
      aggregate(x,FUN=mean,by=list(wateryears.obs),na.rm=T)[,2]})
    wy.tmean.data.obs <- apply(wy.tmean.data.obs,2,function(x){
      return(x[2:(length(x)-1)])})#drop first and last to only keep full water years
    
    wy.tmean.data.sim <- apply(tmean.site.data.sim,2,function(x){
      aggregate(x,FUN=mean,by=list(wateryears.sim),na.rm=T)[,2]})
    wy.tmean.data.sim <- apply(wy.tmean.data.sim,2,function(x){
      return(x[2:(length(x)-1)])})#drop first and last to only keep full water years
  }
  
  
  # average precipitation & temperature at-basin #
  {
    # average precipitation #
    mean.prcp.site <- as.matrix(apply(prcp.site.data,1,mean))
    mean.prcp.site.sim <- as.matrix(apply(prcp.site.data.sim,1,mean))
    
    # average temperature #
    mean.tmean.site <- as.matrix(apply(tmean.site.data,1,mean))
    mean.tmean.site.sim <- as.matrix(apply(tmean.site.data.sim,1,mean))
    
    # average wy precipitation #
    mean.wy.prcp.site <- as.matrix(apply(wy.data.obs,1,mean))
    mean.wy.prcp.site.sim <- as.matrix(apply(wy.data.sim,1,mean))
    
    # average wy temperature #
    mean.wy.tmean.site <- as.matrix(apply(wy.tmean.data.obs,1,mean))
    mean.wy.tmean.site.sim <- as.matrix(apply(wy.tmean.data.sim,1,mean))
    
  }
  
  # Figure 1:
  # 1) precipitation extremes, wet-spells, flood stats
  {
    # label1 <- 'Obs [1950-2013]' # change wrt your observed dates/meteohydro
    # label2 <- 'Sim (WGEN: baseline)'
    label11 <- 'Obs [GEV]'
    label12 <- 'Obs [GPD]'
    
    file.name.fig <- paste0("lines_prcp.GEV.GPD.NEP_dur.maxima_wet.spells_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 14; hh.mom <- 14
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfrow=c(2,2),mai=c(1,1,1,0.5))
    
    # (a)
    #----------------------------------------- #
    ## -- GPD threshold -- ##
    #----------------------------------------- #
    my.u.cluster <- 0.99
    num.per.block <- 365
    #my.ret.per <- c(2,5,10,20,50,100,500,1000)
    my.ret.per <- seq(1.001,10000)
    # obs annual max
    PRCP.input <- mean.prcp.site
    years.input <- wateryears.weather
    site.data.info <- data.frame(cbind(PRCP.input,years.input))
    colnames(site.data.info) <- c("PRCP","YEARS")
    site.data.info$time <- 1:nrow(site.data.info)
    my.u.threshold <- quantile(site.data.info$PRCP[site.data.info$PRCP>0],probs=my.u.cluster)
    site.data.info$lambda <- my.u.threshold
    idx.pot <- which(site.data.info$PRCP>my.u.threshold,arr.ind = T)#extra
    years.p.max <- site.data.info$YEARS[idx.pot]#extra
    lambda.occurrence_counts <- mean(table(years.p.max))#extra
    
    #PDS is the data frame with annual maxima
    PDS <- site.data.info
    my.thre <- unique(PDS$lambda)
    #fit competing stationary models based on PRCP
    # static.gpd <- ismev::gpd.fit(PDS$PRCP,
    #                              show=TRUE,threshold=my.thre,
    #                              npy=num.per.block,
    #                              sigl=NULL,ydat=NULL,method="SANN")
    # gpd.diag(static.gpd)
    # mrl.plot(PDS$PRCP,umin = min(PDS$PRCP), umax = max(PDS$PRCP)- 0.1,
    #          conf = 0.95, nint = 100)
    static.gpd <- fevd(PDS$PRCP,threshold = my.thre,type='GP')
    my.rls <- return.level(static.gpd,do.ci=TRUE,
                           return.period = my.ret.per)
    my.rl.vals <- as.matrix(my.rls)
    rownames(my.rl.vals) <- NULL
    
    
    #----------------------------------------- #
    ## -- GEV, NEP -- ##
    #----------------------------------------- #
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
      # ymax <- max(dat, ymax, q+ 1.96 * sqrt(v))
      # ymin <- min(dat, q- 1.96 * sqrt(v))
      ymax <- max(dat, ymax, q+sqrt(v))
      ymin <- min(dat, q-sqrt(v))
      plot(-1/log(f), q, log = "x", type = "n", xlim = c(0.1, 5000), panel.first=grid(equilogs=TRUE,lty=1,col="gray80",lwd=0.5),
           ylim = c(ymin, ymax), xlab = "Return Period [year]",font.main = 1, 
           ylab = "Precipitation [return level]",xaxt = 'n',
           cex=2,main='',
           cex.axis=1.75,cex.lab=1.5,cex.main=2, frame=FALSE)
      # newx <- -1/log(f)
      # newy <- q
      # polygon(c(rev(newx), newx), c(rev(newy),newy), density=20, col = rgb(0,1,0,0.3), border = NA)
      
      lines(-1/log(f), q, col='red',lty=1, lwd=1.5)
      lines(-1/log(f), q + 1.96 * sqrt(v), col = alpha("red",.4), lty=2, lwd=1)
      lines(-1/log(f), q - 1.96 * sqrt(v), col = alpha("red",.4), lty=2, lwd=1)
      points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat),col='red', pch=4)
    }
    
    fct.prcp.sim.annmax <- function(prcp,my.years,my.range=NULL,mc.cores=1){
      uniq.years <- unique(my.years)
      d <- length(prcp)
      n.sites <- dim(prcp[[1]])[2]
      if(is.null(my.range)) {my.range <- 1:dim(prcp[[1]])[1]}
      stat.rp.rl.list <- mclapply(X=1:d,mc.preschedule=TRUE,mc.cores=mc.cores,FUN=function(i){  
        cur.prcp <- as.matrix(prcp[[i]][my.range,])
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
    psite.obs.amax <- aggregate(mean.prcp.site,FUN=max,by=list(wateryears.obs))[,2]
    # sim annual max
    prcp.sim.annmax0 <- fct.prcp.sim.annmax(prcp=list(mean.prcp.site.sim),my.years=wateryears.sim)
    prcp.sim.annmax <- apply(prcp.sim.annmax0,2,FUN=as.vector)
    ymax <- max(c(psite.obs.amax,prcp.sim.annmax),na.rm=T)
    
    require(ismev)
    my.gev.fit1 <- gev.fit(psite.obs.amax,show = F)
    my.gev.rl(my.gev.fit1$mle, my.gev.fit1$cov, my.gev.fit1$data,ymax,site.num)
    dat <- prcp.sim.annmax
    dat <- dat[!is.na(dat)]
    points(-1/log((1:length(dat))/(length(dat) + 1)), sort(dat),
           col="black",cex=0.75,pch=16)
    
    lines(my.ret.per,my.rl.vals[,2],col='blue',lwd=2,lty=1)
    lines(my.ret.per,my.rl.vals[,1],col='blue',lwd=1,lty=3)
    lines(my.ret.per,my.rl.vals[,3],col='blue',lwd=1,lty=3)
    
    myTicks = axTicks(1)
    axis(1, at = myTicks, labels = formatC(myTicks, format = 'd'),cex.axis=1.75)
    legend("topleft", legend=c(label1,label11,label12,label2),
           col=c("red", "red",'blue',"black"), 
           lty=c(NA,1,1,NA),cex=1.5,pch=c(4,NA,NA,16),
           lwd=c(NA,1.5,2,NA),
           bty = "n",horiz = FALSE,
           title='at-basin',title.col = 'purple')

    
    # (b)
    #----------------------------------------- #
    ## -- correlation between extremes -- ##
    #----------------------------------------- #
    resampled.date.sim_sfx <- "/resampled.date.sim.sim"
    load(paste0(dir.to.sim.files,resampled.date.sim_sfx,simulated.file.run.model.saved,".RData"))
    resampled.dates <- match(as.Date(resampled.date.sim[[1]]),dates.weather)
    quantile.nz <- function(x,p=qq){quantile(x[x!=0],p,na.rm=T)}
    thshd.prcp <- apply(prcp.site.data,2,FUN=quantile.nz)
    pot.idx <- sapply(1:n.sites,function(x,v){(prcp.site.data[,x]>v[x])
    },v=thshd.prcp)
    pot.idx.boot <- pot.idx[resampled.dates,]
    prcp.site.boot <- prcp.site.data[resampled.dates,]
    prcp.site.data.sim.ex <- prcp.site.data.sim
    prcp.site.data.sim.ex <- prcp.site.data.sim.ex[pot.idx.boot]
    prcp.site.boot.ex.specific <- prcp.site.boot[pot.idx.boot]
    vall <- round(cor(prcp.site.boot.ex.specific,prcp.site.data.sim.ex),2)
    
    plot(prcp.site.data.sim.ex,prcp.site.boot.ex.specific,
         xlim=c(min(prcp.site.data.sim.ex),
                max(prcp.site.data.sim.ex)),
         ylim=c(min(prcp.site.data.sim.ex),
                max(prcp.site.data.sim.ex)),
         ylab=paste0(label1,'(boot), max=',round(max(prcp.site.boot.ex.specific),1)),
         xlab=paste0(label2,', max=',round(max(prcp.site.data.sim.ex),1)),
         cex.axis=1.75,cex.lab=1.5,cex.main=2,font.main=1,cex.sub=1.15,
         pch=19,col=NULL,cex=1.5,
         frame=T,main=paste0('Precipitation Extremes'),
         sub=paste0('cor=',vall),col.sub='purple')
    abline(a=0,b=1, col='gray75',lwd=1.5, lty=2)
    points(prcp.site.data.sim.ex,prcp.site.boot.ex.specific,
           pch=19,col=alpha('darkblue',0.25),cex=1.5)
    legend("topleft", legend=c('peak-over-threshold events'),
           col=c('darkblue'),cex=1.5,pch=c(19),
           bty = "n",horiz = FALSE,
           title='at-site',title.col = 'purple')
    
    # (c)
    #----------------------------------------- #
    ## -- short/long-term extreme -- ##
    #----------------------------------------- #
    # less than a year: n-days/years precip max
    prcp.max.extreme.diagnostics <- function(prcp,my.years,my.range=NULL,mc.cores=1) {
      d <- length(prcp)
      n.sites <- dim(prcp[[1]])[2]
      ma.options.rolling <- c(7,10,30,90,180,seq(1,1)*365) # n-day moving averages
      stat.names2 <- c(paste0(c(7,10),'.day.max'),
                       paste0(c(1,3,6),'.month.max'),
                       paste0(seq(1,1),'.year.max'))
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
          my.stats2[,s] <- apply(cur.prcp,2,function(x){max(rollapply(data=x, FUN = mean, width = ma.options.rolling[s],na.rm=T))})
        }
        colnames(my.stats2) <- stat.names2
        return(my.stats2)
      })
      my.stats2 <- unlist(stat.list2)
      dim(my.stats2) <- c(n.sites,n.stats,d)    
      dimnames(my.stats2)[[2]] <- stat.names2
      return(my.stats2)
    }
    
    stat.names2 <-c(paste0(c(7,10),'-day max'),
                    paste0(c(1,3,6),'-month max'),
                    paste0(seq(1,1),'-year max'))
    years.sim <- wateryears.sim
    years.obs <- wateryears.obs
    obs.diagnostics.min <- prcp.max.extreme.diagnostics(prcp=list(mean.prcp.site),my.years=years.weather)
    obs.array <-  obs.diagnostics.min.extreme <- as.matrix(obs.diagnostics.min[,,1])
    sim.diagnostics.min <- prcp.max.extreme.diagnostics(prcp=list(mean.prcp.site.sim),my.years=years.sim)
    sim.array <- sim.diagnostics.min[1,,]
    my.pch <- 15:20
    nd.selec <- c(1,2,3,4,5,6,7,8,10)
    my.x.y.min.max <- c(min(sim.array,obs.array),max(sim.array,obs.array))
    
    plot(sim.array,obs.array,
         xlim=my.x.y.min.max,ylim=my.x.y.min.max,
         pch=my.pch,col=rainbow(length(nd.selec)),frame=F,
         main='Short/Long-term Maxima',sub='',cex.sub=1.5,xlab=label2,ylab='',
         cex.axis=1.75,cex.lab=1.5,cex.main=2,cex=3.5,font.main=1)
    mtext(label1, side=2, line=3, cex=1.5)
    abline(a=0,b=1, col='gray75',lwd=1.5, lty=2)
    legend('topleft',legend=paste(stat.names2),pch=my.pch,cex=1.5,col=rainbow(length(nd.selec)),
           inset=c(0.025,0), xpd=TRUE,horiz = FALSE,title='at-basin',title.col = 'purple')
    
    
    # (d)
    #----------------------------------------- #
    ## -- wet spells -- ##
    #----------------------------------------- #
    spells.labels <- c('mean.wet.spell','max.wet.spell',
                       'mean.dry.spell','max.dry.spell')
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
    my.sim.col <- c('black')
    my.obs.col <- c('red')
    obs.spell.stats <- fct.spell.stats(prcp.site.data)
    sim.spell.stats <- fct.spell.stats(prcp.site.data.sim)
    x1 <- obs.spell.stats[1,]
    x2 <- sim.spell.stats[1,]
    x11 <- obs.spell.stats[2,]
    x22 <- sim.spell.stats[2,]
    
    boxplot(x1,horizontal = T,main='',
            col=alpha(my.obs.col,.2),
            boxcol=my.obs.col,medcol=my.obs.col,
            whiskcol=my.obs.col,staplecol=my.obs.col,outcol=my.obs.col,
            boxwex=.35,xlab='Mean length (days)',
            cex.axis=1.75,cex.lab=1.5,cex.main=2,font.main=1)
    boxplot(x2,add=T,horizontal = T,
            col=alpha(my.sim.col,.2),
            boxcol=my.sim.col,medcol=my.sim.col,
            whiskcol=my.sim.col,staplecol=my.sim.col,outcol=my.sim.col,
            boxwex=.2,
            cex.axis=1.75,cex.lab=1.5,cex.main=2,font.main=1)
    text(mean(c(x1,x2)),0.5,cex=1.5,col='blue',
         paste0('max. Obs=',round(max(x11),1),
                ', max. Sim=',round(max(x22),1)))
    mtext('Wet Spells',3,cex=1.75,line = 1)
    legend("topleft", legend=c(label1,label2),
           fill=c("red", alpha(my.sim.col,.5)),cex=1.5,
           bty = "n",horiz = FALSE,
           title='at-site',title.col = 'purple')
    
    dev.off()
  }
  
  
  # Figure 2:
  # 2) precipitation cdf, mean, standard variation (yearly, year-to-year)
  {
    
    file.name.fig <- paste0("lines_prcp.mean_cdf_sd_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 14; hh.mom <- 14
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfcol=c(2,2),mai=c(1,1,1,0.5))
    
    
    # (a)
    #------------------------------------------ #
    ## -- cdf of precipitation distribution -- ##
    #------------------------------------------ #
    annual.prcp <- aggregate(mean.prcp.site,FUN=sum,by=list(wateryears.obs))[,2]
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
    prcp.basin.sim <- mean.prcp.site.sim
    annual.prcp.sim <- aggregate(prcp.basin.sim,FUN=sum,by=list(wateryears.sim),na.rm=T)[,2]
    annual.prcp.sim <- annual.prcp.sim[2:(length(annual.prcp.sim)-1)] #drop first and last to only keep full water years
    
    
    y.hat <- annual.prcp
    EP.hat <- (1:length(y.hat))/(length(y.hat)+1)
    y <- annual.prcp.sim
    EP <- (1:length(y))/(length(y)+1)
    plot(EP.hat,sort(y.hat),log="y",type='l',frame=F,
         col="red",lwd=5,xlab='Exceedance Probability',ylab='Precipitation [WY]',
         cex.axis=2,cex.main=4,cex.lab=2, font.main=1)
    lines(EP,sort(y),lty=1, col='gray50',lwd=6)
    polygon(x=c(EP.hat,rev(EP.hat)),y=c(boot.05,rev(boot.95)),
            col=alpha('red',0.15),border=NA)
    legend('topleft',c(paste0(label1,' (boot)'),paste0(label2)),
           horiz = F,
           # cex=1.75,
           cex=2.25,
           # title=paste0('at-basin'),
           # title.col = 'purple',title.cex = 2.5,
           fill=c('red','black'),bty = 'n')
    
    # (b)
    #----------------------------------------- #
    ## -- year-to-year SD of precipitation -- ##
    #----------------------------------------- #
    annual.prcp <- aggregate(mean.prcp.site,FUN=sum,by=list(wateryears.obs))[,2]
    annual.prcp <- annual.prcp[2:(length(annual.prcp)-1)]   #drop first and last to only keep full water years
    annual.prcp.sd <- array(NA,1000)
    for (i in 1:1000) {
      annual.prcp.sd[i] <- sd(sample(annual.prcp,size=length(annual.prcp),replace=TRUE))
    }
    #get annual total prcp from simulation
    prcp.basin.sim <- mean.prcp.site.sim
    annual.prcp.sim <- aggregate(prcp.basin.sim,FUN=sum,by=list(wateryears.sim),na.rm=T)[,2]
    annual.prcp.sim <- annual.prcp.sim[2:(length(annual.prcp.sim)-1)] #drop first and last to only keep full water years
    
    boxplot(annual.prcp.sd,
            col=alpha('red',0.2),boxwex=0.5,boxcol=alpha('red',0.5),
            outcol=alpha('red',0.5),outcex=1.15, outpch=21,medcol=alpha('red',0.5),
            outbg=alpha('red',0.5),whiskcol=alpha('red',0.5),whisklty=1,staplecol=alpha('red',0.5),
            frame=F,main='',
            xlab='',ylab='annual stdev [Total Precipitation]',cex.axis=2,cex.main=2.5,cex.lab=1.65, font.main=1)
    points(sd(annual.prcp),col="red",pch=16,cex=3)
    points(sd(annual.prcp.sim),col="black",pch=16,cex=3)
    legend('topleft',c(paste0(label1,' (boot)'),paste0(label2)),
           horiz = F,
           # cex=1.75,
           cex=2.25,
           # title=paste0('at-basin'),
           # title.col = 'purple',title.cex = 2.5,
           fill=c('red','black'),bty = 'n')
    
    
    # (c...)
    #--------------------------------------- #
    ## -- whiskers of long trace metrics -- ##
    #--------------------------------------- #
    lst.annual.obs <- lst.annual.sim <- list()
    lst.annual.obs[[1]] <- aggregate(prcp.site.data,FUN=mean,by=list(wateryears.obs))[,-1]
    lst.annual.sim[[1]] <- aggregate(prcp.site.data.sim,FUN=mean,by=list(wateryears.sim))[,-1]
    lst.annual.obs[[2]] <- aggregate(prcp.site.data,FUN=sd,by=list(wateryears.obs))[,-1]
    lst.annual.sim[[2]] <- aggregate(prcp.site.data.sim,FUN=sd,by=list(wateryears.sim))[,-1]
    
    labs.metrics <- c('mean','stdev')
    
    for (k in 1:length(labs.metrics)){
      annual.obs <- lst.annual.obs[[k]]
      annual.sim <- lst.annual.sim[[k]]
      my.sim.y <- annual.sim
      my.obs.y <- annual.obs
      
      my1.x0 <- apply(my.sim.y,2,FUN=median)
      my1.y0 <- apply(my.obs.y,2,FUN=quantile,0.05)
      my1.y1 <- apply(my.obs.y,2,FUN=quantile,0.95)
      
      my2.x0 <- my2.x1 <- apply(my.obs.y,2,FUN=median)
      my2.y0 <- apply(my.sim.y,2,FUN=quantile,0.05)
      my2.y1 <- apply(my.sim.y,2,FUN=quantile,0.95)
      x_y_lim <- c(min(c(my1.y0,my2.y0)),max(c(my1.y1,my2.y1)))
      
      plot(my1.x0, my2.x0, xlim=x_y_lim, ylim=x_y_lim,
           pch=1, cex=2.25, frame=T,col=NULL,
           xlab='',ylab='',
           cex.axis=2,
           main=paste0(labs.metrics[k]),sub='|5th, median, 95th|',col.sub=alpha('blue',0.5),
           cex.lab=2, cex.main=3, font.main=1,cex.sub=1.75)
      mtext(paste0(label2),1,col='gray60',line=2.25, cex=1.75)
      mtext(label1,2,col='red',line=2.25, cex=1.75)
      arrows(x0=my2.y0,
             y0=my2.x0, x1=my2.y1, y1=my2.x0, code=3, angle=90, length=0.025, col=alpha("gray40",0.75), lwd=0.5)
      arrows(x0=my1.x0,
             y0=my1.y0, x1=my1.x0, y1=my1.y1, code=3, angle=90, length=0.025, col=alpha("red",0.75), lwd=0.5)
      abline(0,1,lwd=1.5,lty=2)
      points(my1.x0, my2.x0, xlim=x_y_lim, ylim=x_y_lim, pch=19, 
             cex=1.5,col=alpha("gray40",0.65))
      legend("topleft",legend=c(NA) ,title=paste0('at-site'),pch=19,cex=2.5,
             col='gray40',title.col = 'purple',inset=0.01,
             horiz = FALSE)
    }
    
    
    dev.off()
  }
  
  # Figure 3:
  # 3) water-year cumulative precipitations, monthly distributions
  {
    # label1 <- 'Obs [1948-2018]'
    # label2 <- 'Sim (WGEN: baseline)'
    
    file.name.fig <- paste0("lines_prcp.cumulative_seasonality_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 15; hh.mom <- 14
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfcol=c(2,1),mai=c(1,1,0.75,0.5))
    
    # (a)
    #--------------------------------------------- #
    ## -- water-year cumulative precipitations -- ##
    #--------------------------------------------- #
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
        single.wy <- prcp.site.data[start_wy[y]:end_wy[y+1],k]
        cumsum.WY.obs.sites[y,k,(1:length(single.wy))] <- cumsum(single.wy)
      }
    }
    psite.obs.mean.metric2 <- apply(cumsum.WY.obs.sites,c(2,3),FUN=mean,na.rm=T)
    psite.obs.metric2 <- apply(cumsum.WY.obs.sites,c(2,3),FUN=median,na.rm=T)
    psite.obs.metric2.10th <- apply(cumsum.WY.obs.sites,c(2,3),FUN=quantile,0.1,na.rm=T)
    psite.obs.metric2.90th <- apply(cumsum.WY.obs.sites,c(2,3),FUN=quantile,0.9,na.rm=T)
    
    mean.psite.obs.metric2 <- as.matrix(apply(psite.obs.metric2, 2, mean))
    mean.psite.obs.metric2.10th <- as.matrix(apply(psite.obs.metric2.10th, 2, mean))
    mean.psite.obs.metric2.90th <- as.matrix(apply(psite.obs.metric2.90th, 2, mean))
    mean.psite.obs.mean.metric2 <- as.matrix(apply(psite.obs.mean.metric2, 2, mean))
    
    start_wy <- which(format(dates.sim,'%m-%d')=='10-01')
    end_wy <- which(format(dates.sim,'%m-%d')=='09-30')
    num.WY.available <- length(start_wy)-1
    cumsum.WY.sim.sites <- array(NA,c(num.WY.available,n.sites,366))
    for(y in 1:num.WY.available){
      extracted.wy.dates <- dates.sim[start_wy[y]:end_wy[y+1]]
      doy.WY.single.wy.idx <- as.numeric(format(extracted.wy.dates,'%j'))
      
      for (k in 1:n.sites){
        
        single.wy <- prcp.site.data.sim[start_wy[y]:end_wy[y+1],k]
        cumsum.WY.sim.sites[y,k,(1:length(single.wy))] <- cumsum(single.wy)
      }
    }
    
    psite.sim.mean.metric2 <- apply(cumsum.WY.sim.sites,c(2,3),FUN=mean,na.rm=T)
    psite.sim.metric2 <- apply(cumsum.WY.sim.sites,c(2,3),FUN=median,na.rm=T)
    psite.sim.metric2.10th <- apply(cumsum.WY.sim.sites,c(2,3),FUN=quantile,0.1,na.rm=T)
    psite.sim.metric2.90th <- apply(cumsum.WY.sim.sites,c(2,3),FUN=quantile,0.9,na.rm=T)
    
    mean.psite.sim.metric2 <- as.matrix(apply(psite.sim.metric2, 2, mean))
    mean.psite.sim.metric2.10th <- as.matrix(apply(psite.sim.metric2.10th, 2, mean))
    mean.psite.sim.metric2.90th <- as.matrix(apply(psite.sim.metric2.90th, 2, mean))
    mean.psite.sim.mean.metric2 <- as.matrix(apply(psite.sim.mean.metric2, 2, mean))
    
    start_months <- ceiling_date(sample.WY.dates %m-% months(1), 'month')
    idx.fst.d.month <- which(sample.WY.dates %in% start_months)
    my.seq.days <- 1:365
    
    # median/mean
    plot(my.seq.days,mean.psite.sim.metric2[my.seq.days],
         type='l',col='black',cex=1.5,ylab='Precipitation [cumulative]',frame=F,
         font.main = 1,sub='',xaxt = "n",xlab='Water Year',lwd=2,
         ylim=c(min(c(mean.psite.sim.metric2.10th[my.seq.days],
                      mean.psite.obs.metric2.10th[my.seq.days])),
                max(c(mean.psite.sim.metric2.90th[my.seq.days],
                      mean.psite.obs.metric2.90th[my.seq.days]))),
         cex.lab=2,cex.axis=1.5,cex.main=3.5,
         )
    
    axis(1,at = my.seq.days[idx.fst.d.month],col='gray50',cex.lab=3,cex.axis=2,
         labels = sample.WY.Months[idx.fst.d.month])
    
    abline(v=my.seq.days[idx.fst.d.month],col='gray90',lty=2)
    
    lines(mean.psite.obs.metric2[my.seq.days],
          col="red",cex=1,lwd=3.5)
    
    lines(mean.psite.sim.metric2.10th[my.seq.days],
          col="gray25",cex=1,lwd=3,lty=2)
    lines(mean.psite.sim.metric2.90th[my.seq.days],
          col="gray25",cex=1,lwd=3,lty=2)
    
    lines(mean.psite.obs.metric2.10th[my.seq.days],
          col=rgb(1,0,0,0.75),cex=2,lwd=3,lty=2)
    lines(mean.psite.obs.metric2.90th[my.seq.days],
          col=rgb(1,0,0,0.75),cex=2,lwd=3,lty=2)
    
    lines(mean.psite.obs.mean.metric2[my.seq.days],
          col="tan1",cex=2,lwd=8, lty=1)
    lines(mean.psite.sim.mean.metric2[my.seq.days],
          col="blue",cex=2,lwd=3.5, lty=1)
    
    legend("topleft", legend=c(paste0(label2," (10,median,90th)"),
                               paste0(label1," (10,median,90th)"),
                               paste0(label2," (mean)"),
                               paste0(label1," (mean)")),
           col=c("black", "red","blue", "tan1"),
           lty=c(1,1,1,1),lwd=c(2,2,2,4),
           # cex=1.5,
           cex=2,
           horiz = FALSE,bty='n',
           # title=paste0('at-basin'),
           # title.col = 'purple',title.cex = 2
    )
    
    # (b)
    #----------------------------------------- #
    ## -- water-year monthly distributions -- ##
    #----------------------------------------- #
    years.months.weather <- format(dates.weather,'%Y-%m')
    unique.years.months.weather <- format(ym(unique(years.months.weather)),'%m')
    
    dummy.long.years.vec <- long.years.vec+2000
    years.months.sim <- paste0(dummy.long.years.vec,'-',long.months.vec)
    unique.years.months.sim <- format(ym(unique(years.months.sim)),'%m')
    
    psite.obs.metric <- aggregate(mean.prcp.site,FUN=sum,by=list(years.months.weather))[,2]
    psite.obs.metric2 <- aggregate(psite.obs.metric,FUN=median,by=list(unique.years.months.weather))[,2]
    psite.obs.metric.q1th <- aggregate(psite.obs.metric,FUN=quantile,0.05,by=list(unique.years.months.weather))[,2]
    psite.obs.metric.q2th <- aggregate(psite.obs.metric,FUN=quantile,0.95,by=list(unique.years.months.weather))[,2]
    
    psite.sim.metric <- aggregate(mean.prcp.site.sim,FUN=sum,by=list(years.months.sim))[,-1]
    psite.sim.metric2 <- aggregate(psite.sim.metric,FUN=median,by=list(unique.years.months.sim))[,-1]
    psite.sim.metric.q1th <- aggregate(psite.sim.metric,FUN=quantile,0.05,by=list(unique.years.months.sim))[,-1]
    psite.sim.metric.q2th <- aggregate(psite.sim.metric,FUN=quantile,0.95,by=list(unique.years.months.sim))[,-1]
    
    ##reordering.dummy.aggre <- c(1,10,11,12,2,3,4,5,6,7,8,9) #stupidty of aggregate for sim
    reordering.dummy.aggre <- c(1,5,6,7,8,9,10,11,12,2,3,4)
    psite.sim.metric2 <- psite.sim.metric2[reordering.dummy.aggre]
    psite.sim.metric.q1th <- psite.sim.metric.q1th[reordering.dummy.aggre]
    psite.sim.metric.q2th <- psite.sim.metric.q2th[reordering.dummy.aggre]
    
    plot(months,psite.sim.metric2,
         pch=4,col='black',cex=4.5,ylab='Precipitation [total]',frame=F,
         font.main = 1,sub='',xlab='calendar month',
         cex.lab=2,cex.axis=1.75,cex.main=2,
         main='',
         ylim=c(min(psite.sim.metric.q1th),max(psite.sim.metric.q2th)))
    arrows(months,psite.obs.metric.q1th
           ,months,psite.obs.metric.q2th,
           length = 0,lwd=4,col=alpha('red',0.65))
    
    arrows(months,psite.sim.metric.q1th
           ,months,psite.sim.metric.q2th,
           length = 0,lwd=28,col=alpha('gray',0.5))
    points(psite.obs.metric2,col="red",cex=3.5,pch=19)
    legend("top", legend=c(paste0(label2," (5,median,95th)"),
                           paste0(label1," (5,median,95th)")),
           col=c("black", "red"), lty=c(1,0),cex=2,
           pch=c(4,19),
           horiz = F,bty='n',
           # title=paste0('at-basin'),
           # title.col = 'purple',title.cex = 2
    )
    
    dev.off()
  }
  
  
  # Figure 4:
  # 4) PBIAS for a list of metrics for precipitation as boxplots
  {
    # label1 <- 'Obs [1948-2018]'
    # label2 <- 'Sim (WGEN: baseline)'
    
    file.name.fig <- paste0("boxplots_prcp.PBIAS_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 16; hh.mom <- 12
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfcol=c(1,1),mai=c(1,1,0.75,0.5),font.main=1)
    
    ##Metrics for PBIAS
    nsites <- n.sites
    metrics <- cbind(
      'SD'=sapply(1:nsites,function(x){return(100*(sd(wy.data.sim[,x])-sd(wy.data.obs[,x]))/sd(wy.data.obs[,x]))}),
      'MEAN'=sapply(1:nsites,function(x){return(100*(mean(wy.data.sim[,x])-mean(wy.data.obs[,x]))/mean(wy.data.obs[,x]))}),
      'MEDIAN'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.5)-quantile(wy.data.obs[,x],0.5))/quantile(wy.data.obs[,x],0.5))}),
      'MIN'=sapply(1:nsites,function(x){return(100*(min(wy.data.sim[,x])-min(wy.data.obs[,x]))/min(wy.data.obs[,x]))}),
      'MIN2'=sapply(1:nsites,function(x){return(100*(min(rollmean(wy.data.sim[,x],2))-min(rollmean(wy.data.obs[,x],2)))/min(rollmean(wy.data.obs[,x],2)))}),
      'MIN3'=sapply(1:nsites,function(x){return(100*(min(rollmean(wy.data.sim[,x],3))-min(rollmean(wy.data.obs[,x],3)))/min(rollmean(wy.data.obs[,x],3)))}),
      'MIN5'=sapply(1:nsites,function(x){return(100*(min(rollmean(wy.data.sim[,x],5))-min(rollmean(wy.data.obs[,x],5)))/min(rollmean(wy.data.obs[,x],5)))}),
      'Q05'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.05)-quantile(wy.data.obs[,x],0.05))/quantile(wy.data.obs[,x],0.05))}),
      'Q90'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.9)-quantile(wy.data.obs[,x],0.9))/quantile(wy.data.obs[,x],0.9))}),
      'Q99'=sapply(1:nsites,function(x){return(100*(quantile(wy.data.sim[,x],0.99)-quantile(wy.data.obs[,x],0.99))/quantile(wy.data.obs[,x],0.99))}),
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
    
    spell.stats.metrics <- t(100*(fct.spell.stats(prcp.site.data.sim) - fct.spell.stats(prcp.site.data))/fct.spell.stats(prcp.site.data))
    pBias.full.set.metrics <- cbind(metrics,spell.stats.metrics)
    colnames(pBias.full.set.metrics) <- c(colnames(metrics),c('MEAN.W.SPELL','MAX.W.SPELL','MEAN.D.SPELL','MAX.D.SPELL'))
    rownames(pBias.full.set.metrics) <- NULL
    
    my.median.values <- paste0(round(apply(pBias.full.set.metrics,2,median),0),"%")
    
    boxplot(pBias.full.set.metrics,boxwex=.65,
            col='wheat',outline=F,frame=F,cex.axis=2.5,cex.lab=2.5,cex.main=3.5,font.main=1,
            xlab='Precipitation [water year]',ylab='PBIAS[%]',xaxt = "n",ylim=c(-60,60),
            border=alpha('black',0.6),
            whisklty = 1,medcol = alpha('red',0.5),medlwd=1)
    abline(v=1:dim(pBias.full.set.metrics)[2],col=alpha('gray',0.45),lty=3)
    abline(h=0,col=alpha('thistle',0.65),lwd=1,lty=1)
    axis(side = 1, labels = FALSE)
    text(x = 1:dim(pBias.full.set.metrics)[2],
         y = par("usr")[3],labels = colnames(pBias.full.set.metrics),
         xpd = NA,adj=c(rep(0,dim(metrics)[2]),rep(0.1,dim(spell.stats.metrics)[2])),
         ## Rotate the labels by 35 degrees.
         srt = 90,cex = 2, col='gray50')
    text(x = 1:dim(pBias.full.set.metrics)[2],
         bg ="lightblue",
         y = 60,labels = my.median.values,
         xpd = NA,adj=0,
         ## Rotate the labels by 35 degrees.
         srt = 45,cex = 1.8, col=alpha('red',0.5))
    text(x = 0,
         y = 60,labels = 'median',
         xpd = NA,adj=0,
         ## Rotate the labels by 35 degrees.
         srt = 90,cex = 2, col=alpha('red',0.5))
    mtext(paste0(label2,' vs.', label1),4,cex=2,col='black')
    legend('top',legend=c(paste0('sites: 1-',nsites)),
           title='at-site',box.col = "white",
           title.col = 'purple',title.cex = 2,
           fill='wheat',cex=1.6,inset = 0.08)
    
    dev.off()
  }
  
  
  # Figure 5:
  # 5) temperature cdf, mean, standard variation (yearly, year-to-year)
  {

    file.name.fig <- paste0("lines_tmean_mean_cdf_sd_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 14; hh.mom <- 14
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfcol=c(2,2),mai=c(1,1,1,0.5))
    
    
    # (a)
    #------------------------------------------ #
    ## -- cdf of temperature distribution -- ##
    #------------------------------------------ #
    annual.tmean <- stats::aggregate(mean.tmean.site,FUN=mean,by=list(wateryears.obs))[,2]
    annual.tmean <- annual.tmean[2:(length(annual.tmean)-1)]   #drop first and last to only keep full water years
    annual.tmean.boot.sort <- array(NA,c(length(annual.tmean),1000))
    annual.tmean.sd <- array(NA,1000)
    for (i in 1:1000) {
      annual.tmean.boot.sort[,i] <- base::sort(sample(annual.tmean,size=length(annual.tmean),replace=TRUE))
      annual.tmean.sd[i] <- stats::sd(sample(annual.tmean,size=length(annual.tmean),replace=TRUE))
    }
    boot.95 <- apply(annual.tmean.boot.sort,1,FUN=quantile,.95)
    boot.05 <- apply(annual.tmean.boot.sort,1,FUN=quantile,.05)
    #get annual average temp from simulation
    tmean.basin.sim <- mean.tmean.site.sim
    annual.tmean.sim <- stats::aggregate(tmean.basin.sim,FUN=mean,by=list(wateryears.sim),na.rm=T)[,2]
    annual.tmean.sim <- annual.tmean.sim[2:(length(annual.tmean.sim)-1)] #drop first and last to only keep full water years
    
    
    y.hat <- annual.tmean
    EP.hat <- (1:length(y.hat))/(length(y.hat)+1)
    y <- annual.tmean.sim
    EP <- (1:length(y))/(length(y)+1)
    plot(EP.hat,sort(y.hat),log="y",type='l',frame=F,
         col="red",lwd=5,
         xlab='Exceedance Probability',ylab='Temperature [WY]',
         cex.axis=2,cex.main=4,cex.lab=2, font.main=1)
    lines(EP,sort(y),lty=1, col='gray50',lwd=6)
    polygon(x=c(EP.hat,rev(EP.hat)),y=c(boot.05,rev(boot.95)),
            col=alpha('red',0.15),border=NA)
    legend('topleft',c(paste0(label1,' (boot)'),paste0(label2)),
           horiz = F,cex=1.75,title=paste0('at-basin'),
           title.col = 'purple',title.cex = 2.5,
           fill=c('red','black'),bty = 'n')
    
    # (b)
    #----------------------------------------- #
    ## -- year-to-year SD of temperature -- ##
    #----------------------------------------- #
    annual.tmean <- aggregate(mean.tmean.site,FUN=mean,by=list(wateryears.obs))[,2]
    annual.tmean <- annual.tmean[2:(length(annual.tmean)-1)]   #drop first and last to only keep full water years
    annual.tmean.sd <- array(NA,1000)
    for (i in 1:1000) {
      annual.tmean.sd[i] <- sd(sample(annual.tmean,size=length(annual.tmean),replace=TRUE))
    }
    #get annual average tmean from simulation
    tmean.basin.sim <- mean.tmean.site.sim
    annual.tmean.sim <- aggregate(tmean.basin.sim,FUN=mean,by=list(wateryears.sim),na.rm=T)[,2]
    annual.tmean.sim <- annual.tmean.sim[2:(length(annual.tmean.sim)-1)] #drop first and last to only keep full water years
    
    boxplot(annual.tmean.sd,
            col=alpha('red',0.2),boxwex=0.5,boxcol=alpha('red',0.5),
            outcol=alpha('red',0.5),outcex=1.15, outpch=21,medcol=alpha('red',0.5),
            outbg=alpha('red',0.5),whiskcol=alpha('red',0.5),whisklty=1,staplecol=alpha('red',0.5),
            frame=F,main='',
            ylim = c(min(min(annual.tmean.sd),sd(annual.tmean),sd(annual.tmean.sim)),
                     max(max(annual.tmean.sd),sd(annual.tmean),sd(annual.tmean.sim))),
            xlab='',ylab='annual stdev [Average Temperature]',
            cex.axis=2,cex.main=2.5,cex.lab=1.65, font.main=1)
    points(sd(annual.tmean),col="red",pch=16,cex=3)
    points(sd(annual.tmean.sim),col="black",pch=16,cex=3)
    legend('topleft',c(paste0(label1,' (boot)'),paste0(label2)),
           horiz = F,cex=2,title=paste0('at-basin'),
           title.col = 'purple',title.cex = 2.5,
           col=c('red','black'),pch=16,bty = 'n')
    
    
    # (c...)
    #--------------------------------------- #
    ## -- whiskers of long trace metrics -- ##
    #--------------------------------------- #
    lst.annual.obs <- lst.annual.sim <- list()
    lst.annual.obs[[1]] <- aggregate(tmean.site.data,FUN=mean,by=list(wateryears.obs))[,-1]
    lst.annual.sim[[1]] <- aggregate(tmean.site.data.sim,FUN=mean,by=list(wateryears.sim))[,-1]
    lst.annual.obs[[2]] <- aggregate(tmean.site.data,FUN=sd,by=list(wateryears.obs))[,-1]
    lst.annual.sim[[2]] <- aggregate(tmean.site.data.sim,FUN=sd,by=list(wateryears.sim))[,-1]
    
    labs.metrics <- c('mean','stdev')
    
    for (k in 1:length(labs.metrics)){
      annual.obs <- lst.annual.obs[[k]]
      annual.sim <- lst.annual.sim[[k]]
      my.sim.y <- annual.sim
      my.obs.y <- annual.obs
      
      my1.x0 <- apply(my.sim.y,2,FUN=median)
      my1.y0 <- apply(my.obs.y,2,FUN=quantile,0.05)
      my1.y1 <- apply(my.obs.y,2,FUN=quantile,0.95)
      
      my2.x0 <- my2.x1 <- apply(my.obs.y,2,FUN=median)
      my2.y0 <- apply(my.sim.y,2,FUN=quantile,0.05)
      my2.y1 <- apply(my.sim.y,2,FUN=quantile,0.95)
      x_y_lim <- c(min(c(my1.y0,my2.y0)),max(c(my1.y1,my2.y1)))
      
      plot(my1.x0, my2.x0, xlim=x_y_lim, ylim=x_y_lim,
           pch=1, cex=2.25, frame=T,col=NULL,
           xlab='',ylab='',
           cex.axis=2,
           main=paste0(labs.metrics[k]),sub='|5th, median, 95th|',col.sub=alpha('blue',0.5),
           cex.lab=2, cex.main=3, font.main=1,cex.sub=1.75)
      mtext(paste0(label2),1,col='gray60',line=2.5, cex=1.75)
      mtext(label1,2,col='red',line=2.5, cex=1.75)
      arrows(x0=my2.y0,
             y0=my2.x0, x1=my2.y1, y1=my2.x0, code=3, angle=90, length=0.025, col=alpha("gray40",0.75), lwd=0.5)
      arrows(x0=my1.x0,
             y0=my1.y0, x1=my1.x0, y1=my1.y1, code=3, angle=90, length=0.025, col=alpha("red",0.75), lwd=0.5)
      abline(0,1,lwd=1.5,lty=2)
      points(my1.x0, my2.x0, xlim=x_y_lim, ylim=x_y_lim, pch=19, 
             cex=1.5,col=alpha("gray40",0.65))
      legend("topleft",legend=c(NA) ,title=paste0('at-site'),pch=19,cex=2.5,
             col='gray40',bty = 'n',title.col = 'purple',
             horiz = FALSE,inset=0.0)
    }
    
    dev.off()
  }
  
  
  # Figure 6:
  # 6) heat and cold waves/stress temperature
  {
    temp.specific.diagnostics <- function(temp.mean,temp.max,temp.min,
                                          my.range=NULL,mc.cores=1)
    {
      
      # temp.min <- temp.max <- temp.mean
      #This function will calculate a series of specific diagnostic statistics for temp, which is a list of independent temp vectors 
      #
      #Arguments:
      #temp = a list of temperature (min,max) time series for each site (i.e., a list of temperature matrices), with list dimension d
      
      # Heatwave definition: Three or more consecutive days over 90 degrees F. (based on the state of MA definition)
      # Coldwave definition: Ten or more consecutive days below 20 degrees F. (User-specific criteria)
      hw.temp <- 32.22 # 32.22 deg. C  == 90 deg. F
      cw.temp <- -7 # -7 deg. C  == 20 deg. F
      min.days.hw <- 3
      min.days.cw <- 5 # this was 10 for MA project
      
      temp.extreme.high.thresh <- 30    # trace extreme high temperature == 86 F
      temp.extreme.low.thresh <- 0    # trace extreme low temperature == 32 F
      ma.options.rolling <- c(3) # for 3-day moving average
      
      rep.row<-function(x,n){
        matrix(rep(x,each=n),nrow=n)
      }
      
      d <- length(temp.min)
      n.sites <- dim(temp.min[[1]])[2]
      n.days <- dim(temp.min[[1]])[1]
      stat.names <- c('number.heatwaves','mean.heatwaves.duration','max.heatwaves.duration',
                      'number.coldwaves','mean.coldwaves.duration','max.coldwaves.duration',
                      'number.heatstress','number.coldstress',
                      'mean.intensity.heatwaves','mean.intensity.coldwaves',
                      'max.intensity.heatwaves','max.intensity.coldwaves')
      n.stats <- length(stat.names)
      if(is.null(my.range)) {my.range <- 1:dim(temp.min[[1]])[1]}
      
      #calculate at-site statistics  
      stat.list <- mclapply(X=1:d,mc.preschedule=TRUE,mc.cores=mc.cores,FUN=function(i){  
        
        my.stats <- array(0,c(n.sites,n.stats))
        
        
        if (sum(is.infinite(temp.min[[i]]))>0){
          idx.dummy <- which(is.infinite(temp.min[[i]]),arr.ind = T)
          temp.min[[i]][idx.dummy[,1],idx.dummy[,2]] <- NA} # checkpoints to not have "Inf" before going ahead
        if (sum(is.infinite(temp.max[[i]]))>0){
          idx.dummy <- which(is.infinite(temp.max[[i]]),arr.ind = T)
          temp.max[[i]][idx.dummy[,1],idx.dummy[,2]] <- NA} # checkpoints to not have "Inf" before going ahead
        
        # cur.temp.min <- temp.min[[i]][my.range,]
        # cur.temp.max <- temp.max[[i]][my.range,]
        # cur.temp <- (cur.temp.min+cur.temp.max)/2
        
        # Number of heatwaves, mean length, max values
        cur.temp <- temp.max[[i]][my.range,]
        
        bin.hw <- (cur.temp>hw.temp)*1
        bin.hw.cur.temp <- cur.temp*bin.hw
        if (sum(bin.hw>0)>0){
          bin.cur.temp <- bin.hw.cur.temp
          mat.my.stats <- sapply(1:n.sites,function(x){
            my.intensity <- array(0,c(5,1))
            obs_time <- data.frame(cbind(bin.cur.temp[,x],my.range))
            colnames(obs_time) <- c('obs','time')
            if(!sum(obs_time$obs)==0){
              
              clust.obs_time_durations <- try(clust(obs_time,u=0,
                                                    tim.cond=0,clust.max=FALSE,
                                                    plot=FALSE))
              site.declust <- c()
              site.declust$intensity <- sapply(1:length(clust.obs_time_durations),function(y){
                return(mean(clust.obs_time_durations[[y]][2,])-hw.temp)})
              site.declust$maxintensity <- sapply(1:length(clust.obs_time_durations),function(y){
                return(max(clust.obs_time_durations[[y]][2,])-hw.temp)})
              site.declust$durations <- sapply(1:length(clust.obs_time_durations),function(y){
                return(length(clust.obs_time_durations[[y]][2,]))})
              if (sum(site.declust$durations>=min.days.hw & site.declust$intensity>=0)>0){
                
                my.intensity[1,1] <- mean(site.declust$intensity[site.declust$durations>=min.days.hw & site.declust$intensity>=0])
                my.intensity[2,1] <- mean(site.declust$maxintensity[site.declust$durations>=min.days.hw & site.declust$intensity>=0])
                my.intensity[3,1] <- mean(site.declust$durations[site.declust$durations>=min.days.hw & site.declust$intensity>=0])
                my.intensity[4,1] <- max(site.declust$durations[site.declust$durations>=min.days.hw & site.declust$intensity>=0])
                my.intensity[5,1] <- sum(site.declust$durations>=min.days.hw & site.declust$intensity>=0)
                
              }
            }
            return(my.intensity)
          })
          my.stats[,1] <- mat.my.stats[5,] # freq
          my.stats[,2] <- mat.my.stats[3,] # mean dur
          my.stats[,3] <- mat.my.stats[4,] # max dur
          my.stats[,9] <- mat.my.stats[1,] # mean intensity
          my.stats[,11] <- mat.my.stats[2,] # max intensity
        }
        
        # Number of coldwaves, mean length, max values
        cur.temp <- temp.min[[i]][my.range,]
        
        bin.cw <- (cur.temp<cw.temp)*1
        bin.cw.cur.temp <- cur.temp*bin.cw
        if (sum(bin.cw>0)>0){
          bin.cur.temp <- bin.cw.cur.temp
          mat.my.stats <- sapply(1:n.sites,function(x){
            my.intensity <- array(0,c(5,1))
            obs_time <- data.frame(cbind(-bin.cur.temp[,x],my.range)) # x -1 for POT
            colnames(obs_time) <- c('obs','time')
            if(!sum(obs_time$obs)==0){
              
              clust.obs_time_durations <- try(clust(obs_time,u=0,
                                                    tim.cond=0,clust.max=FALSE,
                                                    plot=FALSE))
              site.declust <- c()
              site.declust$intensity <- sapply(1:length(clust.obs_time_durations),function(y){
                return(mean(clust.obs_time_durations[[y]][2,])-abs(cw.temp))})
              site.declust$maxintensity <- sapply(1:length(clust.obs_time_durations),function(y){
                return(max(clust.obs_time_durations[[y]][2,])-abs(cw.temp))})
              site.declust$durations <- sapply(1:length(clust.obs_time_durations),function(y){
                return(length(clust.obs_time_durations[[y]][2,]))})
              if (sum(site.declust$durations>=min.days.cw & site.declust$intensity>=0)>0){
                my.intensity[1,1] <- -mean(site.declust$intensity[site.declust$durations>=min.days.cw & site.declust$intensity>=0]) # x -1 for POT
                my.intensity[2,1] <- -mean(site.declust$maxintensity[site.declust$durations>=min.days.cw & site.declust$intensity>=0]) # x -1 for POT
                my.intensity[3,1] <- mean(site.declust$durations[site.declust$durations>=min.days.cw & site.declust$intensity>=0])
                my.intensity[4,1] <- max(site.declust$durations[site.declust$durations>=min.days.cw & site.declust$intensity>=0])
                my.intensity[5,1] <- sum(site.declust$durations>=min.days.cw & site.declust$intensity>=0)
                
              }
            }
            return(my.intensity)
          })
          my.stats[,4] <- mat.my.stats[5,] # freq
          my.stats[,5] <- mat.my.stats[3,] # mean dur
          my.stats[,6] <- mat.my.stats[4,] # max dur
          my.stats[,10] <- mat.my.stats[1,] # mean intensity 
          my.stats[,12] <- mat.my.stats[2,] # max intensity
        }
        
        
        # cur.temp.min <- temp.min[[i]][my.range,]
        # cur.temp.max <- temp.max[[i]][my.range,]
        # cur.temp <- (cur.temp.min+cur.temp.max)/2
        
        #heatstress
        cur.temp <- temp.max[[i]][my.range,]
        #extreme.high.temp.val <- apply(cur.temp,2,FUN=quantile,temp.extreme.high.thresh,na.rm=T)
        mat.extreme.high.temp.val <- rep.row(temp.extreme.high.thresh,n.days-ma.options.rolling[1]+1)
        y1 <- apply(cur.temp,2,function(x){rollapply(data=x, FUN = mean, width = ma.options.rolling[1])})
        
        #coldstress
        cur.temp <- temp.min[[i]][my.range,]
        #extreme.low.temp.val <- apply(cur.temp,2,FUN=quantile,temp.extreme.low.thresh,na.rm=T)
        mat.extreme.low.temp.val <- rep.row(temp.extreme.low.thresh,n.days-ma.options.rolling[1]+1)
        y2 <- apply(cur.temp,2,function(x){rollapply(data=x, FUN = mean, width = ma.options.rolling[1])})
        
        #my.stats[,7] <- max(c(0,sum(y>mat.extreme.high.temp.val)))
        my.stats[,7] <- apply(y1,2,function(x){sum(x>mat.extreme.high.temp.val)})
        
        #my.stats[,8] <- max(c(0,sum(y<mat.extreme.low.temp.val)))
        my.stats[,8] <- apply(y2,2,function(x){sum(x<mat.extreme.low.temp.val)})
        
        return(my.stats)
      })
      my.stats <- unlist(stat.list)
      dim(my.stats) <- c(n.sites,n.stats,d)    
      dimnames(my.stats)[[2]] <- stat.names    
      
      return(my.stats)
      
    }
    
    temp.specific.obs <- temp.specific.diagnostics(list(tmean.site.data),
                                                   list(tmax.site.data),
                                                   list(tmin.site.data))
    temp.specific.obs <- data.frame(temp.specific.obs[,,1])
    # converting to per year metrics:
    temp.specific.obs[,c(1,4,7,8)] <- temp.specific.obs[,c(1,4,7,8)]/length(unique(years.weather))
    temp.specific.sim <- temp.specific.diagnostics(list(tmean.site.data.sim),
                                                   list(tmax.site.data.sim),
                                                   list(tmin.site.data.sim))
    temp.specific.sim <- data.frame(temp.specific.sim[,,1])
    # converting to per year metrics:
    temp.specific.sim[,c(1,4,7,8)] <- temp.specific.sim[,c(1,4,7,8)]/total.num.yrs
    
    temp.specific.obs[is.infinite(as.matrix(temp.specific.obs))] <- NA
    temp.specific.sim[is.infinite(as.matrix(temp.specific.sim))] <- NA
    
    # label1 <- 'Obs [1948-2018]'
    # label2 <- 'Sim (WGEN: baseline)'
    labs.metrics <- c('Heat Wave Frequency','Heat Wave Duration (mean)','Heat Wave Duration (max)',
                      'Cold Wave Frequency','Cold Wave Duration (mean)','Cold Wave Duration (max)',
                      'Heat Stress Frequency','Cold Stress Frequency',
                      'Heat Wave Intensity (mean)','Cold Wave Intensity (mean)',
                      'Heat Wave Intensity (max)','Cold Wave Intensity (max)')
    unit.metric <- c('#','(days)','(days)',
                     '#','(days)','(days)',
                     '#','#',
                     '(above-threshold, C)','(below-threshold, C)',
                     '(above-threshold, C)','(below-threshold, C)')
    
    file.name.fig <- paste0("scatters_tmax_tmin_heat.cold.specific_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 18; hh.mom <- 10
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfrow=c(2,4),mar=c(5.1, 6.1, 4.1, 2.1))
    
    for (k in c(2,3,9,11,5,6,10,12)){
      my1.x0 <- temp.specific.obs[,k]
      my2.x0 <- temp.specific.sim[,k]
      x_y_lim <- c(min(c(my1.x0,my2.x0),na.rm=T),
                   max(c(my1.x0,my2.x0),na.rm=T))
      
      plot(my2.x0,my1.x0, xlim=x_y_lim, ylim=x_y_lim,
           pch=1, cex=2, frame=F,col=NULL,
           xlab='',ylab='',
           cex.axis=1.5,
           main=paste0(labs.metrics[k]),sub=paste0(unit.metric[k]),col.sub=alpha('blue',0.5),
           cex.lab=1.25, cex.main=2, font.main=1,cex.sub=1.5)
      mtext(paste0(label2),1,col='gray10',line=2.35)
      mtext(paste0(label1),2,col='red',line=2.35)
      
      abline(a=0,b=1, col='gray75',lwd=1.5, lty=2)
      points(my2.x0,my1.x0, xlim=x_y_lim, ylim=x_y_lim, pch=16, 
             cex=2.5,col=alpha("gray40",0.65))
      legend("topleft",pch=16,
             col='gray30',
             horiz = FALSE,inset=0.01,text.col='purple',
             cex=1.5,legend=paste0('at-site'),
             title.col = 'black',title.cex = 2.5,bty = 'n')
    }
    dev.off()
  }
  
  
  # Figure 7:
  # 7) precipitation minima, dry-spells, drought stats
  {

    file.name.fig <- paste0("hist_prcp.n.year.droughts_minima_dry.spells_",n.sites,".files_s",
                            num.states,"_with_",num.iter,"_ens.png")
    fname.fig <- paste0("./Figures/",file.name.fig)
    ww.mom <- 26; hh.mom <- 14
    png(fname.fig,width=ww.mom,height=hh.mom,units="in",res=300)
    par(mfcol=c(2,4),mai=c(0.75,1,1,0.5))
    
    # (a)
    #----------------------------------------- #
    ## -- worst/long-term droughts -- ##
    #----------------------------------------- #
    prcp.min.extreme.diagnostics <- function(prcp,my.range=NULL,mc.cores=1) {
      nyears <- 10
      
      d <- length(prcp)
      n.sites <- dim(prcp[[1]])[2]
      # ma.options.rolling <- seq(1,nyears)*365 # n-day moving averages
      ma.options.rolling <- seq(1,nyears) # n-year moving averages
      stat.names2 <- paste0(seq(1,nyears),'.year.min')
      #calculate at-site statistics: min
      n.stats <- length(stat.names2)
      
      
      if(is.null(my.range)) {my.range <- 1:dim(prcp[[1]])[1]}
      stat.list2 <- mclapply(X=1:d,mc.preschedule=TRUE,mc.cores=mc.cores,FUN=function(i){  
        my.stats2 <- array(NA,c(n.sites,n.stats))
        cur.prcp <- as.matrix(prcp[[i]][my.range,])
        if (sum(is.infinite(cur.prcp))>0){
          idx.dummy <- which(is.infinite(cur.prcp),arr.ind = T)
          #cur.prcp[idx.dummy[,1],idx.dummy[,2]] <-.Machine$double.eps # checkpoints to not have "Inf" before going ahead
          cur.prcp[idx.dummy[,1],idx.dummy[,2]] <- NA
        }
        
        for (s in 1:n.stats) {
          my.stats2[,s] <- apply(cur.prcp,2,function(x){min(rollapply(data=x, FUN = mean, width = ma.options.rolling[s],na.rm=T))})
          print(paste('stats:',stat.names2[s]))
        }
        colnames(my.stats2) <- stat.names2
        return(my.stats2)
      })
      my.stats2 <- unlist(stat.list2)
      dim(my.stats2) <- c(n.sites,n.stats,d)    
      dimnames(my.stats2)[[2]] <- stat.names2
      return(my.stats2)
    }
    
    n.dig <- 5
    #individual site diagnostics
    nyears <- 10
    stat.names2 <- paste0(seq(1,nyears),'.year')
    stat.names22 <- paste0(seq(1,nyears),'-yr')
    obs.diagnostics.min <- prcp.min.extreme.diagnostics(list(mean.wy.prcp.site))
    obs.diagnostics.min.extreme <- as.matrix(obs.diagnostics.min[,,1])
    obs.stats.min <- round(obs.diagnostics.min.extreme,digits = n.dig)
    
    sim.diagnostics.min <- prcp.min.extreme.diagnostics(list(mean.wy.prcp.site.sim))
    sim.diagnostics.min.extreme <- as.matrix(sim.diagnostics.min[,,1])
    
    n.yr.obs.diagnostics.min.extreme.all <- obs.diagnostics.min.extreme
    n.yr.sim.diagnostics.min.all <- sim.diagnostics.min.extreme
    
    min.obs.sim <- min(c(n.yr.obs.diagnostics.min.extreme.all,n.yr.sim.diagnostics.min.all))
    max.obs.sim <- max(c(n.yr.obs.diagnostics.min.extreme.all,n.yr.sim.diagnostics.min.all))
    plot(n.yr.sim.diagnostics.min.all,
         n.yr.obs.diagnostics.min.extreme.all,
         main='Worst Droughts', font.main = 1,
         cex.axis=3,cex.lab=3.5,cex.main=4,
         font.main=1,cex.sub=1.15,
         frame=F,
         xlab='',cex=6,xlim=c(min.obs.sim,
                              max.obs.sim),
         ylim=c(min.obs.sim,
                max.obs.sim),
         ylab='',col=seq(1,nyears),pch=seq(1,nyears))
    abline(h=n.yr.obs.diagnostics.min.extreme.all,
           lty=2,col=alpha(seq(1,nyears),0.25))
    abline(v=n.yr.sim.diagnostics.min.all,
           lty=2,col=alpha(seq(1,nyears),0.25))
    abline(a=0,b=1, col='gray75',lwd=1.5, lty=2)
    mtext(paste0(label2),1,col='gray30',line=3.5, cex=1.75)
    mtext(paste0(label1),2,col='red',line=3.5, cex=1.75)
    legend("bottomright", legend=stat.names22,pch=seq(1,nyears),
           cex=2.5,
           col=seq(1,nyears),title ='WY Total',
           horiz = FALSE,inset=c(0,0.02),ncol=2)
    #mtext(paste0('worst droughts'),4,cex=1.85,col='black')
    mtext(paste0('at-basin'),4,adj=0.925, cex=2,col='purple')
    
    
    # (b)
    #----------------------------------------- #
    ## -- dry spells -- ##
    #----------------------------------------- #
    spells.labels <- c('mean.wet.spell','max.wet.spell',
                       'mean.dry.spell','max.dry.spell')
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
    my.sim.col <- c('black')
    my.obs.col <- c('red')
    obs.spell.stats <- fct.spell.stats(prcp.site.data)
    sim.spell.stats <- fct.spell.stats(prcp.site.data.sim)
    x1 <- obs.spell.stats[3,]
    x2 <- sim.spell.stats[3,]
    x11 <- obs.spell.stats[4,]
    x22 <- sim.spell.stats[4,]
    
    min_max_xs <- base::range(c(x1,x2))
    boxplot(x1,horizontal = T,main='',ylim=min_max_xs,
            col=alpha(my.obs.col,.2),outline=T,whisklty = 1,
            boxcol=my.obs.col,medcol=my.obs.col,
            whiskcol=my.obs.col,staplecol=my.obs.col,outcol=my.obs.col,
            boxwex=.35,xlab='Mean length (days)',
            cex.axis=2.75,cex.lab=2.5,cex.main=4,font.main=1)
    boxplot(x2,add=T,horizontal = T,ylim=min_max_xs,
            col=alpha(my.sim.col,.2),outline=T,whisklty = 1,
            boxcol=my.sim.col,medcol=my.sim.col,
            whiskcol=my.sim.col,staplecol=my.sim.col,outcol=my.sim.col,
            boxwex=.2,
            cex.axis=2.75,cex.lab=2.5,cex.main=5,font.main=1)
    text(median(c(x1,x2)),0.65,cex=2.75,col='blue',
         paste0('max. Obs=',round(max(x11),1),
                ',\nmax. Sim=',round(max(x22),1)))
    mtext('Dry Spells',3,cex=2.5,line = 1)
    legend("topleft", legend=c(label1,label2),
           fill=c("red", alpha(my.sim.col,.5)),cex=2.5,
           bty = "n",horiz = FALSE,
           title='at-site',title.col = 'purple',title.cex = 3.5)
    
    # (c...)
    #------------------------------------------------- #
    ## -- frequency of annual rolling mean drought -- ##
    #------------------------------------------------- #
    drought.rol.mean <- c(1,2,3,5,10,30)
    labs.metrics <- paste0(drought.rol.mean,'-yr')
    annual.prcp <- apply(wy.data.obs,1,mean)
    annual.prcp.sim <- apply(wy.data.sim,1,mean)
    for (my.opt in 1:length(drought.rol.mean)){
      obs.HRU <- rollmean(annual.prcp,drought.rol.mean[my.opt])
      wgen.trace <- rollmean(annual.prcp.sim,drought.rol.mean[my.opt])
      x_min_max <- c(min(obs.HRU,wgen.trace),max(obs.HRU,wgen.trace))
      hgA <- hist(wgen.trace,plot = FALSE,freq = T,
                  main=paste0(labs.metrics[my.opt]),font.main = 1,
                  cex.axis=2,cex.lab=1.5,cex.main=4,
                  xlab='Precipitation [WY]',xlim=x_min_max,
                  ylab='frequency [#WYs]',col='gray60',lty="blank",breaks = 25)
      hgB <- hist(obs.HRU,plot = FALSE,freq = T,
                  main=paste0(labs.metrics[my.opt]),breaks = 25)
      plot(hgA, col = alpha('gray',0.4),lty="blank",
           main=paste0(labs.metrics[my.opt]),font.main = 1,
           cex.axis=2.25,cex.lab=2.5,cex.main=4,
           xlab='Precipitation',xlim=x_min_max,
           ylab='# WYs') # Plot 1st histogram using a transparent color
      plot(hgB, col = alpha('red',0.4), add = TRUE,lty="blank") # Add 2nd histogram using different color
      abline(v=median(wgen.trace),
             col='black',lwd=2,lty=1)
      abline(v=median(obs.HRU),
             col='red',lwd=3,lty=1)
      abline(v=min(wgen.trace),
             col='black',lwd=1.5,lty=2)
      abline(v=min(obs.HRU),
             col='red',lwd=2,lty=2)
      abline(v=max(wgen.trace),
             col='black',lwd=1.5,lty=2)
      abline(v=max(obs.HRU),
             col='red',lwd=2,lty=2)
      if (my.opt==3){
        legend("top",
               legend=c('Sim [min,median,max]','Obs [min,median,max]'),
               col=c('gray','red'),
               cex=2.5,
               fill=c('gray','red'),
               border = c(NA,NA),
               horiz = FALSE,inset=0.01)
        mtext(paste0('at-basin'),4,adj=0.925, cex=2,col='purple')
      }
    }
    
    dev.off()
  }
  
  print(paste0("-->> figures were saved at: ", './Figures/'))
  
  # The end #
}