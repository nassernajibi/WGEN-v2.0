wgen.simulator <- function(weather.state.assignments,mc.sim,
                           prcp.basin,dates.weather,
                           first.month,last.month,dates.sim,
                           months,window.size,min.thresh) {
  
  #This function uses a simulation of weather states to conduct 
  #the primary daily weather generation, which is based on a block-bootstrap
  #procedure. It returns the simulated weather regimes it used, the resampled
  #dates from the block bootstrap, and matricies of resampled data for prcp, ...
  #at all sites in the basin.
  
  #Arguments:
  #weather.state.assignments = observed weather regimes for the historical period
  #mc.sim = a simulation of weather regimes to use in the daily simulation
  #prcp.basin = a time series of basin averaged precipitation (same length as weather.state.assignments)
  #dates.weather = a vector of dates associated with the historical record (same length as weather.state.assignments)
  #first.month = first month of the season
  #last.month = last month of the season
  #dates.sim = a vector of dates associated with the simulation (same length as mc.sim)
  #months = a vector of the calendar months that define the season (max length == 12) 
  #window.size = the window size, by calendar month, to use for bootstrapping based on days into the season
  #min.thresh = lower threshold below which day is considered dry


  
  
  
  #month sequence for the simulation period
  months.sim <- as.numeric(format(dates.sim,"%m"))
  
  n.sim <- length(dates.sim)  #length of simulation
  #################################
  
  
  #########Observed Data for Resampling##########
  #here is the observed cluster runs
  months.weather <- as.numeric(format(dates.weather,"%m"))
  clust.run.obs <- cluster.runs(my.clust=weather.state.assignments,my.date=dates.weather,first.month=first.month,last.month=last.month)
  #properties of observed clusters  
  num.clust.run.obs <- dim(clust.run.obs)[1]
  num.missing.prcp.basin <- array(NA,num.clust.run.obs)
  avg.prcp.basin <- array(NA,num.clust.run.obs)
  for (j in 1:num.clust.run.obs) {
    num.missing.prcp.basin[j] <- length(which(is.na(prcp.basin[clust.run.obs$first.run[j]:clust.run.obs$last.run[j]])))
    avg.prcp.basin[j] <- mean(prcp.basin[clust.run.obs$first.run[j]:clust.run.obs$last.run[j]],na.rm=T)
  }
  #observed lag-1 basin prcp and lag1 occurrence
  lag1.pos <- clust.run.obs$first.run-1
  lag1.pos <- lag1.pos[which(lag1.pos>=1)]
  lag1.prcp.basin.obs <- c(0,prcp.basin[lag1.pos])
  lag1.prcp.basin.obs.occur <- array(0,length(lag1.prcp.basin.obs))
  lag1.prcp.basin.obs.occur[lag1.prcp.basin.obs>min.thresh] <- 1
  
  clust.run.obs <- cbind(clust.run.obs,'lag1.prcp'=lag1.prcp.basin.obs,'lag1.prcp.occur'=lag1.prcp.basin.obs.occur,'missing'=num.missing.prcp.basin,'basin.avg'=avg.prcp.basin)
  #################################
  
  
  ###########Simulation Loop######################
  #here is the simulated cluster runs
  clust.run.sim <- cluster.runs(my.clust=mc.sim,my.date=dates.sim,first.month=first.month,last.month=last.month)
  num.clust.run.sim <- dim(clust.run.sim)[1]
  
  #initialize simulation
  sampled.date.sim <- array(NA,n.sim)
  prcp.basin.sim <- array(0,n.sim)
  
  #sample from paired observations above based on similar state sequences in simulated chain  
  for (i in 1:num.clust.run.sim) {
    #foreach (i=1:num.clust.run.sim) %dopar% {
    #initial characteristics of current WR block to be simulated
    y1 <- clust.run.sim$first.run[i]
    cur.state.sim <- clust.run.sim$state[i]
    cur.month.start.sim <- clust.run.sim$month.first[i]
    cur.month.last.sim <- clust.run.sim$month.last[i]
    cur.days.into.season.sim <- clust.run.sim$avg.day.into.season[i]
    cur.run.length.sim <- clust.run.sim$run.length[i]
    lag1.prcp.basin <- prcp.basin.sim[clust.run.sim$first.run[i]-1]
    lag1.prcp.basin.occur <- 0
    if((clust.run.sim$first.run[i]-1)<=0) {lag1.prcp.basin <- 0; lag1.prcp.basin.occur <- 0}    #for day one of simulation
    if(lag1.prcp.basin>min.thresh) {lag1.prcp.basin.occur <- 1}
    
    #this makes sure we enter the second while loop at least once,
    # AND there is 'at least one similar state' within the "entire" simulation;
    # otherwise we go to the next iteration.
    if (sum(clust.run.obs$state==cur.state.sim)>0){
      #we define the remaining length of a simulated state that requires resampling
      remaining.run.length.sim <- cur.run.length.sim
      while(remaining.run.length.sim>0) {
        #the current seasonality window to look for blocks to resample
        my.window <- window.size[which(months==months.sim[y1])]  #in days on either side of the current day of simulation
        
        #subset the data based on the state and starting month of run, and make sure there are no missing basin prcp values
        #we require that there be some observations in the subset and that there are runs with length 1 
        #so that we are guarenteed to fill remaining sim run length
        #therefore, we relax the seasonality window constraint if we need to get more samples
        cur.sub <- 0
        while(!sum(cur.sub)>0){
          cur.obs.runs <- clust.run.obs
          cur.sub <- which(cur.obs.runs$state==cur.state.sim &
                             abs(as.numeric(cur.obs.runs$avg.day.into.season)-as.numeric(cur.days.into.season.sim))%%(366-my.window) < my.window &
                             cur.obs.runs$lag1.prcp.occur==lag1.prcp.basin.occur &
                             cur.obs.runs$missing==0
          )

          #here we redefine the window to relax the constraint in the subset.
          if (!sum(cur.sub)>0){
            my.window <- my.window+3
          }
          else{
            cur.obs.runs <- clust.run.obs[cur.sub,]
          }
        }
        ########Here we choose which observed block to resample########
        #find the observed runs that are shorter than the remaining simulated run length
        #NOTE: THIS ASSUMES THAT FOR ALL MONTHS, THERE ARE STATES OF RUN LENGTH 1 TO FILL IN GAPS AS NEEDED
        order.iwd <- 1
        if (!sum(cur.sub)>0){remaining.run.length.sim <- 0}
        else{
          run.dif <-  cur.obs.runs$run.length - remaining.run.length.sim
          #randomly sample a block closest in length to the target length if there are more than 1 to choose from
          if (nrow(cur.obs.runs)>1) {
            samp.prob <- (1/(.Machine$double.eps+abs(run.dif)))^order.iwd/(.Machine$double.eps+sum((1/(.Machine$double.eps+abs(run.dif)))^order.iwd))
            my.samp <- sample(1:nrow(cur.obs.runs),size=1,replace=TRUE,prob=samp.prob)
          } else {
            my.samp <- 1
          }
          
          ########Here we update values in the simulation########
          #update remaining run length for current state run
          # condition on cutting the length if is extra
          if (run.dif[my.samp]>0)
          {
            cur.obs.runs$run.length[my.samp] <- remaining.run.length.sim
            if (rbinom(1,1,prob=c(0.5))>0){
              cur.obs.runs$last.run[my.samp] <- cur.obs.runs$last.run[my.samp]-run.dif[my.samp] #correcting last run, w.r.t. chopped-day arbitrary run-length
            } else{
              cur.obs.runs$first.run[my.samp] <- cur.obs.runs$first.run[my.samp]+run.dif[my.samp] #correcting first run, w.r.t. chopped-day arbitrary run-length
            }
            my.samp.run.length <- remaining.run.length.sim
          }
          else{
            my.samp.run.length <- cur.obs.runs$run.length[my.samp]
          }
          remaining.run.length.sim <- remaining.run.length.sim - my.samp.run.length
          #resample dates, basin prcp for the simulation
          x1 <- cur.obs.runs$first.run[my.samp]
          x2 <- cur.obs.runs$last.run[my.samp]
          x0 <- length(x1:x2) #length of segment to be filled in the dates matrix
          y2 <- y1 + x0 - 1  #ending location in simulation time series where current resample ends
          sampled.date.sim[y1:y2] <- dates.weather[x1:x2]
          prcp.basin.sim[y1:y2] <- prcp.basin[x1:x2]
          
          #update the starting location in the simulation time series for next iteration
          y1 <- y2+1
          
          #update the lag-1 and lag-2 basin prcp in the simulation in the loop for this state run
          lag1.prcp.basin <- prcp.basin.sim[y2]
          lag1.prcp.basin.occur <- 0
          if(lag1.prcp.basin>min.thresh) {lag1.prcp.basin.occur <- 1}
        }
      } #end of while loop on remaining.run.length.sim
    } #end of if condition on there is at least one with a similiar state
  } #end of for loop on num.clust.run.sim
  #################################
  
  #mean(prcp.basin.sim)
  
  #resample weather variables for all sites
  sampled.date.loc.sim <- match(as.Date(sampled.date.sim),dates.weather)
  
  return(list(mc.sim,#1
              sampled.date.sim,#2
              sampled.date.loc.sim #3
              ))
  
}