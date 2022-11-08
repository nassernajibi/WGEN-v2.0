cluster.runs <- function (my.clust,my.date,first.month,last.month) {
  
  #This function will return a data.frame with information on cluster runs, including
  #the position of the first and last day of each run, the run length, the month
  #of the first and last day of the run, and the average day of the season for the run
  
  #Arguments:
  #my.clust = a daily vector of states 
  #my.date = a daily vector of dates of the same length as my.clust
  #first.month = the first month of the season
  #last.month = the last month of the season
  
  my.month <- as.numeric(format(my.date,"%m"))
  
  n.clust <- length(my.clust)
  
  first.of.run <- c(1,which(diff(my.clust)!=0) + 1)
  last.of.run <- c(which(diff(my.clust)!=0),n.clust)
  
  #find those events that cross over a year
  cross <- which(my.month[first.of.run]==last.month & my.month[last.of.run]==first.month)
  while(length(cross)!=0) {
    cur.cross <- cross[1]
    months.of.cross <- my.month[first.of.run[cur.cross]:last.of.run[cur.cross]]
    first.step <- which(diff(months.of.cross)!=0)
    last.step <- length(months.of.cross)-first.step
    first.of.run <- c(first.of.run[1:(cur.cross-1)],first.of.run[cur.cross],first.of.run[cur.cross]+first.step,first.of.run[(cur.cross+1):length(first.of.run)])
    last.of.run <- c(last.of.run[1:(cur.cross-1)],last.of.run[cur.cross]-last.step,last.of.run[cur.cross],last.of.run[(cur.cross+1):length(last.of.run)])
    cross <- which(my.month[first.of.run]==last.month & my.month[last.of.run]==first.month)
  }  
  
  #run lengths and type
  run.length <- last.of.run-first.of.run + 1
  clust.type <- my.clust[first.of.run]
  
  days.into.season.avg <- sapply(1:length(first.of.run),function(y){
    return(yday(median(my.date[first.of.run[y]:last.of.run[y]])))
  })

  cluster.run <- data.frame('state'=clust.type,'first.run'= first.of.run, 'last.run'=last.of.run,'run.length'=run.length,
                            'month.first'=my.month[first.of.run],'month.last'=my.month[last.of.run],'avg.day.into.season'=days.into.season.avg)
  
  return(cluster.run)
  
}  