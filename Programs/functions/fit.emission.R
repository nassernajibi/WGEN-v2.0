
fit.emission <- function(prcp.site,months,months.weather,n.sites,thshd.prcp) {
  
  #this function fits gamma and GPD distributions to the prcp.site data, split by a threshold
  
  #prcp.site = a matrix with precipitation through time (rows) across sites (columns)
  #months = a vector of unique months over which data is available
  #months.weather = a time series of months of the same length as prcp.site
  #n.sites = the number of sites
  #thshd.prcp = the threshold for non-zero precipitation that separates the gamma and GPD distributions
  
  #identify indexes for below/above threshold per site for gamma and gpd fits  
  idx1 <- sapply(1:n.sites,function(x,v){
    (prcp.site[,x]!=0 & !is.na(prcp.site[,x]) & prcp.site[,x]<=v[x])
  },v=thshd.prcp)
  idx2 <- sapply(1:n.sites,function(x,v){
    (prcp.site[,x]!=0 & !is.na(prcp.site[,x]) & prcp.site[,x]>v[x])
  },v=thshd.prcp)
  
  #1: gamma
  prcp.site11 <- prcp.site*idx1 # tagging those ones for fitting gamma (0 or values below threshold)
  emission.fit.site1 <- apply(prcp.site11,2,function(x,m,mm) {
    sapply(m,function(m) {fitdistr(as.numeric(x[x!=0 & mm==m & !is.na(x)]),'gamma')$estimate})
  },m=months, mm=months.weather)
  
  dim(emission.fit.site1) <- c(2,length(months),n.sites)   #reformat dimensions so that shape/rate down each row, months in each column, 3rd dimension indexes the sites
  
  #2: gpd
  prcp.site22 <- prcp.site*idx2 # tagging those ones for fitting gpd (0 or values above threshold)
  emission.fit.site2 <- sapply(1:n.sites,function(x) {
    y <- prcp.site22[,x]
    my.fit <- eva::gpdFit(as.numeric(y[y!=0 & !is.na(y)]),threshold = thshd.prcp[x])
    return(my.fit$par.ests)
    })
  
  emission.fits.site <- list("Gamma"=emission.fit.site1,
                             "GPD"=emission.fit.site2)
  return(emission.fits.site)
}