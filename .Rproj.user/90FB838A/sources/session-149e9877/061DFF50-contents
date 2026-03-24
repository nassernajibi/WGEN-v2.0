quantile.mapping <- function(prcp.site,Sbasin,thshd.prcp,perc.q,emission.old1,emission.old2,emission.new1,n.sites,months,dates.sim,cur.jitter) {
  
  #this function takes simulated precipitation at n.sites sites, and then adjusts the values through quantile mapping.
  #in the mapping "new" (i.e., adjusted) gamma distributions for each site are used (these could be CC-scaled from original).
  #also, there is an option to adjust the non-exeedance probabilities using the conditional variance of the sites based on the basin average
  
  #Arguments:
  #(1): gamma; (2): gpd
  #prcp.site =  a matrix (n.sim x n.sites) of simulated prcp values for the sites 
  #Sbasin = the spearman correlation matrix (n.site x nsite) between the precipitation sites
  #thshd.prcp = threshold for gpd fit
  #perc.q = the percentage to scale the upper tail of the distributions
  #emission.old1 = 'gamma', an (n.sim x n.sites x n.parameters) matrix containing time series of original emission parameters for each site
  #emission.old2 = 'gpd', an (n.sim x n.sites x n.parameters) matrix containing time series of original emission parameters for each site
  #emission.new1 = 'gamma', an (n.sim x n.sites x n.parameters) matrix containing time series of adjusted emission parameters for each site
  #n.sites = the number of individual sites
  #months = a vector of calendar months for the season
  #dates.sim = a time series of the simulated dates
  #cur.jitter = TRUE/FALSE indicating whether to perturb the non-exceedance probabilities
  
  months.sim <- as.numeric(format(dates.sim,"%m")) # t-by-1
  
  #create the conditional correlation matrix for the sites
  std.S.cond <- diag(rep(0.4,n.sites)) # lambda == 0.1,...,0.9 etc here (0.4)
  SIGMA <- std.S.cond%*%Sbasin%*%t(std.S.cond)
  S.cond <- SIGMA
  
  n <- dim(prcp.site)[1]
  
  idx1 <- sapply(1:n.sites,function(x,v){
    which(prcp.site[,x]!=0 & !is.na(prcp.site[,x]) & prcp.site[,x]<v[x])
  },v=thshd.prcp)
  idx2 <- sapply(1:n.sites,function(x,v){
    which(prcp.site[,x]!=0 & !is.na(prcp.site[,x]) & prcp.site[,x]>=v[x])
  },v=thshd.prcp)
  months.sim.1 <- sapply(1:n.sites,function(x) {months.sim[idx1[[x]]]
  })
  months.sim.2 <- sapply(1:n.sites,function(x) {months.sim[idx2[[x]]]
  })
  
  #//(1) gamma:
  #calculate the U's using old emission (1) distribution for precip less than threshold
  prcp.site.u1 <- sapply(1:n.sites,function(x,z,m,mm) {
    pgamma(prcp.site[idx1[[x]],x],'shape'=z[1,match(m[idx1[[x]]],mm),x],'rate'=z[2,match(m[idx1[[x]]],mm),x])
  }, z=emission.old1,m=months.sim,mm=months
  )

  #replace any probability==1 (to avoid getting 'inf' instances in next 'sapply') with a very close value to 1.
  check.probs.one <- sapply(1:n.sites,function(x) {sum(prcp.site.u1[[x]]==1)>0})
  if (sum(check.probs.one==1)>0){
    idx.ones <- lapply(prcp.site.u1,function(x){which(x==1)})
    prcp.site.u1 <- sapply(1:n.sites,function(x) {
      replace(prcp.site.u1[[x]],idx.ones[[x]],1-.Machine$double.eps)
    })
  }
  
  #transform new U's to new precipitation using new emission distributions
  prcp.site.new1 <- sapply(1:n.sites,function(x,z,m,mm) {
    qgamma(prcp.site.u1[[x]],'shape'=z[1,match(m[[x]],mm),x],'rate'=z[2,match(m[[x]],mm),x])
  },z=emission.new1,m=months.sim.1,mm=months
  )
  #round to limit memory size
  prcp.site.new1 <- sapply(1:n.sites,function(x) {round(prcp.site.new1[[x]],2)})
  
  
  #//(2) gpd:
  #calculate the U's using old emission (2) distribution for precip greater than threshold
  prcp.site.u2 <- sapply(1:n.sites,function(x,z) {
    evmix::pgpd(prcp.site[idx2[[x]],x],'xi'=z[2,x],'sigmau'=z[1,x],'u'=thshd.prcp[x])
  }, z=emission.old2)
  
  prcp.site.u2.new <- prcp.site.u2
  if (cur.jitter) {   # if we want to perturb U's, jump into Z space to do the jittering :)
    #then get Zs
    prcp.site.z2 <- lapply(prcp.site.u2,function(x) {qnorm(x)})
    #simulate new Zs
    jitter.samp <- rmvnorm(n,rep(0,n.sites),S.cond)
    jitter.samp2 <- sapply(1:n.sites,function(x) {jitter.samp[idx2[[x]],x]})
    
    prcp.site.z2.proposed <- sapply(1:n.sites,function(x) {prcp.site.z2[[x]]+jitter.samp2[[x]]})

    #transform new Z's to new U's
    prcp.site.u2.proposed <- lapply(prcp.site.z2.proposed,function(x){pnorm(x)})

    #accept proposed u with probability proportional to ratio of old and new exceedance ratios
    set.seed(1)
    samp.cf <- rnorm(n,0,1)
    samp2.cf <- sapply(1:n.sites,function(x) {samp.cf[idx2[[x]]]})
    coin.flip <- lapply(samp2.cf,function(x){pnorm(x)})
    
    prcp.site.u2.new <- sapply(1:n.sites,function(i,x,y,cf) {
      final <- y[[i]]

      lower <- which(x[[i]]<y[[i]])
      chance.of.replace <- x[[i]][lower]/y[[i]][lower]
      replace <- cf[[i]][lower]<chance.of.replace
      final[lower][replace] <- x[[i]][lower][replace]

      higher <- which(x[[i]]>y[[i]])
      chance.of.replace <- (1-x[[i]][higher])/(1-y[[i]][higher])
      replace <- cf[[i]][higher]<chance.of.replace
      final[higher][replace] <- x[[i]][higher][replace]

      return(final)
    },x=prcp.site.u2.proposed,y=prcp.site.u2,cf=coin.flip)
  }
  
  #replace any probability==1 (to avoid getting 'inf' instances in next 'sapply') with a very close value to 1.
  check.probs.one <- sapply(1:n.sites,function(x) {sum(prcp.site.u2.new[[x]]==1)>0})
  if (sum(check.probs.one==1)>0){
    idx.ones <- lapply(prcp.site.u2.new,function(x){which(x==1)})
    prcp.site.u0.new <- sapply(1:n.sites,function(x) {
      replace(prcp.site.u2.new[[x]],idx.ones[[x]],1-.Machine$double.eps)
    })
    prcp.site.u2.new <- prcp.site.u0.new
  }
  
  #transform new U's to new precipitation using old emission distributions then multiply the values by (1.07)^d.T
  prcp.site.new2 <- sapply(1:n.sites,function(x,z) {
    perc.q*evmix::qgpd(prcp.site.u2.new[[x]],'xi'=z[2,x],'sigmau'=z[1,x],'u'=thshd.prcp[x])
  },z=emission.old2)
  #round to limit memory size
  prcp.site.new2 <- sapply(1:n.sites,function(x) {round(prcp.site.new2[[x]],2)})
  
  
  #// place the estimated precip back to the original time-series based on their indexes
  prcp.site.new <- prcp.site
  for(i in 1:n.sites){
    prcp.site.new[idx1[[i]],i] <- prcp.site.new1[[i]]
    prcp.site.new[idx2[[i]],i] <- prcp.site.new2[[i]]
  }
  return(prcp.site.new)
  
}