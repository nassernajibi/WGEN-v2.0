
fit.NHMM <- function(my.nstates,
                     my.synoptic.pcs,
                     my.covar.trans,
                     fo,
                     n.eofs) {
  
  #'  ---------- # s-NHMM for Annual WRs Simulations -------------------------------
  #'  Fits a non-homogeneous HMM with Seasonality  using 
  #'  depmixs4 package and simulates 
  #'  markov chains and response models with Seasonality 
  #'  
  #'  Arguments: 
  #' 
  #'  my.nstates : Number of hidden states of Weather Regimes
  #'  my.synoptic.pcs : "Observed" series of EOFs 
  #'  my.covar.trans : a matrix of covariates to use in the multinomial regression of the NHMM
  #'  fo : the formula for the multinomial regression of the NHMM
  #'  n.eofs : number of PCs of geopotential height used (i.e., # of response variables) in the NHMM


  
  # #initialize NHMM based on HMM
  # hhpars <- c(unlist(getpars(fit.mod.HMM)))
  # hhconMat <- fit.mod.HMM@conMat
  # init.pars <- list()
  # init.pars[['transition']] <- lapply(fit.mod.HMM@transition,
  #                                     function(x) x@parameters$coefficients)
  # init.pars[['prior']]  <- fit.mod.HMM@prior@parameters$coefficients #  prob.kmeans.list[[p]]  #
  # init.pars[['response']] <- lapply(fit.mod.HMM@response,
  #                                   function(x) lapply(x, function(y) unlist(y@parameters)))
  # init.pars[['conMat']] <- hhconMat

  
  #A list of transition models, each created by a call to transInit. The length
  #of this list should be the number of states of the model.
  transition <- list() 
  for(i in 1:my.nstates){
    transition[[i]] <- transInit(formula = fo, 
                                 data = my.covar.trans, 
                                 nstates = my.nstates)
  }  
  
  my.prior <- transInit(~1,ns = my.nstates,data=data.frame(1),
                        ps = runif(my.nstates))
  
  my.response.models <- list() 
    for(i in 1:my.nstates){
      my.response.models[[i]] <- list()
      for (eof in 1:n.eofs) {
        fo.eof <- as.formula(paste("PC",eof,"~1",sep=""))
        my.response.models[[i]][[eof]] <- GLMresponse(formula = fo.eof,data = data.frame(my.synoptic.pcs ),family = gaussian())
      }
    }  


  #create the model
  mod <- makeDepmix(response=my.response.models,
                    transition=transition,
                    prior=my.prior,
                    ntimes= nrow(my.synoptic.pcs),
                    homogeneous=FALSE)  
  
  
  
  ########---------model fitting-------##############
  
  #fit model 10 times, pick the best one
  tmp.mod.list <- list()
  for(j in 1:10){ 
    tmp.mod.list[[j]] <- fit(mod,emc = em.control(random = TRUE),verbose = F)
  }
  
  
  #check for convergence  
  all.msgs <- sapply(tmp.mod.list, function(x){
    if(class(x) == "depmix.fitted"){return(x@message)} else {return(c())}
  })
  #identify which ones converged
  index.converged <- stringr::str_match(all.msgs, "Log likelihood converged") 
  #find log likelihood for all models, and set to 99999 for those that down converge
  logLike.list <- as.numeric(unlist(lapply(tmp.mod.list,logLik)))
  logLike.list[which(is.na(index.converged))] <- 99999
  #find the best model (i.e., the smallest negative log likelihood)
  mod.num <- which.min(logLike.list) 
  fmod.depmix <- tmp.mod.list[[mod.num]]

  #get historical state sequence
  state.seq <- posterior(fmod.depmix)$state
  #############################################
  
  return(list(fitted.model = fmod.depmix,
              viterbi.seq = state.seq))
  
}
