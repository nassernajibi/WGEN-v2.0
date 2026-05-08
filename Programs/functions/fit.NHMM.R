
fit.NHMM <- function(my.nstates,
                     my.synoptic.pcs,
                     my.covar.trans,
                     fo,
                     n.eofs) {
  
  #'  ---------- # s-NHMM for Annual WRs Simulations -------------------------------
  #'  Fits a non-homogeneous HMM with Seasonality using depmixs4 package and 
  #'  simulates markov chains and response models with Seasonality 
  #'  
  #'  Arguments: 
  #' 
  #'  my.nstates : Number of hidden states of Weather Regimes
  #'  my.synoptic.pcs : "Observed" series of EOFs 
  #'  my.covar.trans : a matrix of covariates to use in the multinomial regression of the NHMM
  #'  fo : the formula for the multinomial regression of the NHMM
  #'  n.eofs : number of PCs of geopotential height used (i.e., # of response variables) in the NHMM

  # A list of transition models, each created by a call to transInit. The length
  # of this list should be the number of states of the model.
  # PRE-ALLOCATE: Create transition models with `my.nstates` size
  transition <- vector("list", length = my.nstates) 
  for(i in 1:my.nstates){
    transition[[i]] <- transInit(formula = fo, 
                                 data = my.covar.trans, 
                                 nstates = my.nstates)
  }  
  
  my.prior <- transInit(~1,
                        ns = my.nstates,
                        data=data.frame(1),
                        ps = runif(my.nstates)
                        )
  
  my.response.models <- vector("list", length = my.nstates)
    for(i in 1:my.nstates){
      my.response.models[[i]] <- list()
      for (eof in 1:n.eofs) {
        fo.eof <- as.formula(paste0("PC",eof, "~1"))
        my.response.models[[i]][[eof]] <- GLMresponse(formula = fo.eof,
                                                      data = data.frame(my.synoptic.pcs ),
                                                      family = gaussian())
      }
    }  


  # create the model
  mod <- makeDepmix(
    response = my.response.models,
    transition = transition,
    prior = my.prior,
    ntimes = nrow(my.synoptic.pcs),
    homogeneous = FALSE
  )
  
  
  
  ########-------- Model fitting ------##############
  
  # fit model n.attempts times, pick the best one
  n.attempts <- 10
  tmp.mod.list <- vector("list", length = n.attempts)
  for(j in 1:n.attempts){ 
    cat('\n NHMM fitting attempt #: ', j, '\n')
    tmp.mod.list[[j]] <- fit(mod,emc = em.control(random = TRUE),
                             verbose = F)
  }
  
  
  # check for convergence  
  all.msgs <- sapply(tmp.mod.list, function(x) {
    if (class(x) == "depmix.fitted") {
      return(x@message)
    } else {
      return(c())
    }
  })
  
  # identify which ones converged, use `grepl` instead of `stringr::str_match`
  index.converged <- grepl( "Log likelihood converged", all.msgs)
  
  # find log likelihood for all models, and set to -Inf for those that didn't converge
  logLike.list <- as.numeric(unlist(lapply(tmp.mod.list, logLik)))
  logLike.list[!index.converged] <- -Inf
  
  # find the best model (i.e., the smallest NEGATIVE log likelihood)
  # that's why `which.max()` is used instead of `which.min()`
  mod.num <- which.max(logLike.list) 
  cat("Selected model", mod.num, "with log-likelihood:", logLike.list[mod.num], "\n")
  
  fmod.depmix <- tmp.mod.list[[mod.num]]

  # get historical state sequence
  state.seq <- posterior(fmod.depmix)$state

  invisible(gc())
    
  return(list(fitted.model = fmod.depmix,
              viterbi.seq = state.seq))
  
}





