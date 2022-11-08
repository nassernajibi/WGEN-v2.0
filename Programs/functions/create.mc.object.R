
create.mc.object <- function(my.NHMM,
                             my.nstates,
                             my.covar.trans.sim) {
  
  #'  ---------- # s-NHMM for Annual WRs Simulations -------------------------------
  #'  Simulate from a non-homogeneous HMM with Seasonality  using depmixs4 package
  #'  
  #'  Arguments: 
  #' 
  #'  my.NHMM = a fitted NHMM
  #'  my.nstates : Number of hidden states of Weather Regimes
  #'  my.covar.trans.sim : covariates to drive the NHMM under simulation
  
  
  #number of days in the coviariates for simulation
  n.sim=dim(my.covar.trans.sim)[[1]]
  #final transition probability matrices for each of the n.sim days, where each matrix (dimension my.nstates X my.nstates) for each day is stored in a list
  final.transition.matrix <- list()
  
  #start the loop over n.sim days
  for (j in 1:n.sim){
    #transition probability matrix for day j, which will be square and based on the number of states 
    transition_matrix=matrix(NA, nrow = my.nstates, ncol = my.nstates)
    
    #look through the states 
    for (i in 1:my.nstates) {
      
      # for each state, we find the multinomial regression parameters that define the probaiblity of going from state i to the other states (columns of the coef matrix), based on the covariates (rows of the coef matrix)
      my.coef <- my.NHMM@transition[[i]]@parameters$coefficients
      
      #probabilities (before normalization) of going from state i to state h, based on multinomial regression and coviarates on day j
      for (h in 1:my.nstates) {
        transition_matrix[i,h] <- sum(my.coef[,h]*my.covar.trans.sim[j,])
      }
      #normalize into actual probabilities by dividing by sum of values by row
      transition_matrix[i,]=exp(transition_matrix[i,]) / sum(exp(transition_matrix[i,]))
    } 
    #save for day j of n.sim total days
    final.transition.matrix[[j]]=transition_matrix
  }
  
  #create markov chain object with time-varying transition probabilities
  mcObject.time.varying <- mclapply(X=1:n.sim,mc.preschedule=TRUE,mc.cores=1,FUN=function(t){
    tr.prob=as.matrix(final.transition.matrix[[t]])
    mcObject.time.varying.out <- new("markovchain", states = paste(1:my.nstates),
                                     transitionMatrix = tr.prob, name = paste("McObject",t,sep=""))    
    return(mcObject.time.varying.out)
  })
  
  mcObject.final <- new("markovchainList",markovchains = mcObject.time.varying, name = "mcObject.nh")
  
  
  return(mcObject.final)
  
}
