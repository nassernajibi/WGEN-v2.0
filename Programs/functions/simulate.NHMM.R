
simulate.NHMM <- function(mcObject.final) {
  
  #'  ---------- # s-NHMM for Annual WRs Simulations -------------------------------
  #'  Simulate from a non-homogeneous HMM
  #'  
  #'  Arguments: 
  #'  mcObject.final = the final set of transition probability matricies to use for WR simulation, as stored in an mc object 

    n.sim <-length(mcObject.final@markovchains)
    mc.sim <- as.numeric(rmarkovchain(n=1,object=mcObject.final[[1]]))
    end.state <- paste(mc.sim[length(mc.sim)])
    for (t in 1:n.sim) {
      mc.sim <- c(mc.sim,as.numeric(rmarkovchain(n=1,object=mcObject.final[[t]],t0=end.state)))
      end.state <- paste(mc.sim[length(mc.sim)])
    }    
    
    #here is the final mc simulation
    final.mc.sim <- mc.sim[2:(n.sim+1)]
    return(final.mc.sim)

}

  

