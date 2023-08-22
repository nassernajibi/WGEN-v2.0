

##########################Linear Program Attempt 1####################
WR_prob_estimator <- function(WR_prob_chunk,WR_prob_change,
                              threshold,piecewise_limit,
                              weather.nonpar.state.assignments) {

  ####notation####
  # i = WR (1,...,K)
  # j = time period (1,...,n)
  
  ####key constraints####
  #the sum of the probabilities for sampling different time periods is unity 
  #p1+...+pn = 1        
  
  #for each WR (i), the weighted sum of the frequencies of WR[i] across time periods should be as close as possible to the target
  #p1*theta[i,1] + ... + pn*theta[i,n] - above[i] - below[i] = target_theta[i]
  
  ####objective function###
  #minimize:
  #the degree to which the weighted sum of frequencies deviate from the targets across WRs, 
  #plus the degree to which the sampling probabilities deviate from uniform probabilities (with a piecewise cost curve)
  #
  #min 
  #sum(above[i] + below[i]) + sum(above[i] + below[i]) +                             #WR targets
  #sum(C1*above_p1_1[j] + C2*above_p1_1[j] + C1*below_p1_1[j] + C2*below_p1_1[j])    #uniform probability targets
  
  #decision variables
  #p1 = probability of selecting the 1st time period in the sampling
  #...
  #pn = probability of selecting the nth time period in the sampling
  #
  #
  #above1 = amount above target probability for WR1
  #below1 = amount below target probability for WR1
  #...
  #aboveK = amount above target probability for WRK
  #belowK = amount below target probability for WRK
  #
  #above_p1 = amount that p1 is above the uniform probability of 1/num_periods
  #below_p1 = amount that p1 is below the uniform probability of 1/num_periods
  #...
  #above_pn = amount that p1 is above the uniform probability of 1/num_periods
  #below_pn = amount that p1 is below the uniform probability of 1/num_periods
  #
  #above_p1_1 = the component of above_p1 that gets the smaller penalty in the piecewise cost
  #above_p1_2 = the component of above_p1 that gets the larger penalty in the piecewise cost
  #below_p1_1 = the component of below_p1 that gets the smaller penalty in the piecewise cost
  #below_p1_2 = the component of below_p1 that gets the larger penalty in the piecewise cost
  #...
  #above_pn_1 = the component of above_pn that gets the smaller penalty in the piecewise cost
  #above_pn_2 = the component of above_pn that gets the larger penalty in the piecewise cost
  #below_pn_1 = the component of below_pn that gets the smaller penalty in the piecewise cost
  #below_pn_2 = the component of below_pn that gets the larger penalty in the piecewise cost  
  
  #dimensions of WR chunks
  num_periods <- ncol(WR_prob_chunk)
  num_WRs <- nrow(WR_prob_chunk)
  
  #get climatological probabilities
  WR_prob_avg <- table(weather.nonpar.state.assignments)/length(weather.nonpar.state.assignments)
  
  #calculate target probabilities for WRs
  target_probabilities <- WR_prob_avg*(1+WR_prob_change)
  target_probabilities <- target_probabilities/sum(target_probabilities)   #need to rescale to sum to 1
  
  #these are the target (uniform) probabilities for each sampling period
  uniform_prob <- 1/num_periods
  

  #cost coefficients associated with decision variables
  penalty_target <- 1
  penalty_uniform1 <- 5
  penalty_uniform2 <- 100
  C <- c(rep(0,num_periods),                                       #no cost on probabilities for different time periods
         rep(penalty_target,2*num_WRs),                            #equal costs for being above and below target WR probabilities
         rep(0,2*num_periods),                                     #we place no cost on the decision variables of being above and below uniform probabilities,
         rep(c(penalty_uniform1,penalty_uniform2),2*num_periods)   #but instead place a piecewise linear cost on their components
         )
  
  # Create constraint matrix A and right hand side of constraints
  A1 <- WR_prob_chunk                     #first part of constraints relates to p1*theta[i,1] + ... + pn*theta[i,n]
  A2 <- array(0,c(num_WRs,2*num_WRs))     #second part of constraints relate to p1*theta[i,1] + ... + pn*theta[i,n]
  A3 <- array(0,c(num_WRs,2*num_periods))     #second part of constraints relate to p1*theta[i,1] + ... + pn*theta[i,n]
  A4 <- array(0,c(num_WRs,4*num_periods))     #second part of constraints relate to p1*theta[i,1] + ... + pn*theta[i,n]
  A <- cbind(A1,A2,A3,A4)
  
  #create constraints for the target probabilities of each WR
  for (j in 1:num_WRs) {
    A[j,(num_periods + 2*(j-1)+1):(num_periods + 2*(j-1)+2)] <- c(-1,1)  #-above, +below
  }
  B <- c(target_probabilities)
  constraints_direction  <- rep("=",length(B))
  
  #constraint so all probabilities sum to 1
  A <- rbind(A,c(rep(1,num_periods),rep(0,2*num_WRs),rep(0,2*num_periods),rep(0,4*num_periods)))    
  B <- c(B,1)
  constraints_direction  <- c(constraints_direction,"=") 
  
  #require that the above and below target decision variables associated with the constraints above are all below some threshold
  for (j in 1:num_WRs) {
    aboves <- rep(0,num_WRs*2)
    aboves[2*(j-1)+1] <- 1
    cur_set <- c(rep(0,num_periods),aboves,rep(0,2*num_periods),rep(0,4*num_periods))
    A <- rbind(A,cur_set)
    B <- c(B,threshold)
    constraints_direction <- c(constraints_direction,"<=")
    
    belows <- rep(0,num_WRs*2)
    belows[2*(j-1)+2] <- 1
    cur_set <- c(rep(0,num_periods),belows,rep(0,2*num_periods),rep(0,4*num_periods))
    A <- rbind(A,cur_set)
    B <- c(B,threshold)
    constraints_direction <- c(constraints_direction,"<=")
  }
  
  #now add a series of constraints setting deviations of the probabilities above and below the uniform probability
  start <- num_WRs*2 + num_periods + 1
  for (j in 1:num_periods) {
    cur_con <- rep(0,ncol(A))
    cur_con[j] <- 1
    cur_con[(start + 2*(j-1)):(start + 2*(j-1) + 1)] <- c(-1,1) #-above, +below
    A <- rbind(A,cur_con)
    B <- c(B,uniform_prob)
    constraints_direction <- c(constraints_direction,"=")
  }
  
  #now add definitions of the piecewise components of above and below deviations from the uniform probabilities
  start <- num_WRs*2 + num_periods + 1
  start2 <- num_WRs*2 + num_periods + 2*num_periods + 1
  for (j in 1:num_periods) {  
    cur_con_above <- rep(0,ncol(A))
    cur_con_above[(start + 2*(j-1))] <- 1
    cur_con_above[(start2 + 4*(j-1)):(start2 + 4*(j-1) + 1)] <- c(-1,-1)
    A <- rbind(A,cur_con_above)
    B <- c(B,0)
    constraints_direction <- c(constraints_direction,"=")
    
    cur_con_below <- rep(0,ncol(A))
    cur_con_below[(start + 2*(j-1) + 1)] <- 1
    cur_con_below[(start2 + 4*(j-1)+2):(start2 + 4*(j-1) + 3)] <- c(-1,-1)
    A <- rbind(A,cur_con_below)
    B <- c(B,0)
    constraints_direction <- c(constraints_direction,"=")
  }
  
  #finally, we constrain the first part of the above_p and below_p components of the piecewise to be less than piecewise_limit 
  start2 <- num_WRs*2 + num_periods + 2*num_periods + 1
  for (j in 1:num_periods) {  
    cur_piecewise_above <- rep(0,ncol(A))
    cur_piecewise_above[(start2 + 4*(j-1))] <- 1
    A <- rbind(A,cur_piecewise_above)
    B <- c(B,piecewise_limit)
    constraints_direction <- c(constraints_direction,"<=")
    
    cur_piecewise_below <- rep(0,ncol(A))
    cur_piecewise_below[(start2 + 4*(j-1) + 2)] <- 1
    A <- rbind(A,cur_piecewise_below)
    B <- c(B,piecewise_limit)
    constraints_direction <- c(constraints_direction,"<=")
    
  }
  
    
  # Find the optimal solution
  optimum <-  lp(direction="min",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constraints_direction,
                 const.rhs = B
  )
  
  # Print status: 0 = success, 2 = no feasible solution
  optimum$status
  
  # Display the optimum values 
  sampling_probabilities <- optimum$solution[1:num_periods]
  deviation_from_WR_target <- optimum$solution[(num_periods+1):(num_periods+2*num_WRs)]
  deviation_from_prob_unif <- optimum$solution[(num_periods+2*num_WRs+1):(num_periods+2*num_WRs+2*num_periods)]
  deviation_from_prob_unif_piecewise <- optimum$solution[(num_periods+2*num_WRs+2*num_periods+1):length(optimum$solution)]
  
  # Check the value of objective function at optimal point
  optimum$objval
  
  output <- list(optimum$status,sampling_probabilities,deviation_from_WR_target,deviation_from_prob_unif,deviation_from_prob_unif_piecewise,optimum$objval)
  
  return(output)
}
#################################################





