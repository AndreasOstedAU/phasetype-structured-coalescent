library(partitions)

#HELPER FUNCTION 1: Generate the state space 
state_space = function(n, D) {
  #n: sample size
  #D: number of subpopulations
  states = list()
  for (i in 2:n){
    
    #gives all the vectors of length D such that their sum is i
    states[[i-1]] <- t(compositions(i, D)) 
  }
  states <- rev(states) #reverse, so higher n at top 
  states_matrix <- do.call(rbind, states) #without absorbing state, binding rows of lists into matrix
  
  states_matrix_full <- rbind(states_matrix, Inf) #adding absorbing state
  return(states_matrix_full)
}



#HELPER FUNCTION 2: Get the type of transition in the markov process
get_transition_type = function(alpha, beta) {
  #alpha and beta are configurations of the sample, i.e. vectors of length D giving the number of lineages in each subpopulation.
  diff = alpha - beta
  
  #If no change => self transition
  if (all(diff == 0)) {
    index = NA
    type = 'self'
  }
  
  #Total of 1 change, difference is 1 => coalescent
  else if ((sum(diff != 0) == 1) & sum(diff == 1)) {
    index = which(diff == 1) #Index of subpopulation coalescent happened in.
    type = 'coalescent'
  }
  
  # If no change is greater than 1, exactly two indexes (subpopulations) change
  # and the sums over the vector entries are the same => migration
  else if ( (all(abs(diff) <= 1)) && (length(diff[diff != 0]) == 2) && (sum(alpha) == sum(beta))) {
    #vector of length 2: index of subpop. migration happens FROM, index of subpop. migration happens TO
    index = c(which(diff == 1),which(diff == -1)) 
    type = 'migration'
  }
  
  #transition to absorbing state (MRCA)
  # diff will be Inf by construction. 
  # check that it only happens from a state with exactly one subpop. of 2. 
  else if ( all(is.infinite(diff)) && (sum(alpha) == 2) && (length(alpha[alpha] != 0) == 1)) {
    index = which(alpha == 2)
    type = 'absorbtion'
  }
  
  else {
    index = NA
    type = 'not possible'
  }
  return(list('sub_pop' = index, 'type' = type))
}



# Function to generate subintensity matrix for a coalescent model with structure, given sample size n and # subpopulations D:
# Implements formula (5.22) in Wakeley, 2009.
generate_subint_matrix = function(n, D, subpop_sizes = rep(1, D), M = matrix(rep(1/(D-1), D*D), nrow = D, ncol = D)) {
  
  # Optional arguments:
    #subpop_sizes: optional argument giving relative sizes of the subpopulations, 
      # Defaults to 1 for each subpop. (i.e. same sizes)
    #M: optional matrix of size D x D giving rates of migration between each subpopulation. 
      # Defaults to 1/D-1 in each entry, i.e. migration rate is 1 between all subpops for D = 2. 

  # generate states 
  states_matrix_full = state_space(n, D)
  
  n_states = dim(states_matrix_full)[1]
  
  matrix_out = matrix(rep(NA, n_states^2), nrow = n_states)
  
  for (i in 1:(n_states-1)) {
    for (j in 1:n_states) {
      
      transition <- get_transition_type(states_matrix_full[i,], states_matrix_full[j,])
      trans_type <- transition$type
      subpop_index <- transition$sub_pop 
      
      # alpha_i and c_i(size and relative size of subpopulation i) only available if subpop_index is not NA, 
      # i.e. if transition is coalescent or migration
      alpha_i <- ifelse(is.na(subpop_index[1]), NA, states_matrix_full[i,][subpop_index[1]])
      c_i <- ifelse(is.na(subpop_index[1]), NA, subpop_sizes[subpop_index[1]])
      
      # fill entries of matrix using formula (5.22) in Wakeley, 2009
      if (trans_type == 'coalescent') {
        matrix_out[i,j] =  choose(alpha_i,2)/c_i
      }
      else if (trans_type == 'migration') {
        matrix_out[i,j] = alpha_i*M[subpop_index[1],subpop_index[2]]/2
        
      }
      else if (trans_type == 'absorbtion') {
        matrix_out[i,j] = 1/c_i
      }
      else {
        matrix_out[i,j] = 0
      }
    }
  }
  
  #set diagonal elements to minus rowsums (works since diagonals are set to 0 initially)
  diag(matrix_out) <- -rowSums(matrix_out) 
  
  #take subintensity part of matrix
  subint_mat = matrix_out[1:n_states-1, 1:n_states-1] 
  
  #last col = exit vector (excluding row with absorbing state)
  exit_vector = matrix_out[1:n_states-1, n_states] 
  return(list('subint_mat' = subint_mat, 'exit_vector' = exit_vector))
}
