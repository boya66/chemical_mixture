library(foreach)
library(MASS)
#numCores <- detectCores()
#numCores

interaction <- function(N, d, gpi, seed = 1){

  set.seed(seed)
  M = lhs::randomLHS(N, d)
  Mprime = lhs::randomLHS(N, d)
  # fist-order sensitivity
  f_sens = first_order_sens(gpi, d, M, Mprime)
  
  # total sensitivity
  t_sens = total_sens(gpi, d, M, Mprime)
  
  # difference between first-order and total sensitivity scores
  I = t_sens - f_sens
  I[I < 0] = 0
  
  return(cbind(f_sens, t_sens, I))
}

interaction_mult <- function(N, d, gpi, idx, seed = 1){
  set.seed(seed)
  M = lhs::randomLHS(N, d)
  Mprime = lhs::randomLHS(N, d)
  # fist-order sensitivity
  f_sens = first_order_sens_mult(gpi, d, M, Mprime,idx)
  
  # total sensitivity
  t_sens = total_sens_mult(gpi, d, M, Mprime,idx)
  
  # difference between first-order and total sensitivity scores
  I = t_sens - f_sens
  I[I < 0] = 0
  
  return(cbind(f_sens, t_sens, I))
}

interaction_thres <- function(X, sd, N, d = ncol(X), upper = 30, seed = 1){
  # The measurement error is tau square times nugget g from Section 5.2.2 on Book "Surrogate" by Grammacy
  # All variables have no interaction and the first column is null
  n = nrow(X)
  set.seed(seed)
  Y_noInter = rowSums(X[,2:d]) + rnorm(n,0,sd)
  gpi_perm = mleHomGP(X = as.matrix(X[,1:d]), Z = Y_noInter, lower = rep(1, d), upper = rep(upper, d),
                      covtype = 'Gaussian', maxit = 500)
  while(gpi_perm$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
    upper = upper + 1
    gpi_perm = mleHomGP(X = as.matrix(X[,1:d]), Z = Y_noInter, lower = rep(1, d), upper = rep(upper, d),
                        covtype = 'Gaussian', maxit = 500)
  }
  sens = interaction(N, d, gpi_perm, seed = 1)
  return(sens[,3])
}

interaction_mc <- function(N, J, d, gpi, seed = 1){
  interaction <- function(N, d, gpi, seed = 1){
    
    set.seed(seed)
    M = lhs::randomLHS(N, d)
    Mprime = lhs::randomLHS(N, d)
    # fist-order sensitivity
    f_sens = first_order_sens(gpi, d, M, Mprime)
    
    # total sensitivity
    t_sens = total_sens(gpi, d, M, Mprime)
    
    # difference between first-order and total sensitivity scores
    I = t_sens - f_sens
    I[I < 0] = 0
    
    return(cbind(f_sens, t_sens, I))
  }
  
  sens = foreach::foreach(i = 1:J, .combine = "rbind", .packages = "hetGP") %dopar% {
    interaction(N, d, gpi, seed + i)
  }
  return(sens)
}

interaction_mult_mc <- function(N,J,d,gpi,idx,seed = 1){
  interaction_mult <- function(N, d, gpi, idx, seed = 1){
    set.seed(seed)
    M = lhs::randomLHS(N, d)
    Mprime = lhs::randomLHS(N, d)
    # fist-order sensitivity
    f_sens = first_order_sens_mult(gpi, d, M, Mprime,idx)
    
    # total sensitivity
    t_sens = total_sens_mult(gpi, d, M, Mprime,idx)
    
    # difference between first-order and total sensitivity scores
    I = t_sens - f_sens
    I[I < 0] = 0
    
    return(cbind(f_sens, t_sens, I))
  }
  sens = foreach::foreach(i = 1:J, .combine = "rbind", .packages = "hetGP") %dopar% {
    interaction_mult(N, d, gpi, idx, seed + i)
  }
  return(sens)
  
}

