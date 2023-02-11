library(foreach)
library(MASS)
#numCores <- detectCores()
#numCores

interaction <- function(N, d, gpi, seed = 1, uk = NULL, p = 0, fixed_val = NULL){
  set.seed(seed)
  M = lhs::randomLHS(N, d)
  Mprime = lhs::randomLHS(N, d)
  
  # map the LHS points to the ones with margin defined by uk
  if(!is.null(uk)){
    for(i in 1:d){
      if(uk$dist[d] == "normal"){
        M[,d] = qnorm(M[,d], uk$param[d,1], uk$param[d,2])
        Mprime[,d] = qnorm(Mprime[,d], uk$param[d,1], uk$param[d,2])
      }else if(uk$dist[d] == "lognormal"){
        M[,d] = qlnorm(M[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
        Mprime[,d] = qlnorm(Mprime[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
      }else if(uk$dist[d] == "bernoulli"){
        M[,d] = ifelse(M[,d] < uk$param[d,1], 0, 1)
        Mprime[,d] = ifelse(Mprime[,d] < uk$param[d,1], 0, 1)
      }else{
        stop("Only support normal, longnormal and bernoulli")
      }
    }
  }
  # If the number of covariates p is larger than 0, set the last p columns to a 
  # fixed value so that they will not join the computation of sensitivity indices.
  if(p > 0){
    for(i in 1:p){
      M = cbind(M,rep(fixed_val[i],N))
      Mprime = cbind(Mprime,rep(fixed_val[i],N))
    }
  }
  # fist-order sensitivity
  f_sens = first_order_sens(gpi, d, M, Mprime)
  
  # total sensitivity
  t_sens = total_sens(gpi, d, M, Mprime)
  
  # difference between first-order and total sensitivity scores
  I = t_sens - f_sens
  I[I < 0] = 0
  
  return(cbind(f_sens, t_sens, I))
}

interaction_mult <- function(N, d, gpi, idx, seed = 1, uk = NULL, p = 0, 
                             fixed_val = NULL){
  set.seed(seed)
  M = lhs::randomLHS(N, d)
  Mprime = lhs::randomLHS(N, d)
  if(!is.null(uk)){
    for(i in 1:d){
      if(uk$dist[d] == "normal"){
        M[,d] = qnorm(M[,d], uk$param[d,1], uk$param[d,2])
        Mprime[,d] = qnorm(Mprime[,d], uk$param[d,1], uk$param[d,2])
      }else if(uk$dist[d] == "lognormal"){
        M[,d] = qlnorm(M[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
        Mprime[,d] = qlnorm(Mprime[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
      }else if(uk$dist[d] == "bernoulli"){
        M[,d] = ifelse(M[,d] < uk$param[d,1], 0, 1)
        Mprime[,d] = ifelse(Mprime[,d] < uk$param[d,1], 0, 1)
      }else{
        stop("Only support normal, longnormal and bernoulli")
      }
    }
  }
  if(p > 0){
    for(i in 1:p){
      M = cbind(M,rep(fixed_val[i],N))
      Mprime = cbind(Mprime,rep(fixed_val[i],N))
    }
  }
  # fist-order sensitivity
  f_sens = first_order_sens_mult(gpi, d, M, Mprime,idx)
  
  # total sensitivity
  t_sens = total_sens_mult(gpi, d, M, Mprime,idx)
  
  # difference between first-order and total sensitivity scores
  I = t_sens - f_sens
  I[I < 0] = 0
  
  return(cbind(f_sens, t_sens, I))
}

interaction_thres <- 
  function(X, sd, N, d = ncol(X), upper = 30, seed = 1, uk = NULL){
  # The measurement error is tau square times nugget g from Section 5.2.2 in the 
  # Book "Surrogate" by Grammacy
  # All variables have no interaction and the first column is null
  n = nrow(X)
  set.seed(seed)
  Y_noInter = rowSums(X[,1:d]) + rnorm(n,0,sd)
  gpi_perm = mleHomGP(X = as.matrix(X[,1:d]), Z = Y_noInter, lower = rep(1, d), 
                      upper = rep(upper, d), covtype = 'Gaussian', maxit = 500)
  while(gpi_perm$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
    upper = upper + 1
    gpi_perm = mleHomGP(X = as.matrix(X[,1:d]), Z = Y_noInter, lower = rep(1, d), 
                        upper = rep(upper, d), covtype = 'Gaussian', maxit = 500)
  }
  sens = interaction(N, d, gpi_perm, seed = seed, uk = uk)
  return(sens[,3])
}

##### ------------ Multiple-run of interaction ------------- ####
interaction_mc <- function(N, J, d, gpi, seed = 1, uk = NULL, p = 0,
                           fixed_val = NULL){
  interaction <- function(N, d, gpi, seed = 1, uk = NULL, p = 0,
                          fixed_val = NULL){
    set.seed(seed)
    source("sensitivity_hetGP.R")
    M = lhs::randomLHS(N, d)
    Mprime = lhs::randomLHS(N, d)
    if(!is.null(uk)){
      for(i in 1:d){
        if(uk$dist[d] == "normal"){
          M[,d] = qnorm(M[,d], uk$param[d,1], uk$param[d,2])
          Mprime[,d] = qnorm(Mprime[,d], uk$param[d,1], uk$param[d,2])
        }else if(uk$dist[d] == "lognormal"){
          M[,d] = qlnorm(M[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
          Mprime[,d] = qlnorm(Mprime[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
        }else if(uk$dist[d] == "bernoulli"){
          M[,d] = ifelse(M[,d] < uk$param[d,1], 0, 1)
          Mprime[,d] = ifelse(Mprime[,d] < uk$param[d,1], 0, 1)
        }else{
          stop("Only support normal, longnormal and bernoulli")
        }
      }
    }
    
    if(p > 0){
      for(i in 1:p){
        M = cbind(M,rep(fixed_val[i],N))
        Mprime = cbind(Mprime,rep(fixed_val[i],N))
      }
    }
    
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
    interaction(N, d, gpi, seed + i, uk = uk, p = p, fixed_val = fixed_val)
  }
  return(sens)
}

interaction_mult_mc <- function(N,J,d,gpi,idx,seed = 1,uk =NULL, p = 0, 
                                fixed_val = NULL){
  interaction_mult <- function(N, d, gpi, idx, seed = 1, uk =NULL, p = 0,
                               fixed_val = NULL){
    set.seed(seed)
    source("sensitivity_hetGP.R")
    M = lhs::randomLHS(N, d)
    Mprime = lhs::randomLHS(N, d)
    if(!is.null(uk)){
      for(i in 1:d){
        if(uk$dist[d] == "normal"){
          M[,d] = qnorm(M[,d], uk$param[d,1], uk$param[d,2])
          Mprime[,d] = qnorm(Mprime[,d], uk$param[d,1], uk$param[d,2])
        }else if(uk$dist[d] == "lognormal"){
          M[,d] = qlnorm(M[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
          Mprime[,d] = qlnorm(Mprime[,d], meanlog = uk$param[d,1], sdlog = uk$param[d,2])
        }else if(uk$dist[d] == "bernoulli"){
          M[,d] = ifelse(M[,d] < uk$param[d,1], 0, 1)
          Mprime[,d] = ifelse(Mprime[,d] < uk$param[d,1], 0, 1)
        }else{
          stop("Only support normal, longnormal and bernoulli")
        }
      }
    }
    
    if(p > 0){
      for(i in 1:p){
        M = cbind(M,rep(fixed_val[i],N))
        Mprime = cbind(Mprime,rep(fixed_val[i],N))
      }
    }
    
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
    interaction_mult(N, d, gpi, idx, seed + i,uk = uk, p = p, fixed_val = fixed_val)
  }
  return(sens)
  
}

##### ------------ Multiple-run of interaction (single-thread version) ------------- ####
interaction_mc_single <- function(N, J, d, gpi, seed = 1, uk = NULL, p = 0,
                                  fixed_val = NULL){
  f_sens <- t_sens <- I <- rep(NA, J*d)
  for(j in 1:J){
    sens = interaction(N, d, gpi, seed + j, uk = uk, p = p, fixed_val = fixed_val)
    f_sens[(j-1)*d + (1:d)] <- sens[,1]
    t_sens[(j-1)*d + (1:d)] <- sens[,2]
    I[(j-1)*d + (1:d)] <- sens[,3]
  }
  return(cbind(f_sens,t_sens,I))
}

interaction_mult_mc_single <- function(N,J,d,gpi,idx,seed = 1, uk = NULL, p = 0,
                                       fixed_val = NULL){
  f_sens <- t_sens <- I <- rep(NA, J)
  for(j in 1:J){
    sens = interaction_mult(N, d, gpi, idx, seed + j, uk = uk, p = p, 
                            fixed_val = fixed_val)
    f_sens[j] <- sens[1]
    t_sens[j] <- sens[2]
    I[j] <- sens[3]
  }
  return(cbind(f_sens,t_sens,I))
  
}
