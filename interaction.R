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
  
  return(c(f_sens, t_sens, I))
}

interaction_mc <- function(N, J, d, gpi, seed = 1){
  sens = foreach::foreach(i = 1:J, .combine = "rbind", .packages = "hetGP") %dopar% {
    interaction(N, d, gpi, seed + i)
  }
  return(sens)
}
