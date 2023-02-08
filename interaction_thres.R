interaction_thres_mc <- 
  function(X, sd, N, J, d = ncol(X), upper = 30, seed = 1, uk = NULL){
  sens = foreach::foreach(i = 1:J, .combine = "rbind", .packages = "hetGP") %dopar% {
    source("interaction.R")
    source("sensitivity_hetGP.R")
    sens = as.vector(interaction_thres(X, sd, N, d, upper = 30, seed = seed+i, uk = uk))
  }
  # threshold = quantile(sens, 1-alpha)
  return(sens)
}

interaction_diff_thres <- function(int_thres, B = 1000, alpha = 0.05){
  D1 = sample(int_thres, B, replace = T)
  D2 = sample(int_thres, B, replace = T)
  D12 = sample(int_thres, B, replace = T)
  
  D_delta = pmax(D1, D2) - D12
  res = quantile(D_delta, 1-alpha)
  return(res)
}
