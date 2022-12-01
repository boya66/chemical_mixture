#### This script contains codes for modifying the synthetic data generator by
#### voiding some variables/clusters of variables for true variable/cluster
#### importance ranking
# true function (the X8 is useless)
## Inputs: n: number of samples
##         d: dimension of the input space
##         voidIdx: the indices for the variables to void out.
## Outputs: a n*4 matrix with the cols y:
##          1st col: true function with all dimensions
##          2nd col: true function without voidIdx dimensions
##          3rd col: true function with all dimensions and measurement error
##          4th col: true function without voidIdx dimensions but with error

f_test_loop_importance = function (n, d=8, voidIdx)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  
  X[,voidIdx] = 0
  Ytrue.void <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y.void <- Ytrue.void + rnorm(n, 0, 0.5)
  return(data.frame(Ytrue, Y, Ytrue.void, Y.void))
}

f_test_chain_importance = function(n, d = 8, voidIdx){
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4])^2 + 
    10*(X[,3] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  
  X[,voidIdx] = 0
  Ytrue.void <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4])^2 + 
    10*(X[,3] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y.void <- Ytrue.void + rnorm(n, 0, 0.5)
  return(data.frame(Ytrue, Y, Ytrue.void, Y.void))
}

f_test_complete_importance = function(n, d = 8, voidIdx)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 
    10 * X[,3] * abs(X[,4] - X[,5]) + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  
  X[,voidIdx] = 0
  Ytrue.void <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 
    10 * X[,3] * abs(X[,4] - X[,5]) + 10*X[,6] + 5*X[,7]
  Y.void <- Ytrue.void + rnorm(n, 0, 0.5)
  return(data.frame(Ytrue, Y, Ytrue.void, Y.void))
}

f_test_3way_importance = function(n, d = 8, voidIdx)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 
    10 *cos(pi*X[,3]) * abs(X[,4] - X[,5])+ 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  
  X[,voidIdx] = 0
  Ytrue.void <- 12*sin(pi*X[,1]*X[,2]) + 
    10 *cos(pi*X[,3]) * abs(X[,4] - X[,5])+ 10*X[,6] + 5*X[,7]
  Y.void <- Ytrue.void + rnorm(n, 0, 0.5)
  return(data.frame(Ytrue, Y, Ytrue.void, Y.void))
}

#### --------- Partial ranking distance ----------####
# This function is to calculate the distance between two partial rankings
# It is a generalization of the Kendall distance.
# Inputs: r1: a list of ranked clusters
#         r2: the list of another ranked clusters
#         p: penalty parameter, default 1/2, corresponding to L1 distance 
#            K-profiles.

get_pos = function(r){
  lens = unlist(sapply(r, length))
  len_cum = c(0, cumsum(lens))
  len_cum = len_cum[-length(len_cum)]
  ranking = len_cum + (lens + 1)/2
  ranking = rep(ranking, lens)
  return(cbind(unlist(r),ranking))
}
gKendall = function(r1, r2, p = 1/2){
  idx = unlist(r1)
  if(any(!(unlist(r2) %in% idx))){
    stop('The elements in the two rankings should be identical')
  }
  
  pairs = get_pairIdx(idx)
  r1_pos = get_pos(r1)
  r2_pos = get_pos(r2)
  
  K = rep(NA, nrow(pairs))
  for(k in 1:nrow(pairs)){
    pair = pairs[k,]
    i = pair[1]
    j = pair[2]
    
    r1_i = r1_pos[r1_pos[,1] == i, 2]
    r1_j = r1_pos[r1_pos[,1] == j, 2]
    
    r2_i = r2_pos[r2_pos[,1] == i, 2]
    r2_j = r2_pos[r2_pos[,1] == j, 2]
    
    if(any(c(r1_i == r1_j, r2_i == r2_j))){
      if(all(c(r1_i == r1_j, r2_i == r2_j))){
        # case 1: i and j are in the same bucket
        K[k] = 0
      }else{
        # case 2: i and j are in the same bucket in one list but not in the other
        K[k] = p
      }
    }else{
      # case 3: i and j are in different buckets
      if((r1_i - r1_j) * (r2_i - r2_j) > 0){
        # same direction
        K[k] = 0
      }else{
        # opposite direction
        K[k] = 1
      }
    }
  }
  return(sum(K))
}

# test
r1 = list(c(1,2), 6, 7, 3:5, 8)
## r2 the same as r1
r2 = r1
gKendall(r1, r2)
## r2 are all singletons
r2 = as.list(c(1,2,6,7,3:5,8))
gKendall(r1, cluster$cluster[order(cluster$cluster.f_sens,decreasing = T)])
