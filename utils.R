get_pairIdx <- function(x){
  n = length(x)
  N = n * (n-1)/2
  idx1 = rep(x[-n], (n-1):1)
  idx2 = x[-1]
  if(n > 2){
    for(i in 2:(n-1)){
      idx2 = c(idx2, x[-(1:i)])
    }
  }
  return(cbind(idx1, idx2))
}

#### ---- Rank the clusters ------ #####
## This function is for identifying important clusters. It ranks the clusters by
## relative importance and cut by explained variance for variable selection.
ExtractCluster <- function(clus, clus_f_sens, clus_sens_diff = NULL,
                           prop.var = 0.90){
  orderIdx = order(clus_f_sens, decreasing = T)
  clus_f_sens = clus_f_sens[orderIdx]
  clus = clus[orderIdx]
  prop.vars = cumsum(clus_f_sens)/sum(clus_f_sens)
  clus_sel = clus[1 : which(prop.vars > prop.var)[1]]
  return(clus_sel)
}

# test
# ExtractCluster(cluster$cluster, cluster$cluster.f_sens, prop.var = 0.95)
# ExtractCluster(cluster_pair$clus, cluster_pair$clus_f_sens, prop.var = 0.95)
# ExtractPIPs(fitkm)

#### ------------ classification measurements ----------------####
# Inputs: true_active: a T/F vector indicating if each element is included in the 
#                      true function.
#         sel_active: a vector of selected dimensions
get_selectionCriteria <- function(true_active, sel_active){
  sel_tf = rep(F, length(true_active))
  sel_tf[sel_active] <- T
  
  conf.m <- table(sel_active, true_active)
  FNR <- conf.m[1,2] / (conf.m[1,1] + conf.m[1,2])
  FPR <- conf.m1[2,1] / (conf.m1[2,1] + conf.m1[2,2])
  sensitivity <- conf.m[2,2] / (conf.m[2,2] + conf.m[1,2])
  specificity <- conf.m[1,1] / (conf.m[1,1] + conf.m[2,1])
  F1 <- conf.m[2,2] / (conf.m[2,2] + 0.5*(conf.m[2,1] + conf.m[1,2]))
  return(data.frame(FNR, FPR, sensitivity, specificity,F1))
}

#### ----------- clustering criteria -------------------####
## Adjusted rand index (ARI), defined on four cases
#   a: # of pairs of points from the same cluster in both A and B
#   b: # of pairs of points from the same cluster in A but different clusters in B
#   c: # of pairs of points from the same cluster in B but the same cluster in A
#   d: # of pairs of points from different clusters in both A and B
# The maximum value is 1. When random assigned with equal probability, the 
# expected value is 0.
##  Inputs: A: a list of clusters of idx, the true cluster
##          B: a list of clusters of idx
get_ARI <- function(A,B){
  idx = unlist(A)
  K = length(idx)
  A_pos = rep(1:length(A), sapply(A, length))
  idx_B = unlist(B)
  B_pos = rep(1:length(B), sapply(B, length))
  if(any(!(idx_B %in% idx))){
    stop("The indices should be identical in two clusterings.")
  }
  # define the counters for the four cases
  a <- b <- c <- d <- 0
  pairs = get_pairIdx(idx)
  for(p in 1:nrow(pairs)){
    i = pairs[p,1]
    j = pairs[p,2]
    
    if(A_pos[idx == i] == A_pos[idx == j]){
      if(B_pos[idx_B == i] == B_pos[idx_B == j]){
        a = a + 1
      }else{
        b = b + 1
      }
    }else{
      if(B_pos[idx_B == i] == B_pos[idx_B == j]){
        c = c + 1
      }else{
        d = d + 1
      }
    }
  }
  
  # ARI
  ARI = ((a+d) * K*(K-1)/2 - ((a+b)*(a+c) + (b+d)*(c+d))) /
    ((K*(K-1)/2)^2 - ((a+b)*(a+c) + (b+d)*(c+d)))
  return(ARI)
}

## Integration, the percentage of data points from given cluster of true partion
# that are in the same cluster in partition B.
## Inputs: A: the list of true partitioning
##         B: the list of partitioning
get_integration <- function(A,B){
  idx_B = unlist(B)
  B_pos = rep(1:length(B), sapply(B, length))
  Intj = sapply(A, 
         function(x){
            B_clus = B_pos[idx_B %in% x]
            B_clus = table(B_clus)
            res = max(B_clus)
            res = res / length(x)
            return(res)
          })
  return(sum(Intj)/length(A))
}

## Acontamination, the percentage of the data in the integrating cluster B are 
# from A. It is a complementary to integration
## Inputs: A: the list of true partitioning
##         B: the list of partitioning
get_acontamination <- function(A,B){
  idx_B = unlist(B)
  B_pos = rep(1:length(B), sapply(B, length))
  Acontj = sapply(A, 
                function(x){
                  B_clus = B_pos[idx_B %in% x]
                  B_clus = table(B_clus)
                  res = max(B_clus)
                  res.name = as.numeric(names(B_clus[which.max(B_clus)]))
                  res = res / length(B[[res.name]])
                  return(res)
                })
  return(sum(Acontj)/length(A))
}

