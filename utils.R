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
  sel_tf[sel_active] <- 
}
