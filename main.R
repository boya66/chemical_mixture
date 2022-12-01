library(lhs)
library(hetGP)
library(gtools) # obtain combinations
library(doParallel)
library(bkmr)
library(bartMachine)

cl <- makeCluster(4)
registerDoParallel(cl)
eps = sqrt(.Machine$double.eps)
source("sensitivity_hetGP.R")
source("data_generator.R")
source("var_cluster.R")
source("var_cluster_pair.R")
source("interaction.R")
source("interaction_thres.R")
source("utils.R")
d = 8
K = 50


for(k in 1:K){
  seed = k
  set.seed(seed)
  data <- f_test_chain(300)
  tStart = Sys.time()
  cluster = var_cluster_exhausted(data = data, d = d, N = 3000, J= 20, seed = seed, alpha = 0.05)
  tEnd = Sys.time()
  tElapse1 = tEnd - tStart
  tStart = Sys.time()
  cluster_pair = var_cluster_pair(data = data, d = d, N = 3000, J = 20, seed = seed, alpha = 0.05)
  tEnd = Sys.time()
  tElapse2 = tEnd - tStart
  fitkm = kmbayes(y= data[,9], Z = data[,1:8], varsel = T)
  fitbart = bartMachine(data[,1:8], data[,9])
  save('cluster', 'cluster_pair', 'fitkm', 'fitbart', 'tElapse1','tElapse2',
       file = paste0("../results/chain_s1_",seed,".Rdata"))
}


#### --------------- Simulation summary ---------------- ####
#### Part 1: clustering
d = 8
K = 50
prop = 0.2
ARI <- integration <- acontamination <- matrix(NA, K, 3)
A = list(1:2, 3:5, 6, 7, 8)
for(k in 1:K){
  load(paste0("../results/chain/chain_s1",k,".Rdata"))
  # exhausted approach
  ARI[k,1] = get_ARI(A, cluster$cluster)
  integration[k,1] = get_integration(A, cluster$cluster)
  acontamination[k,1] = get_acontamination(A, cluster$cluster)
  
  # paired approach
  ARI[k,2] = get_ARI(A, cluster_pair$cluster)
  integration[k,2] = get_integration(A, cluster_pair$cluster)
  acontamination[k,2] = get_acontamination(A, cluster_pair$cluster)
  
  # BART
  fitbart = bartMachine(data.frame(cluster$gpi$X0), cluster$gpi$Z0)
  int_sel_bart = interaction_investigator(fitbart)
  pair_info = 
    data.frame(
      get_pairIdx(1:d),
      count.avg = c(t(int_sel_bart$interaction_counts_avg))[c(lower.tri(int_sel_bart$interaction_counts_avg, diag = F))],
      count.sd = c(t(int_sel_bart$interaction_counts_sd))[c(lower.tri(int_sel_bart$interaction_counts_sd, diag = F))])
  
  pair_info = pair_info[order(pair_info$count.avg, decreasing = T),]
  cuts = pair_info$count.avg[1] - c(1,5)*pair_info$count.sd[1]
  cuts = c(sapply(cuts, function(x) which(pair_info$count.avg < x)[1]-1),
           round(d*(d-1)/2 * prop))
  for(j in 1:length(cuts)){
    cut = cuts[j]
    pair = pair_info[1:cut, 1:2]
    pair_idx = unique(c(pair$idx1, pair$idx2))
    if(length(pair_idx) <d){
      cluster_bart = list(setdiff(1:d,unique(c(pair$idx1, pair$idx2))))
    }else{
      cluster_bart = list()
    }
    clus_tmp = list()
    clus_tmp[[1]] = c(pair$idx1[1],pair$idx2[1])
    for(p in 2:nrow(pair)){
      pair_tmp = c(pair$idx1[p], pair$idx2[p])
      if(length(clus_tmp) == 1){
        if(any(pair_tmp %in% clus_tmp[[1]])){
          clus_tmp[[1]] = unique(c(clus_tmp[[1]],pair_tmp))
        }else{
          new_clus = list()
          new_clus[[1]] = pair_tmp
          clus_tmp = c(clus_tmp,new_clus)
        }
      }else{
        if(any(sapply(clus_tmp, function(x) any(pair_tmp %in% x)))){
          clus_num = which(sapply(clus_tmp, function(x) any(pair_tmp %in% x)))
          clus_tmp[[clus_num]] = unique(c(clus_tmp[[clus_num]],pair_tmp))
        }else{
          new_clus = list()
          new_clus[[1]] = pair_tmp
          clus_tmp = c(clus_tmp,new_clus)
        }
      }
    }
    cluster_bart = c(cluster_bart, clus_tmp)
    
  }
}
#### --------------- Application ------------- ####
data = read.csv("../data/working data/TSH.csv")
head(data)
d = ncol(data) - 1

idx_comp = which(!is.na(data[,d+1]))
data = data[idx_comp,]
cluster = var_cluster_thres(data = data, d = 8, N = 3000, J = 20, seed = 518)

# BKMR with variable selection
library(bkmr)
fitkm = kmbayes(y= data[,9], Z = data[,1:8], varsel = T)
summary(fitkm)
ExtractPIPs(fitkm)

# BART with variable selection
fitbart = bartMachine(data[,1:8], data[,9])
var_sel_bart = var_selection_by_permute(fitbart, plot = F)
sum(var_sel_bart$var_true_props_avg)
print(var_sel_bart$important_vars_local_names)
print(var_sel_bart$important_vars_global_max_names)
print(var_sel_bart$important_vars_global_se_names)

var_imp_bar = investigate_var_importance(fitbart)
int_sel_bart = interaction_investigator(fitbart)
int_sel_bart$interaction_counts_avg
int_sel_bart$interaction_counts_sd
