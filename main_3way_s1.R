library(lhs)
library("hetGP", lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.1")
library(gtools) # obtain combinations
library(doParallel)
library(bkmr, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.1")
#library(bartMachine)

cl <- makeCluster(16)
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

for(k in 5:K){
  seed = k
  set.seed(seed)
  data <- f_test_3way(300)
  tStart = Sys.time()
  cluster = var_cluster_exhausted(data = data, d = d, N = 3000, J= 20, seed = seed, alpha = 0.05)
  tEnd = Sys.time()
  tElapse1 = tEnd - tStart
  tStart = Sys.time()
  cluster_pair = var_cluster_pair(data = data, d = d, N = 3000, J = 20, seed = seed, alpha = 0.05)
  tEnd = Sys.time()
  tElapse2 = tEnd - tStart
  fitkm = kmbayes(y= data[,9], Z = data[,1:8], varsel = T)
#  fitbart = bartMachine(data[,1:8], data[,9])
  save('cluster', 'cluster_pair', 'fitkm', 'tElapse1','tElapse2',
       file = paste0("results/order3_s1_",seed,".Rdata"))
}