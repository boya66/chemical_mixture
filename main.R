library(lhs)
library(laGP)
library(hetGP)
library(gtools) # obtain combinations
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
eps = sqrt(.Machine$double.eps)
source("sensitivity_hetGP.R")
source("var_cluster.R")
source("interaction.R")
d = 8
K = 50


tStart = Sys.time()
for(k in 1:K){
  seed = k
  set.seed(seed)
  data <- f_test(300)
  cluster = var_cluster_thres(data = data, d = d, N = 3000, J= 20, seed = seed)
  save('cluster', file = paste0("../results/cluster",seed,".Rdata"))
}
tEnd = Sys.time()
tUse = difftime(tEnd, tStart, units = "mins")/K
save('tUse','clus.test','N','J',file = "../results/clus_N3000J20_max_thres_50_all8.Rdata")

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
var_sel_bart = var_selection_by_permute(fitbart)
sum(var_sel_bart$var_true_props_avg)
print(var_sel_bart$important_vars_local_names)
print(var_sel_bart$important_vars_global_max_names)
print(var_sel_bart$important_vars_global_se_names)

var_imp_bar = investigate_var_importance(fitbart)
int_sel_bart = interaction_investigator(fitbart)
int_sel_bart$interaction_counts_avg
int_sel_bart$interaction_counts_sd
