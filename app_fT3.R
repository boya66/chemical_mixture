library(lhs)
library("hetGP", lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.1")
library(gtools) # obtain combinations
library(doParallel)
# library(bkmr, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.1")
#library(bartMachine)

cl <- makeCluster(14)
registerDoParallel(cl)
eps = sqrt(.Machine$double.eps)
source("sensitivity_hetGP.R")
source("var_cluster.R")
source("var_cluster_thres_pair.R")
source("interaction.R")
source("interaction_thres.R")
source("utils.R")

# data = read.csv("~/Project9_chemical_sens/data/working data/data_rescale_v0.csv")
visit = 0
biomarker = "fT3"
data = read.csv("data/data_rescale_v",visit,".csv")

d = 10
N = 9000

g_max = sd(data[,biomarker])

PFAS = names(data)[1:8]
cov_var = c("Age","PreBMI", "ENRRACE_fm002", "INFSEXR1_fm024", "Education", "nulliparity",
            "GestationalAge", "Cotinine","TotalLipid")
data_ahei = data[,c(PFAS,"AHEI","GDM",biomarker,cov_var)]
data_dash = data[,c(PFAS,"DASH","GDM",biomarker,cov_var)]
data_amed = data[,c(PFAS,"aMED","GDM",biomarker,cov_var)]

uk = list(dist = c(rep("normal",d-1),"bernoulli"),
          param = cbind(c(colMeans(data_ahei[,1:(d-1)]),49/50),
                        c(apply(data_ahei[,1:(d-1)],2,sd),NA)))
cluster_ahei = var_cluster_exhausted(data = data_ahei, d = d, N = N, J = 25, 
                                     g_max = g_max,seed = 518, uk = uk)
# cluster_pair_ahei = var_cluster_thres_pair(data = data_ahei, d = d, N = N, J = 25, 
#                                            g_max = g_max, seed = 518, uk = uk)
save("cluster_ahei",file = paste0("data/v0/", biomarker, "_ahei_cov.Rdata"))

uk = list(dist = c(rep("normal",d-1),"bernoulli"),
          param = cbind(c(colMeans(data_dash[,1:(d-1)]),49/50),
                        c(apply(data_dash[,1:(d-1)],2,sd),NA)))
cluster_dash = var_cluster_exhausted(data = data_dash, d = d, N = N, J = 25, 
                                     g_max = g_max,seed = 518, uk = uk)
# cluster_pair_dash = var_cluster_thres_pair(data = data_dash, d = d, N = N, J = 25, 
#                                            g_max = g_max, seed = 518)
save("cluster_dash", file = paste0("data/v0/", biomarker, "_imp_dash_cov.Rdata"))

uk = list(dist = c(rep("normal",d-1),"bernoulli"),
          param = cbind(c(colMeans(data_amed[,1:(d-1)]),49/50),
                        c(apply(data_amed[,1:(d-1)],2,sd),NA)))
cluster_amed = var_cluster_exhausted(data = data_amed, d = d, N = N, J = 25, 
                                     g_max = g_max,seed = 518, uk =uk)
# cluster_pair_amed = var_cluster_thres_pair(data = data_amed, d = d, N = N, J = 25, 
#                                            g_max = g_max, seed = 518)
save("cluster_amed", file = paste0("data/v0/",biomarker,"_imp_amed_cov.Rdata"))

