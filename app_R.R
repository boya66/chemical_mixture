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
source("data_generator.R")
source("var_cluster.R")
source("var_cluster_thres_pair.R")
source("interaction.R")
source("interaction_thres.R")
source("utils.R")

# data = read.csv("~/Project9_chemical_sens/data/working data/v0/R_fT3_fT4_imp.csv")
data = read.csv("data/v0/R_fT3_fT4_imp.csv")
data_ahei = data[,-c(10,11)]
data_dash = data[,-c(9,11)]
data_amed = data[,-c(9,10)]

d = ncol(data) - 3
N = 9000
g_max = sd(data[,ncol(data)])
uk = list(dist = c(rep("normal",d-1),"bernoulli"),
          param = cbind(c(colMeans(data_ahei[,1:(d-1)]),49/50),
                        c(apply(data_ahei[,1:(d-1)],2,sd),NA)))
cluster_ahei = var_cluster_exhausted(data = data_ahei, d = d, N = N, J = 25, 
                                     g_max = g_max,seed = 518, uk = uk)
# cluster_pair_ahei = var_cluster_thres_pair(data = data_ahei, d = d, N = N, J = 25, 
#                                            g_max = g_max, seed = 518, uk = uk)
save("cluster_ahei",file = "data/v0/R_fT3_fT4_imp_ahei.Rdata")

uk = list(dist = c(rep("normal",d-1),"bernoulli"),
          param = cbind(c(colMeans(data_dash[,1:(d-1)]),49/50),
                        c(apply(data_dash[,1:(d-1)],2,sd),NA)))
cluster_dash = var_cluster_exhausted(data = data_dash, d = d, N = N, J = 25, 
                                     g_max = g_max,seed = 518, uk = uk)
# cluster_pair_dash = var_cluster_thres_pair(data = data_dash, d = d, N = N, J = 25, 
#                                            g_max = g_max, seed = 518)
save("cluster_dash", file = "data/v0/R_fT3_fT4_imp_dash.Rdata")

uk = list(dist = c(rep("normal",d-1),"bernoulli"),
          param = cbind(c(colMeans(data_amed[,1:(d-1)]),49/50),
                        c(apply(data_amed[,1:(d-1)],2,sd),NA)))
cluster_amed = var_cluster_exhausted(data = data_amed, d = d, N = N, J = 25, 
                                     g_max = g_max,seed = 518, uk =uk)
# cluster_pair_amed = var_cluster_thres_pair(data = data_amed, d = d, N = N, J = 25, 
#                                            g_max = g_max, seed = 518)
save("cluster_amed", file = "data/v0/R_fT3_fT4_imp_amed.Rdata")

