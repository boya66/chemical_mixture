###############################################################################
## The exhausted procedure of to obtain variable clusters
## Variables in each cluster has interaction effects on the outcome
## Input: 
##   data: an observation * (d+1) matrix/data.frame, outcome in the last col
##   d: number of predictors
##   N: number of grid points for evaluating GP model
##   J: maximum number of variables in one cluster, default NULL
##   var_perc: the cut-off for cutting out unimportant variables in single
##                  selection, default NULL, no cut-off
## Output: a list of
##         var: a list of variables ordered by the f_sens of the leading var
##         sens: a list of all sensitivity indices for variables in each cluster
###############################################################################
var_cluster_thres <- function(data, d, N = 3000, J= 20, seed = 123, 
                              alpha = 0.05, saveData = F, ...){
    Y = data[, d+1]
    upper = 30
    gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), upper = rep(upper, d),
                covtype = 'Gaussian', maxit = 500)
    while(gpi$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
      upper = upper + 1
      gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), upper = rep(upper, d),
                     covtype = 'Gaussian', maxit = 500)
    }
    
    # The measurement error is tau square times nugget g from Section 5.2.2 on Book "Surrogate" by Grammacy
    # All variables have no interaction and the first column is null
    Y_noInter = rowSums(data[,1:d]) + rnorm(length(Y),0,sd = sqrt(gpi$nu_hat * gpi$g))
    if(saveData){
      save('data','Y_noInter',file = paste0("../simulation_data/data",seed,".Rdata"))
    }
    upper = 30
    gpi_perm = mleHomGP(X = as.matrix(data[,1:d]), Z = Y_noInter, lower = rep(1, d), upper = rep(upper, d),
                        covtype = 'Gaussian', maxit = 500)
    while(gpi_perm$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
      upper = upper + 1
      gpi_perm = mleHomGP(X = as.matrix(data[,1:d]), Z = Y_noInter, lower = rep(1, d), upper = rep(upper, d),
                       covtype = 'Gaussian', maxit = 500)
    }
    #upper;gpi_perm
    # image(seq(0,1,len=50), seq(0,1,len=50), matrix(y.pred$mean,ncol = 50))
    # Single variable importance ordering
    # repeat the calculation for J times and take the median as the value
    system.time({sens = interaction_mc(N, J, d, gpi)})
    #### -------------------- stopped here -----------------####
    # calculate the scores as the medians
    f_sens = apply(f_sens.rep, 2, median)
    t_sens = apply(t_sens.rep, 2, median)
    sens_diff = apply(sens_diff.rep, 2, median)
    
    sens_diff_thres = quantile(sens_diff_perm.rep, 1-alpha)
  
    # sort the total sensitivity
    sens_idx = 1:d
    idx_single = sens_idx[sens_diff<=sens_diff_thres]
    idx_interaction = setdiff(sens_idx, idx_single)
    clus = as.list(idx_single)
    clus_f_sens = f_sens[idx_single]
    clus_t_sens = t_sens[idx_single]
    sens_idx_ord = idx_interaction[order(sens_diff[idx_interaction], decreasing = T)]
    
    clus_num = length(clus) + 1
    idx_sel = sens_idx_ord[1]
    sel_diff = sens_diff[idx_sel]
    
    d_pool = length(idx_interaction) - 1
    while(d_pool > 0){
      comb = cbind(matrix(rep(idx_sel,each = d_pool),ncol = length(idx_sel)), idx_interaction[which(!(idx_interaction %in% idx_sel))]) 
      comb_f_sens.rep <- comb_t_sens.rep <- comb_diff.rep <- matrix(NA, J, d_pool)
      for(j in 1:J){
        set.seed(seed + j)
        M = lhs::randomLHS(N, d)
        Mprime = lhs::randomLHS(N, d)
        comb_f_sens.rep[j,] <- apply(comb,1,function(j) first_order_sens_mult(gpi,d,M,Mprime,j))
        comb_t_sens.rep[j,] <- apply(comb,1,function(j) total_sens_mult(gpi,d,M,Mprime,j))
        I = comb_t_sens.rep[j,] - comb_f_sens.rep[j,]
        comb_diff.rep[j,] <- ifelse(I < 0, 0, I)
      }
      comb_f_sens = apply(comb_f_sens.rep, 2, median)
      comb_t_sens = apply(comb_t_sens.rep, 2, median)
      #boxplot(comb_diff.rep)
      comb_diff = apply(comb_diff.rep, 2, quantile, 0.5)
      
      sel_new = comb[which.min(comb_diff),ncol(comb)]
      comb_diff = min(comb_diff)
      comb_f_sens = comb_f_sens[which.min(comb_diff)]
      comb_t_sens = comb_t_sens[which.min(comb_diff)]
      if(comb_diff < max(c(sens_diff[sel_new], sel_diff))){
      # if(comb_diff < sel_diff){
        idx_sel = c(idx_sel,sel_new)
        d_pool = d_pool - 1
        sel_diff = comb_diff
        if(d_pool == 0){
          clus[[clus_num]] = idx_sel
          clus_f_sens[clus_num] = comb_f_sens
          clus_t_sens[clus_num] = comb_t_sens
          }else{
            if(comb_diff < sens_diff_thres){
              clus[[clus_num]] = idx_sel
              clus_f_sens[clus_num] = comb_f_sens
              clus_t_sens[clus_num] = comb_t_sens
              clus_num = clus_num + 1
              
              idx_interaction = idx_interaction[!(idx_interaction %in% idx_sel)]
              sens_idx_ord = sens_idx_ord[!(sens_idx_ord %in% idx_sel)]
              
              idx_sel = sens_idx_ord[1]
              d_pool = length(idx_interaction) - 1
              sel_diff = sens_diff[idx_sel]
              if(d_pool == 0){
                clus[[clus_num]] = idx_sel
                clus_f_sens[clus_num] = comb_f_sens
                clus_t_sens[clus_num] = comb_t_sens
              }
              
            }
        }
      }else{
        clus[[clus_num]] = idx_sel
        clus_f_sens[clus_num] = comb_f_sens
        clus_t_sens[clus_num] = comb_t_sens
        clus_num = clus_num + 1
        
        idx_interaction = idx_interaction[!(idx_interaction %in% idx_sel)]
        sens_idx_ord = sens_idx_ord[!(sens_idx_ord %in% idx_sel)]
        d_pool = length(idx_interaction) - 1
        
        idx_sel = sens_idx_ord[1]
        sel_diff = sens_diff[idx_sel]
        if(d_pool == 0){
          clus[[clus_num]] = idx_sel
          clus_f_sens[clus_num] = comb_f_sens
          clus_t_sens[clus_num] = comb_t_sens
        }
      }
    }
    return(list(cluster = clus, cluster.f_sens = clus_f_sens, cluster.t_sens = clus_t_sens,
                gpi = gpi, settings = data.frame(N =N, J = J, d=d, upper = upper), data = data))
}
  