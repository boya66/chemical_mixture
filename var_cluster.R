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
var_cluster_exhausted <- function(data, d, N = 3000, J= 20, seed = 123, 
                              alpha = 0.05, uk = NULL, saveData = F, g_max = NULL,...){
    Y = data[, d+1]
    upper = 30
    if(!is.null(g_max)){
      gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), 
                     upper = rep(upper, d), covtype = 'Gaussian', maxit = 500,
                     noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), g_max)))
      while(gpi$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
        upper = upper + 1
        gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), 
                       upper = rep(upper, d), covtype = 'Gaussian', maxit = 500,
                       noiseControl = list(g_bounds = c(sqrt(.Machine$double.eps), g_max)))
      }
    }else{
      gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), 
                     upper = rep(upper, d), covtype = 'Gaussian', maxit = 500)
      while(gpi$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
        upper = upper + 1
        gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), 
                       upper = rep(upper, d), covtype = 'Gaussian', maxit = 500)
      }
    }
    
    
    int_thres = interaction_thres_mc(data[,1:d], sqrt(gpi$nu_hat * gpi$g), N = 3000, 
                                     J, d, seed = seed, uk = uk) #85s/4cores
    #upper;gpi_perm
    # image(seq(0,1,len=50), seq(0,1,len=50), matrix(y.pred$mean,ncol = 50))
    # Single variable importance ordering
    # repeat the calculation for J times and take the median as the value
    sens = interaction_mc(N , J = J, d, gpi, seed = seed, uk = uk) #53s/4cores
    # calculate the scores as the medians
    f_sens = apply(matrix(sens[,1], nrow = d), 1, median)
    t_sens = apply(matrix(sens[,2], nrow = d), 1, median)
    sens_diff = apply(matrix(sens[,3], nrow = d), 1, median)
    
    sens_diff_thres = min(quantile(int_thres, 1-alpha),1e-3)
  
    # sort the total sensitivity
    sens_idx = 1:d
    idx_single = sens_idx[sens_diff<=sens_diff_thres]
    idx_interaction = setdiff(sens_idx, idx_single)
    while(length(idx_interaction) <= 1 & sens_diff_thres > 1e-8){
      sens_diff_thres = sens_diff_thres/10
      idx_single = sens_idx[sens_diff<=sens_diff_thres]
      idx_interaction = setdiff(sens_idx, idx_single)
    }
    if(length(idx_interaction) == 1){
      idx_interaction = sens_idx[sens_diff > 0]
      idx_single = setdiff(sens_idx, idx_interaction)
      if(length(idx_interaction == 1)){
        idx_single = sens_idx
        idx_interaction = setdiff(sens_idx, idx_single)
      }
    }
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
      comb_f_sens <- comb_t_sens <- comb_diff <- rep(NA, d_pool)
      for(p in 1:d_pool){
        comb_tmp = interaction_mult_mc(N, J, d, gpi, comb[p,], seed + p, uk = uk)
        comb_f_sens[p] = median(comb_tmp[,1])
        comb_t_sens[p] = median(comb_tmp[,2])
        comb_diff[p] = median(comb_tmp[,3])
      }
      
      sel_new = comb[which.min(comb_diff),ncol(comb)]
      comb_f_sens = comb_f_sens[which.min(comb_diff)]
      comb_t_sens = comb_t_sens[which.min(comb_diff)]
      comb_diff = min(comb_diff)
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
                clus_f_sens[clus_num] = f_sens[idx_sel]
                clus_t_sens[clus_num] = t_sens[idx_sel]
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
          clus_f_sens[clus_num] = f_sens[idx_sel]
          clus_t_sens[clus_num] = t_sens[idx_sel]
        }
      }
    }
    
    # Add the pairwise info
    pairInfo = vector("list", length=length(clus))
    for(clus_num in 1: length(clus)){
      if(length(clus[[clus_num]]) > 1){
        pair_idx = get_pairIdx(clus[[clus_num]])
        pair_info = matrix(NA, nrow(pair_idx), 7)
        pair_info[,1:2] = pair_idx
        
        for(i in 1: nrow(pair_idx)){
          idx = pair_idx[i,]
          #interaction_mult(N,d,gpi,idx,seed = seed)
          sens = interaction_mult_mc(N,J,d,gpi,idx,seed = seed)
          pair_info[i, 3] = median(sens[,1])
          pair_info[i, 4] = median(sens[,2])
          pair_info[i, 5] = median(sens[,3])
        }
        
        for(i in 1:nrow(pair_idx)){
          idx = pair_idx[i,]
          pair_info[i, 6] =  max(sens_diff[c(idx[1], idx[2])]) - 
            pair_info[i, 5]
          pair_info[i, 7] = - f_sens[idx[1]] - f_sens[idx[2]] +
            pair_info[i, 3]
        }
        pairInfo[[clus_num]] = pair_info
      }
    }
    return(list(cluster = clus, 
                cluster.f_sens = clus_f_sens, 
                cluster.t_sens = clus_t_sens,
                pairInfo = pairInfo,
                singleInfo = data.frame(
                  f_sens = f_sens,
                  t_sens = t_sens,
                  I = sens_diff
                ),
                sens_diff_thres = sens_diff_thres,
                gpi = gpi, 
                settings = data.frame(N =N, J = J, d=d, upper = upper)))
}
  