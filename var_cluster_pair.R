###############################################################################
## The pairwise procedure of to obtain variable clusters
## Variables in each cluster has some pairwise interaction effects on the outcome
## Input: 
##   data: an observation * (d+1) matrix/data.frame, outcome in the last col
##   d: number of predictors
##   N: number of grid points for evaluating GP model
##   J: maximum number of variables in one cluster, default NULL
##   alpha: the significance level of identified clusters
## Output: a list of
##         pair_diff: a matrix of pairwise sensitivity difference
##         clus: a list of variables, each element contains vars in one cluster
##         clus_f_sens: the first-order sensitivity of all clusters
##         clus_t_sens: the total sensitivity of all clusters
##         pair_f_sens: a matrix of all pairwise first-order sensitivity
##         pair_t_sens: a matrix of all pairwise second_order sensitivity
###############################################################################
var_cluster_pair <- function(data, d, N = 3000, J= 20, seed = 123, 
                              alpha = 0.05,  ...){
  Y = data[, d+1]
  upper = 30
  gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), upper = rep(upper, d),
                 covtype = 'Gaussian', maxit = 500)
  while(gpi$msg != "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"){
    upper = upper + 1
    gpi = mleHomGP(X = as.matrix(data[,1:d]), Z = Y, lower = rep(1, d), upper = rep(upper, d),
                   covtype = 'Gaussian', maxit = 500)
  }
  
  int_thres = interaction_thres_mc(data[,1:d], sqrt(gpi$nu_hat * gpi$g), N, J,
                                   d, seed = seed) #85s/4cores
  #upper;gpi_perm
  # image(seq(0,1,len=50), seq(0,1,len=50), matrix(y.pred$mean,ncol = 50))
  # Single variable importance ordering
  # repeat the calculation for J times and take the median as the value
  sens = interaction_mc(N, J, d, gpi, seed = seed) #53s/4cores
  # calculate the scores as the medians
  f_sens = apply(matrix(sens[,1], nrow = d), 1, median)
  t_sens = apply(matrix(sens[,2], nrow = d), 1, median)
  sens_diff = apply(matrix(sens[,3], nrow = d), 1, median)
  
  sens_diff_thres = interaction_diff_thres(int_thres, B = 1000, alpha = alpha)
  
  sens_idx = 1:d
  
  # setup an NA matrix to store pairwise information
  pair_idx = get_pairIdx(sens_idx)
  pair_info = matrix(NA, nrow(pair_idx), 7)
  pair_info[,1:2] = pair_idx
  
  for(i in 1: nrow(pair_idx)){
    idx = pair_idx[i,]
    #interaction_mult(N,d,gpi,idx,seed = seed)
    sens = interaction_mult_mc(N,J,d,gpi,idx,seed = seed)
    pair_info[i, 3] = median(sens[,1])
    pair_info[i, 4] = median(sens[,2])
    pair_info[i, 5] = median(sens[,3])
    pair_info[i, 6] =  max(sens_diff[c(idx[1], idx[2])]) - 
      pair_info[i, 5]
    pair_info[i, 7] = - f_sens[idx[1]] - f_sens[idx[2]] +
      pair_info[i, 3]
  }
  colnames(pair_info) <- c("idx1","idx2","f_sens","t_sens","sens_diff","delta",
                           "inter_f_sens")
  
  pair_sel = pair_info[pair_info[,"delta"] > sens_diff_thres, 1:2]

  # variable indices with interaction
  idx_interaction = unique(c(pair_sel))
  
  # variable indices without interaction
  idx_single = setdiff(sens_idx, idx_interaction)
  
  # put single variables into clusters
  clus = as.list(idx_single)
  clus_f_sens = f_sens[idx_single]
  clus_t_sens = t_sens[idx_single]
  clus_sens_diff = sens_diff[idx_single]
  
  npair = 0
  # setup a list to store pairwise information of the final clusters
  if(length(idx_interaction) == 0){
    clus_pairInfo = pair_info
  }else{
    clus_pairInfo = clus
    clus_pairInfo[1:length(idx_single)] = NA
    clus_num = length(clus) + 1
    if(length(idx_interaction) == 2){
      clus_tmp = idx_interaction
      pair_info_tmp = pair_info[pair_info[,"delta"] > sens_diff_thres,,drop=F ]
      clus[[clus_num]] = clus_tmp
      clus_pairInfo[[clus_num]] = pair_info_tmp
      clus_f_sens[clus_num] = pair_info_tmp[, 3]
      clus_t_sens[clus_num] = pair_info_tmp[, 4]
      clus_sens_diff[clus_num] = pair_info_tmp[, 5]
    }else{
      clus_tmp = pair_sel[1,]
      pair_sel = pair_sel[-1,,drop = F]
      npair = nrow(pair_sel)
    }
  } 
  
  # test
  # clus = list()
  # clus_num = 1
  # clus_tmp = c(1,2)
  # pair_sel = matrix(c(2,4,3,3,5,6), ncol = 2)
  # npair = nrow(pair_sel)
  
  while(npair > 0){
    idx = NULL
    for(i in 1:npair){
      if(any(pair_sel[i,] %in% clus_tmp)){
        clus_tmp = c(clus_tmp, pair_sel[i,])
        clus_tmp = unique(clus_tmp)
        idx = c(idx, i)
      }
    }
    clus[[clus_num]] = clus_tmp
    if(length(clus_tmp) == 2){
      pair_info_tmp = pair_info[pair_info[,"idx1"] %in% clus_tmp &
                                  pair_info[,"idx2"] %in% clus_tmp,, drop = F]
      clus_pairInfo[[clus_num]] = pair_info_tmp
      clus_f_sens[clus_num] = pair_info_tmp[, 3]
      clus_t_sens[clus_num] = pair_info_tmp[, 4]
      clus_sens_diff[clus_num] = pair_info_tmp[, 5]
      clus_tmp = pair_sel[1, ]
      clus_num = clus_num + 1
      pair_sel = pair_sel[-1,,drop = F]
      npair = nrow(pair_sel)
    }else{
      sens_tmp = interaction_mult_mc(N, J, d, gpi, clus_tmp, seed)
      clus_f_sens[clus_num] = median(sens_tmp[,1])
      clus_t_sens[clus_num] = median(sens_tmp[,2])
      clus_sens_diff[clus_num] = median(sens_tmp[,3])
      clus_pairInfo[[clus_num]] = pair_info[pair_info[,"idx1"] %in% clus_tmp &
                                              pair_info[,"idx2"] %in% clus_tmp,]
      clus_num = clus_num + 1
      pair_sel <- pair_sel[-idx,,drop = F]
      if(nrow(pair_sel) == 0){
        break
      }else{
        clus_tmp = pair_sel[1, ]
        pair_sel <- pair_sel[-1,,drop=F]
        npair = nrow(pair_sel)
      }
    }
    if(npair == 0){
      clus[[clus_num]] = clus_tmp
      pair_info_tmp = pair_info[pair_info[,"idx1"] %in% clus_tmp &
                                  pair_info[,"idx2"] %in% clus_tmp,,drop = F]
      clus_pairInfo[[clus_num]] = pair_info_tmp
      clus_f_sens[clus_num] = pair_info_tmp[, 3]
      clus_t_sens[clus_num] = pair_info_tmp[, 4]
      clus_sens_diff[clus_num] = pair_info_tmp[, 5]
    }
  }
  return(list(clus = clus,
              sens_diff_thres = sens_diff_thres,
              clus_pairInfo = clus_pairInfo,
              clus_t_sens = clus_t_sens,
              clus_f_sens = clus_f_sens,
              clus_sens_diff = clus_sens_diff, 
              pair_info = pair_info,
              gpi = gpi,
              settings = data.frame(N =N, J = J, d=d, upper = upper)))
}
