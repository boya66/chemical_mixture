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
var_cluster_thres_pair <- function(data, d, N = 3000, J= 20, seed = 123, 
                             alpha = 0.05, uk = NULL, g_max=NULL,...){
    # The response variable is following the first d columns of chemicals.
    Y = data[, d+1]
    # The initial value of upper limit for optimization
    upper = 30
    
    # fit the GP model: 
    #   X the first d columns of X, 
    #   Z  = Y, the response, the d + 1 column of data
    #   lower = 1, the lower bound for optimization
    #   upper = upper, the initial upper bound for optimization
    #   covtype = "Gaussian", the gaussian kernel will be used
    #   maxit = 500, the maximum iterations, set to be 500.
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
    
    #upper;gpi_perm
    # image(seq(0,1,len=50), seq(0,1,len=50), matrix(y.pred$mean,ncol = 50))
    # Single variable importance ordering
    # repeat the calculation for J times and take the median as the value
    sens = interaction_mc(N, J, d, gpi, seed = seed, uk = uk) #53s/4cores
    # calculate the scores as the medians
    f_sens = apply(matrix(sens[,1], nrow = d), 1, median)
    t_sens = apply(matrix(sens[,2], nrow = d), 1, median)
    sens_diff = apply(matrix(sens[,3], nrow = d), 1, median)
    
    int_thres = interaction_thres_mc(data[,1:d], sqrt(gpi$nu_hat * gpi$g), N = 3000, 
                                     J, d, seed = seed, uk = uk) #85s/4cores
    
    sens_diff_thres = quantile(int_thres, 1-alpha)
    sens_diff_thres = min(1e-3, sens_diff_thres)
    # IQR method for extreme large value outliers
    # sens_diff_thres = 2.5*quantile(int_thres,0.75) -
    #     1.5*quantile(int_thres,0.25)
    
    sens_idx = 1:d
    idx_single = sens_idx[sens_diff<=sens_diff_thres]
    idx_interaction = setdiff(sens_idx, idx_single)
    
    # setup an NA matrix to store pairwise information
    if(length(idx_interaction) == 0){
        pair_idx = get_pairIdx(sens_idx)
    }else{
        pair_idx = get_pairIdx(idx_interaction,d = d)
    }
    
    pair_info = matrix(NA, nrow(pair_idx), 7)
    pair_info[,1:2] = pair_idx
    
    for(i in 1: nrow(pair_idx)){
        idx = pair_idx[i,]
        #interaction_mult(N,d,gpi,idx,seed = seed)
        sens = interaction_mult_mc(N,J,d,gpi,idx,seed = seed, uk = uk)
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
    
    pair_sel = pair_info[pair_info[,"delta"] > 0, 1:2, drop = F]
    
    # variable indices with interaction
    idx_interaction = unique(c(pair_sel))
    
    # variable indices without interaction
    idx_single = setdiff(sens_idx, idx_interaction)
    
    # put single variables into clusters
    clus = as.list(idx_single)
    clus_f_sens = f_sens[idx_single]
    clus_t_sens = t_sens[idx_single]
    clus_sens_diff = sens_diff[idx_single]
    
    # setup a list to store pairwise information of the final clusters
    if(length(idx_interaction) == 0){
        clus_pairInfo = pair_info
    }else{
        ss = length(clus)
        clus_pairInfo = list()
        clus_pairInfo[1:length(idx_single)] = NA
        
        clus_tmp <- pairInfo_tmp <- list()
        clus_tmp[[1]] = c(pair_sel[1,1],pair_sel[1,2])
        pairInfo_tmp[[1]] = 
            pair_info[pair_info[,1] == pair_sel[1,1] & 
                          pair_info[,2] == pair_sel[1,2], ,drop = F]
        
        if(nrow(pair_sel) >= 2){
            for(p in 2:nrow(pair_sel)){
                pair_tmp = c(pair_sel[p,1], pair_sel[p,2])
                if(length(clus_tmp) == 1){
                    if(any(pair_tmp %in% clus_tmp[[1]])){
                        clus_tmp[[1]] = unique(c(clus_tmp[[1]],pair_tmp))
                        pairInfo_tmp[[1]] = 
                            rbind(pairInfo_tmp[[1]],
                                  pair_info[pair_info[,1] == pair_tmp[1] & 
                                                pair_info[,2] == pair_tmp[2],, drop = F])
                    }else{
                        new_clus <- new_pairInfo <- list()
                        new_clus[[1]] = pair_tmp
                        new_pairInfo[[1]] = pair_info[pair_info[,1] == pair_tmp[1] & 
                                                          pair_info[,2] == pair_tmp[2],, drop = F]
                        clus_tmp = c(clus_tmp,new_clus)
                        pairInfo_tmp = c(pairInfo_tmp, new_pairInfo)
                    }
                }else{
                    if(any(sapply(clus_tmp, function(x) any(pair_tmp %in% x)))){
                        clus_num = which(sapply(clus_tmp, function(x) any(pair_tmp %in% x)))
                        if(length(clus_num) == 2){
                            clus_tmp[[clus_num[1]]] = unique(c(clus_tmp[[clus_num[1]]],
                                                               clus_tmp[[clus_num[2]]],
                                                               pair_tmp))
                            clus_tmp = clus_tmp[-clus_num[2]]
                            pairInfo_tmp[[clus_num[1]]] = 
                                rbind(pairInfo_tmp[[clus_num[1]]],pairInfo_tmp[[clus_num[2]]],
                                      pair_info[pair_info[,1] == pair_tmp[1] & 
                                                    pair_info[,2] == pair_tmp[2],, drop = F])
                            pairInfo_tmp = pairInfo_tmp[-clus_num[2]]
                        }else{
                            clus_tmp[[clus_num]] = unique(c(clus_tmp[[clus_num]],pair_tmp))
                            pairInfo_tmp[[clus_num]] = 
                                rbind(pairInfo_tmp[[clus_num]],
                                      pair_info[pair_info[,1] == pair_tmp[1] & 
                                                    pair_info[,2] == pair_tmp[2],, drop = F])
                        }
                    }else{
                        new_clus <- new_pairInfo <- list()
                        new_clus[[1]] = pair_tmp
                        new_pairInfo[[1]] = pair_info[pair_info[,1] == pair_tmp[1] & 
                                                          pair_info[,2] == pair_tmp[2],, drop = F]
                        clus_tmp = c(clus_tmp,new_clus)
                        pairInfo_tmp = c(pairInfo_tmp, new_pairInfo)
                    }
                }
            }
        }
        clus = c(clus, clus_tmp)
        clus_pairInfo = c(clus_pairInfo, pairInfo_tmp)
        
        for(s in (ss+1):length(clus)){
            clus_tmp = clus[[s]]
            sens_tmp = interaction_mult_mc(N, J, d, gpi, clus_tmp, seed, uk = uk)
            clus_f_sens[s] = median(sens_tmp[,1])
            clus_t_sens[s] = median(sens_tmp[,2])
            clus_sens_diff[s] = median(sens_tmp[,3])
        }
    } 
    
    return(list(cluster = clus,
                sens_diff_thres = sens_diff_thres,
                clus_pairInfo = clus_pairInfo,
                clus_t_sens = clus_t_sens,
                clus_f_sens = clus_f_sens,
                clus_sens_diff = clus_sens_diff, 
                pair_info = pair_info,
                gpi = gpi,
                settings = data.frame(N =N, J = J, d=d, upper = upper)))
}
