library(hetGP)
library(lhs)
library(gtools) # obtain combinations
eps = sqrt(.Machine$double.eps)
##########################################################################
# single variable version
# return the 1st order/total sensitivity for all variables respectively (a d-dim vector)
# sec 8.2.3 of https://bookdown.org/rbg/surrogates/chap8.html#chap8sens
##########################################################################

# first order sensitivity
first_order_sens = function(gpi, d, M, Mprime){
  N <- nrow(M)
  pM <- predict(gpi,M)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  S <- EE2j <- rep(NA, d)
  # for-loop is faster than sapply
  for(j in 1:d) {
    Mjprime <- Mprime
    Mjprime[,j] <- M[,j]
    pMprime <- predict(gpi, Mjprime)
    EE2j[j] <- (t(pM$mean) %*% pMprime$mean)/(N - 1)
    S[j] <- (EE2j[j] - Ey^2)/Vary
  }
  # since numerator and denominator are all variances in the definition and they
  # should be greater than 0.
  S[S<0] <- 0
  return(S)
}


total_sens = function(gpi, d, M, Mprime){
  N <- nrow(M)
  pM <- predict(gpi, M)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  T_sens <- EE2mj <- rep(NA, d)
  for(j in 1:d) {
    Mj <- M
    Mj[,j] <- Mprime[,j]
    pMj <- predict(gpi, Mj)
    EE2mj[j] <- (t(pM$mean) %*% pMj$mean)/(N - 1)
    T_sens[j] <- 1 - (EE2mj[j] - Ey^2)/Vary
  }
  # since numerator and denominator are all variances in the definition and they
  # should be greater than 0.
  T_sens[T_sens<0] <- 0
  return(T_sens)
}

# true function (the X8 is useless)
f_test = function (n, d=8)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

f_test_no3 = function(n,d = 8){
  X <- randomLHS(n, d)
  Ytrue <- 15*sin(pi*X[,1]*X[,2]) + 20*(X[,3] + X[,4])^2 + 
    20*(X[,3] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}
# # generate data
#   data <- f_test(300)
#   d = 8
# 
# 
# 
# ### test
#   N = 1000 
#   M <- randomLHS(N, d)
#   
#   Mprime <- randomLHS(N, d)
#   system.time({f_sens = first_order_sens(gpi, d, M, Mprime)})
#   
#   t_sens = total_sens(gpi, d, M, Mprime)
#   t_sens - f_sens


##########################################################################
# arbitrary number of variable version
# return the 1st order/total sensitivity for the group with var_indices (scalar)
##########################################################################
# first order sensitivity
first_order_sens_mult = function(gpi, d, M, Mprime, var_indices){
  N <- nrow(M)
  pM <- predict(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  S <- EE2j <- rep(NA, d)
  
  Mjprime <- Mprime
  Mjprime[,var_indices] <- M[,var_indices]
  pMprime <- predict(gpi, Mjprime, lite=TRUE, nonug=TRUE)
  EE2<- (t(pM$mean) %*% pMprime$mean)/(N - 1)
  S <- (EE2 - Ey^2)/Vary
  # since numerator and denominator are all variances in the definition and they
  # should be greater than 0.
  S[S<0] <- 0
  return(S)
}

total_sens_mult = function(gpi, d, M, Mprime, var_indices){
  N <- nrow(M)
  pM <- predict(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  
  Mj <- M
  Mj[,var_indices] <- Mprime[,var_indices]
  pMj <- predict(gpi, Mj, lite=TRUE, nonug=TRUE)
  EE2m <- (t(pM$mean) %*% pMj$mean)/(N - 1)
  T_sens <- 1 - (EE2m - Ey^2)/Vary
  # since numerator and denominator are all variances in the definition and they
  # should be greater than 0.
  T_sens[T_sens<0] <- 0
  return(T_sens)
}
