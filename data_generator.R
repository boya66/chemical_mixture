#### This script contains functions to simulate different types of interactions,
#### including complete graph, chain, only higher-order, loop 
####  Different scenarios are considered,
####    1. number of null variables can increase freely.
####    2. explanatory variables with/without correlations.
####    3. explanatory variables with/without main effects
####    4. change the noise level

# true function (the X8 is useless)
f_test_loop = function (n, d=8)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  return(data.frame(X, Y, Ytrue))
}

f_test_chain = function(n, d = 8){
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4])^2 + 
    10*(X[,3] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  return(data.frame(X, Y, Ytrue))
}

f_test_complete = function(n, d = 8)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 10*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 
    10 * X[,3] * abs(X[,4] - X[,5]) + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  return(data.frame(X, Y, Ytrue))
}

f_test_3way = function(n, d = 8)
{
  X <- randomLHS(n, d)
  Ytrue <- 12*sin(pi*X[,1]*X[,2]) + 
    10 *cos(pi*X[,3]) * abs(X[,4] - X[,5])+ 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 0.5)
  return(data.frame(X, Y, Ytrue))
}
