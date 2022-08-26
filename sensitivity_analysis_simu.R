library(lhs)
library(laGP)
library(gtools) # obtain combinations
eps = sqrt(.Machine$double.eps)

################## the original code in Bobby's book ##################
fried <- function (n, m=6)
{
  if(m < 5) stop("must have at least 5 cols")
  X <- randomLHS(n, m)
  Ytrue <- 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] - 0.5)^2 + 10*X[,4] + 5*X[,5]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

data <- fried(250)
gpi <- newGPsep(as.matrix(data[,1:6]), data$Y, d=0.1, 
                g=var(data$Y)/10, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=rep(eps, 2), 
                tmax=c(10, var(data$Y)))

N <- 1000 #change to 1e4 for descent results
G <- 30  
m <- q1 <- q2 <- matrix(NA, ncol=6, nrow=G)
grid <- seq(0, 1, length=G)
XX <- matrix(NA, ncol=6, nrow=N)

for(j in 1:6) {
  for(i in 1:G) {
    XX[,j] <- grid[i]
    XX[,-j] <- randomLHS(N, 5)
    p <- predGPsep(gpi, XX, lite=TRUE, nonug=TRUE)
    m[i,j] <- mean(p$mean)
    q1[i,j] <- mean(qnorm(0.05, p$mean, sqrt(p$s2)))
    q2[i,j] <- mean(qnorm(0.95, p$mean, sqrt(p$s2)))
  }
}

plot(0, xlab="grid", ylab="main effect", xlim=c(0,1), 
     ylim=range(c(q1,q2)), type="n")
for(j in 1:6) { 
  lines(grid, m[,j], col=j, lwd=2)
  lines(grid, q1[,j], col=j, lty=2)
  lines(grid, q2[,j], col=j, lty=2) 
}
legend("bottomright", paste0("x", 1:6), fill=1:6, horiz=TRUE, cex=0.75)


##########################################################################
# single variable version
# return the 1st order/total sensitivity for all variables respectively (a d-dim vector)
# sec 8.2.3 of https://bookdown.org/rbg/surrogates/chap8.html#chap8sens
##########################################################################

# first order sensitivity
first_order_sens = function(gpi, d, M, Mprime){
  pM <- predGPsep(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  S <- EE2j <- rep(NA, d)
  for(j in 1:d) {
    Mjprime <- Mprime
    Mjprime[,j] <- M[,j]
    pMprime <- predGPsep(gpi, Mjprime, lite=TRUE, nonug=TRUE)
    EE2j[j] <- (t(pM$mean) %*% pMprime$mean)/(N - 1)
    S[j] <- (EE2j[j] - Ey^2)/Vary
  }
  return(S)
}


total_sens = function(gpi, d, M, Mprime){
  pM <- predGPsep(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  T_sens <- EE2mj <- rep(NA, d)
  for(j in 1:d) {
    Mj <- M
    Mj[,j] <- Mprime[,j]
    pMj <- predGPsep(gpi, Mj, lite=TRUE, nonug=TRUE)
    EE2mj[j] <- (t(pM$mean) %*% pMj$mean)/(N - 1)
    T_sens[j] <- 1 - (EE2mj[j] - Ey^2)/Vary
  }
  return(T_sens)
}

# true function (the X8 is useless)
f_test = function (n, d=8)
{
  X <- randomLHS(n, d)
  Ytrue <- 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

# generate data
data <- fried(250)
d = 6

# fit a surrogate
gpi <- newGPsep(as.matrix(data[,1:d]), data$Y, d=0.1, 
                g=var(data$Y)/10, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=rep(eps, 2), 
                tmax=c(10, var(data$Y))) # use MLE inference to update gpi

### test
N = 10000 
M <- randomLHS(N, d)
Mprime <- randomLHS(N, d)
system.time({f_sens = first_order_sens(gpi, d, M, Mprime)})

t_sens = total_sens(gpi, d, M, Mprime)
t_sens - f_sens


##########################################################################
# arbitrary number of variable version
# return the 1st order/total sensitivity for the group with var_indices (scalar)
##########################################################################
# first order sensitivity
first_order_sens_mult = function(gpi, d, M, Mprime, var_indices){
  pM <- predGPsep(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2
  S <- EE2j <- rep(NA, d)
  
  Mjprime <- Mprime
  Mjprime[,var_indices] <- M[,var_indices]
  pMprime <- predGPsep(gpi, Mjprime, lite=TRUE, nonug=TRUE)
  EE2<- (t(pM$mean) %*% pMprime$mean)/(N - 1)
  S <- (EE2 - Ey^2)/Vary
  return(S)
}


total_sens_mult = function(gpi, d, M, Mprime, var_indices){
  pM <- predGPsep(gpi, M, lite=TRUE, nonug=TRUE)
  Ey <- mean(pM$mean)
  Vary <- (t(pM$mean) %*% pM$mean)/N - Ey^2

  Mj <- M
  Mj[,var_indices] <- Mprime[,var_indices]
  pMj <- predGPsep(gpi, Mj, lite=TRUE, nonug=TRUE)
  EE2m <- (t(pM$mean) %*% pMj$mean)/(N - 1)
  T_sens <- 1 - (EE2m - Ey^2)/Vary
    
  return(T_sens)
}


# true function
f_test = function (n, d=8)
{
  X <- randomLHS(n, d)
  Ytrue <- 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] + X[,4] - X[,5] - 0.5)^2 + 10*X[,6] + 5*X[,7]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

f_test_linear  = function(n, d = 5)
{
  X <- randomLHS(n, d)
  Ytrue <- 10*X[,1] - 10 * X[,2] + 8 * X[,1] * X[,2] + 
    12*X[,3] + 6*X[,3]*X[,4] - 10*X[,5]*(1 - X[,4] + X[,3] * X[,4])
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

firOrd_single <- firOrd_twoWay <- firOrd_threeWay <- 
  tot_single <- tot_twoWay <- tot_threeWay <- 
  single <- twoWay <- threeWay <-
  coef <- pvalue <- NULL

for(rep in 1:20){
  set.seed(1229 + rep)
  # generate linear data
  data <- f_test_linear(500)
  d = 5
  gpi <- newGPsep(as.matrix(data[,1:d]), data$Y, d=0.1, 
                  g=var(data$Y)/10, dK=TRUE)
  mle <- mleGPsep(gpi, param="both", tmin=rep(eps, 2), 
                  tmax=c(10, var(data$Y)))
  
  N = 10000 # change to bigger number for a descent result
  M <- randomLHS(N, d)
  Mprime <- randomLHS(N, d)
  
  # single effect
  xind <- 1 : d
  firOrd_singlei <- sapply(xind, function(x) first_order_sens_mult(gpi, d, M, Mprime, x))
  tot_singlei <- sapply(xind, function(x) total_sens_mult(gpi, d, M, Mprime, x))
  singlei <- tot_singlei - firOrd_singlei
  firOrd_single <- cbind(firOrd_single,firOrd_singlei)
  tot_single <- cbind(tot_single, tot_singlei)
  singlei<- cbind(single, singlei)
  
  IND <- sapply(xind, function(x){
    lm.fit <- lm(data$Y~data[,x])
    summary(lm.fit)$coefficients[2,c(1,4)]
  })
  
  coef <- cbind(coef,IND[1,])
  pvalue <- cbind(pvalue,IND[2,])
  
  # Two-way interaction
  xind <- combinations(d,2)
  firOrd_twoWayi <- apply(xind, 1, function(x) first_order_sens_mult(gpi, d, M, Mprime, x)) 
  tot_twoWayi <- apply(xind, 1, function(x) total_sens_mult(gpi, d, M, Mprime, x))
  twoWayi <- tot_twoWayi - firOrd_twoWayi
  firOrd_twoWay <- cbind(firOrd_twoWay, firOrd_twoWayi)
  tot_twoWay <- cbind(tot_twoWay, tot_twoWayi)
  twoWay <- cbind(twoWay, twoWayi)
  
  # Three-way interaction
  xind <- combinations(d,3)
  firOrd_threeWayi <- apply(xind, 1, function(x) first_order_sens_mult(gpi, d, M, Mprime, x))
  tot_threeWayi <- apply(xind, 1, function(x) total_sens_mult(gpi, d, M, Mprime, x))
  threeWayi <- tot_threeWayi - firOrd_threeWayi
  firOrd_threeWay <- cbind(firOrd_threeWay, firOrd_threeWayi)
  tot_threeWay <- cbind(tot_threeWay, tot_threeWayi)
  threeWay <- cbind(threeWay, threeWayi)
}

plot.data <- function(result, d = 1, D = 5){
  xind = combinations(D, d)
  value = matrix(result, ncol = 1)
  dimension = apply(xind,2, function(x) rep(x,ncol(result)))
  data  = data.frame(value, dimension)
  return(data)
}
# plot the single effects
par(mfrow = c(1,3))
plotData <- plot.data(firOrd_single)
boxplot(value~dimension, plotData, 
        ylab = "First Order Sensitivity",
        xlab = "Dimension")
plotData <- plot.data(tot_single)
boxplot(value ~ dimension, plotData,
        ylab = "Total Sensitivity",
        xlab = "Dimension")
plotData <- plot.data(single)
boxplot(value ~ dimension, plotData,
        ylab = "Difference",
        xlab = "Dimension")

par(mfrow = c(1,2))
plotData <- plot.data(pvalue)
boxplot(-log10(value) ~ dimension, plotData,
        ylab = "-log10(p-value)",
        xlab = "Dimension")
plotData <- plot.data(coef)
boxplot(value ~ dimension, plotData,
        ylab = "Coefficients",
        xlab = "Dimension")

# plot the two-way interaction
plotSingle <- plot.data(single)
plotData <- plot.data(twoWay,d = 2)
plotData <- plotData[plotData$X1 == 1,]
plotData$dimension = paste(plotData$X1, plotData$X2,sep = ":")
plotData <- rbind(plotSingle, plotData[,c("value", "dimension")])
boxplot(value ~ dimension, plotData,
        xlab = "Dimension",
        ylab = "Difference")
# plot the three-way interaction
plotTwoway <- plot.data(twoWay,2)
plotTwoway$dimension = paste(plotTwoway$X1, plotTwoway$X2, sep = ":")
plotTwoway <- plotTwoway[plotTwoway$dimension == "1:2", c("value", "dimension")]

plotData <- plot.data(threeWay,3)
plotData$dimension = paste(plotData$X1, plotData$X2, sep = ":")
plotData <- plotData[plotData$dimension == "1:2",]
plotData$dimension = paste(plotData$dimension, plotData$X3, sep = ":")
plotData <- plotData[,c("value","dimension")]
plotData <- rbind(plotTwoway, plotData)

boxplot(value ~ dimension, plotData,
        xlab = "Dimension",
        ylab = "Difference")
par(mfrow = c(1,2))
# plot the two-way interaction
plotSingle <- plot.data(single)
plotSingle <- plotSingle[!(plotSingle$dimension %in% c(1,2)),]
plotData <- plot.data(twoWay,d = 2)
plotData <- plotData[plotData$X1 == 3,]
plotData$dimension = paste(plotData$X1, plotData$X2,sep = ":")
plotData <- rbind(plotSingle, plotData[,c("value", "dimension")])
boxplot(value ~ dimension, plotData,
        xlab = "Dimension",
        ylab = "Difference")
# plot the three-way interaction
plotTwoway <- plot.data(twoWay,2)
plotTwoway$dimension = paste(plotTwoway$X1, plotTwoway$X2, sep = ":")
plotTwoway <- plotTwoway[plotTwoway$dimension == "3:5", c("value", "dimension")]

plotData <- plot.data(threeWay,3)
plotData$dimension = paste(plotData$X1, plotData$X2, sep = ":")
plotData <- plotData[plotData$dimension == "3:4",]
plotData$dimension = paste(plotData$dimension, plotData$X3, sep = ":")
plotData <- plotData[,c("value","dimension")]
plotData <- rbind(plotTwoway, plotData)

boxplot(value ~ dimension, plotData,
        xlab = "Dimension",
        ylab = "Difference")
# generate data
data <- f_test(550)
d = 8
gpi <- newGPsep(as.matrix(data[,1:d]), data$Y, d=0.1, 
                g=var(data$Y)/10, dK=TRUE)
mle <- mleGPsep(gpi, param="both", tmin=rep(eps, 2), 
                tmax=c(10, var(data$Y)))

N = 20000 # change to bigger number for a descent result
M <- randomLHS(N, d)
Mprime <- randomLHS(N, d)
# f_sens = first_order_sens_mult(gpi, d, M, Mprime, c(1,2))
# t_sens = total_sens_mult(gpi, d, M, Mprime, c(1,2))


- first_order_sens_mult(gpi, d, M, Mprime, c(1)) + total_sens_mult(gpi, d, M, Mprime, c(1))
- first_order_sens_mult(gpi, d, M, Mprime, c(2)) + total_sens_mult(gpi, d, M, Mprime, c(2))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,2)) + total_sens_mult(gpi, d, M, Mprime, c(1,2))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,3)) + total_sens_mult(gpi, d, M, Mprime, c(1,3))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,4)) + total_sens_mult(gpi, d, M, Mprime, c(1,4))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,5)) + total_sens_mult(gpi, d, M, Mprime, c(1,5))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,6)) + total_sens_mult(gpi, d, M, Mprime, c(1,6))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,7)) + total_sens_mult(gpi, d, M, Mprime, c(1,7))
- first_order_sens_mult(gpi, d, M, Mprime, c(1,8)) + total_sens_mult(gpi, d, M, Mprime, c(1,8))

# issues:
# the mult version simply changes single column swap to multiple column swap
# the approx accuracy/resolution would decrease as var_indices grows larger

#TODO: finish the serial framework to cluster vars (bad idea!)
# 1. start from the var with highest 1st order sens
# 2. iterate the rest of vars, find the one that interact with the first var mostly
# 3. iterate the rest of vars, find the one that interact with the first & 2nd vars mostly

# Note:
# essentially, the goal is to find a partition of the vars, which will minimize the sum of (T-S)_i,
# where i = 1..,#of groups