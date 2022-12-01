source("variable_importance.R")
N = 10000
ydiff <- ydiff.noRand <- matrix(NA, N, 8)
for(i in 1:8){
  y <- f_test_complete_importance(N,8,i)
  ydiff.noRand[,i] <- sum((y$Ytrue - y$Ytrue.void)^2)/var(y$Ytrue)
}
cluster <- list(1:2,3:5,6,7,8)
ydiff.cluster <- ydiff.noRand.cluster <- matrix(NA,N,5)
for(i in 1:5){
  y <- f_test_complete_importance(N,8,cluster[[i]])
  ydiff.noRand.cluster[,i] <- sum((y$Ytrue - y$Ytrue.void)^2)/var(y$Ytrue)
}

boxplot(abs(ydiff.noRand), 
        xlab = "Variables",
        ylab = "Abolute difference in the response variable")
boxplot(abs(ydiff.noRand.cluster),xaxt = "n", 
        xlab = "Clusters",
        ylab = "Abolute difference in the response variable")
axis(1, at = 1:5,labels = c("(1,2)", "(3,4,5)","6","7","8"))


# True rank:
rank_loop <- c(1,1,5,5,5,3,4,8)
rank_chain <- c(4,4,1,1,1,6,7,8)
rank_complete <- c(1,1,5,5,5,3,4,8)
rank_3way <- c(1,1,5,5,5,3,4,8)

# True selection:
true_active = c(rep(T,7),8)
# boxplot(abs(ydiff))
# boxplot(abs(ydiff.cluster))
# 
# boxplot(ydiff^2)
# boxplot(ydiff.cluster^2)
