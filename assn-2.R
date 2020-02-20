#################################################################
##              YSC4224: Data Science Accelerator              ##
##                Assignment 2: myroc Function                 ##
##                    by Hans Toby Limanto                     ##
#################################################################
library(pROC)

##============================================
##  1. Importing Logistic Regression Model   =
##============================================
# Importing logistic regression model, rPhi, xtest and ytest
glmresult <- read.csv("tidy_glmresult.csv")
rPhi <- as.matrix(read.csv("rPhi.csv"))[,-1]
beta.Z <- as.matrix(glmresult$estimate[-1])
beta0 <- glmresult$estimate[1] # getting the intercept
beta.X <- rPhi %*% beta.Z

# importing xtest and ytest
xtest <- as.matrix(read.csv("Xtest.csv"))[,-1]
ytest <- read.table("ytest.txt")$x

# getting the linear predictors eta and apply it to invlogit to get the probabilities
hateta <- beta0 + xtest %*% beta.X 
hatprob <- as.numeric(invlogit(hateta))

# using roc function from pROC
g <- roc(ytest ~ hatprob)
plot(g, legacy.axes = TRUE)
coords(g, "best",transpose=TRUE)

### The best threshold is 0.5829
delta <- coords(g, "best",transpose=TRUE)[1]
haty_new <- rep(0,40)
haty_new[hatprob>delta] <- 1
conf_mat <- table(ytest,haty_new)
conf_mat
auc(g)

##======================================
##  2. myroc Function Implementation   =
##======================================
# setting up the thresholds
from <- 0.4006967
to <- 0.7149774
length.out <- 80
myroc <- function(ytest, hat_prob, from, to, length.out){
  # Part 1: Threshold Selection
  thresholds <- seq(from=from, to=to, length.out=length.out)
  specificity <- vector()
  sensitivity <- vector()
  best_acc <- 0
  best_thres <- 0
  best_spec <- 0
  best_sen <- 0
  for (t in thresholds){
    # processing a single threshold and getting the specificity and sensitivity
    delta <- t
    haty_new <- rep(0, length(ytest))
    haty_new[hatprob>delta] <- 1
    conf_mat <- table(ytest,haty_new)
    specificity <- c(specificity, conf_mat[1,1]/sum(1-ytest))
    sensitivity <- c(sensitivity, conf_mat[2,2]/sum(ytest))
    
    # computing the accuracy
    accuracy <- sum(diag(conf_mat))/sum(conf_mat)
    if (accuracy > best_acc){
      best_acc <- accuracy
      best_thres <- t
      best_spec <- conf_mat[1,1]/sum(1-ytest)
      best_sen <- conf_mat[2,2]/sum(ytest)
    }
  }
  
  # Part 2: Area Under the Curve Calculation
  x <- 1 - specificity
  y <- sensitivity
  area <- 0
  for (i in 2:length(x)){
    small_area <- (x[i-1]-x[i]) * (y[i-1]+y[i]) / 2
    area <- area + small_area
  }
  return(list("best_threshold" = best_thres, 
              "specificity" = best_spec, 
              "sensitivity" = best_sen, 
              "auc" = area))
}
print(myroc(ytest, hat_prob, from, to, length.out))
# $best_threshold
# [1] 0.5836956

# $specificity
# [1] 0.9259259

# $sensitivity
# [1] 0.9230769

# $auc
# [1] 0.9358974

##================================
##  3. Comparison of Functions   =
##================================
# As we can see from the results in the previous two sections, the threshold from pROC and
# that obtained from function myroc is not very different (0.5829759 vs 0.5836956).
# As a result, the specificity and sensitivity are very similar to each other.

# However, there is a notable difference between the area under the curve of myroc and roc.
print(abs(myroc(ytest, hat_prob, from, to, length.out)$auc - g$auc))
# [1] 0.03846154

# This is because the method of computing is very primitive which is calculating the area
# of small trapezoid and summing them together instead of doing a proper integration, which
# I believe will give more accurate area under the curve.