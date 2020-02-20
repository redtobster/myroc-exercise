library(arm)
library(ggplot2)
# load testing data
test_seeds <- read.table("seed_test.txt", header = TRUE)
seeds <- read.table("seeds.txt", header=TRUE)
seeds.glm3 <- glm(SURV ~  SPECIES + CAGE + log(LIGHT) + factor(LITTER),  data=seeds, family=binomial)

# make prediction
# run a fake linear regression where the y is the test and the x is the testing data.
# in the glm func, add the option X=TRUE
test_seeds <- read.table("seeds_test.txt", header = TRUE)
ytest <- test_seeds$SURV
seeds.glm3test <- glm(SURV ~ SPECIES + CAGE + log(LIGHT) + factor(LITTER), data = test_seeds, family = binomial, x=TRUE)
is(seeds.glm3test) # covariate matrix that is used for the linear regression.

# Why to call a fake glm? so that the LITTER will have a dummy representation of the variable

hatbeta <- coefficients(seeds.glm3) # from original data
ntest <- dim(test_seeds)[1]
Xtest <- seeds.glm3test$x
hatetatest <- Xtest%*%hatbeta
myhatprob <- invlogit(hatetatest)
plot(hatporb, myhatprob);abline(0, 1, col = "red",lwd = 2)

# create a function myroc
myroc <- function(delta, hatprob, ytest, add_point=TRUE){
  # start with all prediction to be 0
  pred <- rep(0, length(ytest))
  
  # add a logic to label with 1 when the probability exceeds threshold delta
  for (i in 1:length(hatprob)){
    if (hatprob[i] > delta){ pred[i] <- 1 }
  }
  print(pred)

  # add a logic to make sure that the threshold is not too small or not to large
  if (min(hatprob) > delta || max(hatprob) < delta) {stop("change threshold! Its too small (or large)")}
  
  # compute the confusion matrix using the table function
  conf <- table(ytest, pred)
  print(conf)
  # accuracy
  accuracy <- sum(diag(conf))/sum(conf)
  cat("accuracy is:", accuracy * 100, "%")
  
  # TPR (sensitivity) true true / all true
  tpr <- conf[2,2]/sum(conf[,2])
  cat("sensitivity is:", tpr * 100, "%")
  
  # FPR (specificity) true false / all false
  fpr <- conf[1,1]/sum(conf[,1])
  cat("specificty is:", fpr * 100, "%")
  
  # add the add points logic
  return(c(tpr, 1-fpr)) 
}

tpr <- NULL
fpr <- NULL
for (delta in seq(0, 1, by = 0.02)){
  xy <- myroc(delta, myhatprob, ytest)
  tpr <- c(tpr, xy[0])
  fpr <- c(fpr, xy[1])
}
