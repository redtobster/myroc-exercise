#################################################################
##              YSC4224: Data Science Accelerator              ##
##                 Assignment 2b: PC Selection                 ##
##                    by Hans Toby Limanto                     ##
#################################################################
source("assn-2.R")

##===============================
##  1. Importing Dependencies   =
##===============================
# loading the data
Xpix <- as.matrix(read.csv("dependencies/Xpix.csv"))[,-1]
print(dim(Xpix))
y <- read.table("dependencies/ytrain.txt")$x
print(length(y))

##==========================
##  1.b Helper Functions   =
##==========================

# 1) helper function to make prediction
glm_predict <- function(glm, Xpix_test, rPhi){
  beta_Z <- as.matrix(glm$coefficients[-1])
  beta0 <- glm$coefficients[1] # getting the intercept
  beta_X <- rPhi[,1:length(glm$coefficients)-1] %*% beta_Z
  
  # getting the linear predictor
  hat_eta <- beta0 + Xpix_test %*% beta_X
  return(as.numeric(invlogit(hat_eta)))
}

# 2) helper function to build the model from the Phi
build_log_reg_model <- function(Xpix_train, y_train, rPhi, num_PC){
  rPhi <- rPhi[, 1:num_PC]
  rZ <- Xpix_train %*% rPhi
  glm_result <- glm(y_train ~ rZ)
  return(glm_result)
}

# 3) helper function to calculate accuracy, specificity and sensitivity
calculate_metric <- function(ytest, yhat, thresh){
  yhat_new <- rep(0, length(yhat))
  yhat_new[yhat > thresh] <- 1
  conf_mat <- table(ytest,yhat_new)
  acc <- 1 - (sum(abs(ytest-yhat_new)) / length(ytest))
  spec <- try(conf_mat["0","0"]/sum(1-ytest))
  if (class(spec) == "try-error") {spec <- 0}
  sens <- try(conf_mat["1","1"]/sum(ytest))
  if (class(sens) == "try-error") {sens <- 0}
  return(list("accuracy" = acc, "sensitivity" = sens, "specificity" = spec))
}

# 4) helper function to preprocess and split the data
preprocess_data <- function(Xpix, y, split_ratio){
  set.seed(42)
  # splitting the data
  all_idx <- seq_len(length(y))
  idx <- sample(all_idx, length(all_idx)*split_ratio, replace = FALSE)
  Xpix_train <- Xpix[idx,]
  Xpix_test <- Xpix[-idx,]
  y_train <- y[idx]
  y_test <- y[-idx]
  
  # centering and scaling train and test data
  Xpix_train_mean <-  apply(Xpix_train, 2, mean)
  Xpix_train_sd <- apply(Xpix_train, 2, sd)
  Xpix_train <- scale(Xpix_train, center = TRUE, scale = TRUE)
  Xpix_test <- scale(Xpix_test, center = Xpix_train_mean, scale = Xpix_train_sd)
  
  return (list("Xpix_train" = Xpix_train, "Xpix_test" = Xpix_test, "y_train" = y_train, "y_test" = y_test))
}

##=====================
##  2. PC Selection   =
##=====================

get_PC_suggestions <- function(Xpix, y, split_ratio = 0.8, max_PC = 100, from = 0.4, to =0.7, length.out = 100){
  # preprocess data
  print("Preprocessing data..")
  data <- preprocess_data(Xpix, y, split_ratio)  
  
  # performing PCA
  print("Performing PCA..")
  pca <- prcomp(Xpix)
  rPhi <- pca$rotation
  
  # metrics prep
  best_acc <- list("score" = 0, "PC" = 0, "threshold" = 0)
  best_auc <- list("score" = 0, "PC" = 0, "threshold" = 0)
  best_sens <- list("score" = 0, "PC" = 0, "threshold" = 0)
  best_spec <- list("score" = 0, "PC" = 0, "threshold" = 0)
  plot <- list("acc" = vector(), "auc" = vector(), "sens" = vector(), "spec" = vector(), "y" = seq_len(max_PC))
  
  # metrics computation and PC selection
  print("Performing computations and PC selection")
  for (PC in seq_len(max_PC)){
    glm_result <- build_log_reg_model(data$Xpix_train, data$y_train, rPhi, PC)
    yhat <- glm_predict(glm_result, data$Xpix_test, rPhi)
    optimal <- myroc(data$y_test, yhat, from, to, length.out)
    thresh <- optimal$best_threshold
    
    metric <- calculate_metric(data$y_test, yhat, thresh)
    auc <- optimal$auc
    specificity <- metric$specificity
    sensitivity <- metric$sensitivity
    accuracy <- metric$accuracy
    if (accuracy > best_acc$score) {
      best_acc$score <- accuracy
      best_acc$PC <- PC
      best_acc$threshold <- thresh
    }
    if (auc > best_auc$score) {
      best_auc$score <- auc
      best_auc$PC <- PC
      best_auc$threshold <- thresh
    }
    if (sensitivity > best_sens$score) {
      best_sens$score <- sensitivity
      best_sens$PC <- PC
      best_sens$threshold <- thresh
    }
    if (specificity > best_spec$score) {
      best_spec$score <- specificity
      best_spec$PC <- PC
      best_spec$threshold <- thresh
    }
    plot$acc <- c(plot$acc, accuracy)
    plot$auc <- c(plot$auc, auc)
    plot$spec <- c(plot$spec, specificity)
    plot$sens <- c(plot$sens, sensitivity)
  }
  return(list("specificity" = best_spec, 
              "sensitivity" = best_sens, 
              "accuracy" = best_acc, 
              "auc" = best_auc,
              "plot" = plot))
}

##=============================
##  3. Plotting and Remarks   =
##=============================

PC_suggestions <- get_PC_suggestions(Xpix, y)

# 1) Accuracy metric
plot(PC_suggestions[["plot"]]$y, PC_suggestions[["plot"]]$acc,
     xlab = "Principal Components", ylab = "Score",
     main = "Accuracy Scores")
abline(v=PC_suggestions[["accuracy"]]$PC, col = "red", lwd = 3)

# Using the accuracy metric, we will need to select 58 principal components to get
# the maximum accuracy. Although the score we get is pretty high, selecting a lot of principal
# components will definitely make interpreting the results impossible. However, since the data
# that we are dealing with is pixel, we do not really need to focus much on interpreting individual
# pixels as a feature of the image.

# 2) Sensitivity metric
plot(PC_suggestions[["plot"]]$y, PC_suggestions[["plot"]]$sens,
     xlab = "Principal Components", ylab = "Score",
     main = "Sensitivity Scores")
abline(v=PC_suggestions[["sensitivity"]]$PC, col = "red", lwd = 3)

# Using the sensitivity (true positive) metric in my opinion is the best when the data is not
# balanced. However, a score of 1 in my opinion is too suspicious. There might be some overfitting
# and this causes predicting the absence of glasses to be poorly done. As a result, there
# will be too many pictures having labelled as glasses

# 3) Specificity metric
plot(PC_suggestions[["plot"]]$y, PC_suggestions[["plot"]]$spec,
     xlab = "Principal Components", ylab = "Score",
     main = "Specificity Scores")
abline(v=PC_suggestions[["specificity"]]$PC, col = "red", lwd = 3)

# Using specificity metric is also rather tricky. This is because when the principal component
# is small, the model does not have enough information to decipher whether an image contains
# glasses or not. So given the imbalanced data, the model will only predict that the image
# does not have glasses which cause specificity score to be "excellent".

# 4) AUC metric
plot(PC_suggestions[["plot"]]$y, PC_suggestions[["plot"]]$auc,
     xlab = "Principal Components", ylab = "Score",
     main = "AUC Scores")
abline(v=PC_suggestions[["auc"]]$PC, col = "red", lwd = 3)

# I am not very sure why the AUC score remains the same regardless of the number of principal
# components considered. But I am guessing that the actual differences within area under 
# the curve is very small. And using an estimate figure might cause the calculations to be all
# very similar

# All in all, I believe that from this exercise, the accuracy metric is the best metric
# to select the number of principal components.