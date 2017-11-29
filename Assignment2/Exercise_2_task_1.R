library(readr)
library(glmnet) 

# the library for data and 'confusion matrix'
library(MASS); library(kernlab); library(caret)


df_train <- read_delim("~/Documents/BD/Assignment2/BreastCancerDataTrain.txt", " ", escape_double = FALSE, trim_ws = TRUE)
df_train_PC <- df_train
df_train_PC[, 1] <- as.numeric(df_train$Diagnosis == "M") # Malignt = 1
#df_train <- df_train[,-1]
PC <- prcomp(t(df_train_PC))
for (int in c(1,2,3,4,5,6,7,8,9,10,11)){
  varName <- paste("PC_", int, sep = "")
  assign(varName, PC$x[,int])
}

df_pc <- data.frame(PC_1, PC_2, PC_3, PC_4, PC_5)
pairs(df_pc) # Pairwise scatter plot for first 5 components

# --------------------------------------
# Task 2 # Testa Cross validations !!!!!!

df_test <- read_delim("~/Documents/BD/Assignment2/BreastCancerDataTest.txt", " ", escape_double = FALSE, trim_ws = TRUE)
lda_model <- lda(df_train$Diagnosis ~ ., df_train)
lda_pred <- predict(lda_model, df_test)

#---------------------------------------
# Task 3 (blev inget med logit-funktionen) !!! Eventuellt fel
# Köra på hela datasettet eller på training sen test
# First merge data test and train
#df_all <- rbind(df_test, df_train)
lr_model <- cv.glmnet(x = as.matrix(df_train[, -1]), 
                      y = as.matrix(df_train$Diagnosis), family = "binomial", 
                      type.measure = "class", nfolds = 10)

pre_res <- predict(lr_model, as.matrix(df_test[, -1]), s = "lambda.min", type = "class")
confusionMatrix(pre_res, df_test$Diagnosis)

#--------------------------------------


