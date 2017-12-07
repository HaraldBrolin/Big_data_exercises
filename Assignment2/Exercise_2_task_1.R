library(readr)
library(glmnet) 

# the library for data and 'confusion matrix'
library(MASS); library(kernlab); library(caret)

set.seed(2017)

BreastCancerDataTrain <- read.csv("BreastCancerDataTrain.txt", header = TRUE, sep = " ")
#df_test <- read.csv("BreastCancerDataTest.txt", header = TRUE, sep = " ")


df_test <- sample(nrow(BreastCancerDataTrain), nrow(BreastCancerDataTrain) * 0.2) # Selects 20% as test data
df_train <- BreastCancerDataTrain[-df_test,]
df_test <- BreastCancerDataTrain[df_test,]

train_resp <- as.numeric(df_train$Diagnosis == "M") # Malignt = 1
test_resp <- as.numeric(df_test$Diagnosis == "M")

# For a better PCA we should scale the data
#df_train_norm <- as.data.frame(scale(df_train[, -1]))
#df_test_norm <- as.data.frame(scale(df_test[, -1]))

#------------- Name is still norm, but nothing has been done
df_train_norm <- as.data.frame(df_train[, -1])
df_test_norm <- as.data.frame(df_test[, -1])

df_train_norm <- cbind(train_resp, df_train_norm)
names(df_train_norm)[1] <- "Diagnosis"

df_test_norm <- cbind(test_resp, df_test_norm)
names(df_test_norm)[1] <- "Diagnosis"



PC <- prcomp(t(BreastCancerDataTrain[,-1])) # Nu kör vi på hela datasettet inte bara train
for (int in c(1,2,3,4,5,6,7,8,9,10)){
  varName <- paste("PC_", int, sep = "")
  assign(varName, PC$x[,int])
}

df_pc <- data.frame(PC_1, PC_2, PC_3, PC_4, PC_5)
pairs(df_pc, cex =1.2) # Pairwise scatter plot for first 5 components
ggpairs(df_pc)

# --------------------------------------
# Task 2 # Testa Cross validations !

#df_test <- read_delim("~/Documents/BD/Assignment2/BreastCancerDataTest.txt", " ", escape_double = FALSE, trim_ws = TRUE)
lda_model <- lda(Diagnosis ~ ., df_train)
lda_pred <- predict(lda_model, df_test)
#table(lda_pred$class, df_test$Diagnosis) # See the accuracy

confusionMatrix(lda_pred$class, df_test$Diagnosis) 

#---------------------------------------
# Task 3 (blev inget med logit-funktionen) !!! Eventuellt fel
# KÃ¶ra pÃ¥ hela datasettet eller pÃ¥ training sen test
# First merge data test and train
#df_all <- rbind(df_test, df_train)
#lr_model <- cv.glmnet(x = as.matrix(df_train[, -1]), 
#                      y = as.matrix(df_train_norm$Diagnosis), family = "binomial", 
#                      type.measure = "class", nfolds = 10)

# Kör crossvalidation. bättre med bara glmnet

lr_model <- glmnet(x = as.matrix(df_train_norm[which(colnames(df_test) != "Diagnosis")]),
                    y = as.matrix(df_train_norm$Diagnosis),
                    family = "binomial") 


pre_res <- predict(lr_model, as.matrix(df_test_norm[which(colnames(df_test)!="Diagnosis")]), type = "class") # Here is where i get error!!!!!!!!!!!!!!!!!!!!!!!!!
confusionMatrix(pre_res, df_test_norm$Diagnosis)

#--------------------------------------
# Task 4 
# KKN

df_train$Diagnosis <- as.numeric(df_train$Diagnosis == "M")
df_test$Diagnosis <- as.numeric(df_test$Diagnosis == "M")
  
#knn_model <- knn3(df_train$Diagnosis ~ ., df_train, k = 5)
knn_model <- knn3Train(df_train, df_test, df_train$Diagnosis, k = 5, l = 0, prob = TRUE, use.all = TRUE)
table(knn_model[1:100], df_test$Diagnosis)
