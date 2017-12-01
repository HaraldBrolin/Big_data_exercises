library(readr); library(glmnet); library(MASS); library(kernlab); library(caret)

set.seed(42)
df <- read.csv("BreastCancerDataTrain.txt", header = TRUE, sep = " ")
df$Diagnosis <- as.numeric(df$Diagnosis == "M")
df_test <- sample(nrow(df), nrow(df) * 0.25)
df_train <- df[-df_test,]
df_test <- df[df_test,]

# ------------------------ Task 1

PC <- prcomp(t(df_train[, -1]))
for (int in c(1,2,3,4,5,6,7,8,9,10)){
  varName <- paste("PC_", int, sep = "")
  assign(varName, PC$x[,int])
}
df_pc <- data.frame(PC_1, PC_2, PC_3, PC_4, PC_5)
pairs(df_pc)

# ---- ------------------ Task 2 # Testa Cross validations !!!!!!

lda_model <- lda(df_train$Diagnosis ~ ., df_train)
lda_pred <- predict(lda_model, df_test)
#table(lda_pred$class, df_test$Diagnosis) # See the accuracy
confusionMatrix(lda_pred$class, df_test$Diagnosis)

#----------------------- Task 3

lr_model <- cv.glmnet(x = as.matrix(df_train[, -1]), 
                      y = as.matrix(df_train$Diagnosis), family = "binomial", 
                      type.measure = "class", nfolds = 10) # Kan ändra nfolds

pre_res <- predict(lr_model, as.matrix(df_test[, -1]), s = "lambda.min", type = "class")
confusionMatrix(pre_res, df_test$Diagnosis)

#----------------------- Task 4

knn_model <- knn3Train(df_train, df_test, df_train$Diagnosis, k = 5, l = 0, prob = TRUE, use.all = TRUE)
#table(knn_model[1:100], df_test$Diagnosis)
confusionMatrix(knn_model, df_test$Diagnosis)

#----------------------- Task 5 Validate and predict training set

df_test_test <- read.csv("BreastCancerDataTest.txt", header = TRUE, sep = " ")
