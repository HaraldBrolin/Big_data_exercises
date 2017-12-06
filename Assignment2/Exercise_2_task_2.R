library(readr); library(caret);library(glmnet); 
#GeneExpressionData <- read_delim("~/Desktop/Big_data/Big_data_exercises/Assignment2/GeneExpressionData.txt", 
#                                 "\t", escape_double = FALSE, trim_ws = TRUE)
GeneExpressionData <- read_delim("GeneExpressionData.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
#MetaData <- read_delim("~/Desktop/Big_data/Big_data_exercises/Assignment2/MetaData.txt", 
#                                 "\t", escape_double = FALSE, trim_ws = TRUE)
MetaData <- read_delim("MetaData.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
set.seed(2017)
GeneExpressionData <- GeneExpressionData[,-1]
#t_gene <- as.data.frame(t(GeneExpressionData))
#names(t_gene) <- rownames(GeneExpressionData)
#rownames(t_gene) <- names(GeneExpressionData)



PC <- prcomp(t(GeneExpressionData))
df_PC <-  as.data.frame(PC$x[,1:20])

df_PC <- cbind(df_PC,MetaData$SubType)
names(df_PC)[ncol(df_PC)] <- "Subtype"
df_PC$Subtype <- as.numeric(df_PC$Subtype == "IDHmut-codel") #Codel = 1

df_test <- sample(nrow(df_PC), nrow(df_PC) * 0.2) # Selscts 20% as test data
df_train <- df_PC[-df_test,]
df_test <- df_PC[df_test,]


#--------------- Create the logistic regression
for (k in 0:17){
  temp_df_test <- df_train[(k*10)+(1:10),]
  temp_df_train <- df_train[-((k*10)+(1:10)),]
  log_model <-  glm(Subtype ~ temp_df_train$PC1, family = binomial("logit"), temp_df_train)
  pre_log <- predict(log_model, as.data.frame(temp_df_test[,-21]), type = "response")
  confusionMatrix(pre_res, data_test$type)
  
  
}
log_model <-  glm(Subtype ~ df_train$PC1, family = binomial("logit"), df_train)
print(log_model)
for (k in 3:21){
  log_model <-  glm(Subtype ~ df_train[,2:k], family = binomial("logit"), df_train)
  print(log_model)


}




#--------------- Create the logistic regression

lg_model <- cv.glmnet(x = as.matrix(df_train[, -1]),
                     y = df_train$Subtype,
                     family = "binomial",
                     type.measure = "class",
                     nfolds = 10)
plot(lg_model)
pre_res <- predict(lg_model, as.matrix(df_test[,-1]), s = "lambda.min", type = "class")
ans <- confusionMatrix(pre_res, df_test$Subtype)


#----------------- Nu anv?nder jag minst two komponenter
for (k in 2:20){
  lg_model <- cv.glmnet(x = as.matrix((df_train[-1])[1:k]),
                        y = df_train$Subtype,
                        family = "binomial",
                        type.measure = "class",
                        nfolds = 10)
  pre_res <- predict(lg_model, as.matrix((df_test[-1])[1:k]), s = "lambda.min", type = "class")
  ans <- confusionMatrix(pre_res, df_test$Subtype)
  cat("Number of PC:", k, "\n")
  print(ans$overall)
}

#------------------ Apply to test set (best 2?)
lg_model <- cv.glmnet(x = as.vector((df_train[-1])[1:2]),
                      y = df_train$Subtype,
                      family = "binomial",
                      type.measure = "class",
                      nfolds = 10)
pre_res <- predict(lg_model, as.vector((df_test[-1])[1:2]), s = "lambda.min", type = "class")



