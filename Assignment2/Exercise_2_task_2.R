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
PC_nam <- names(df_PC)

df_PC <- sample(df_PC)
df_PC <- sample(df_PC)
df_PC <- sample(df_PC)

df_test <- sample(nrow(df_PC), nrow(df_PC) * 0.2) # Selects 20% as test data
df_train <- df_PC[-df_test,]
df_test <- df_PC[df_test,]


#--------------- Create the logistic regression
pred_tot <- data.frame()
for (i in 1:20){
  comps <- paste(PC_nam[1:i], sep="", collapse = "+")
  pred_df <- data.frame()
  for (k in 0:17){
    temp_df_test <- df_train[(k*10)+(1:10),]
    temp_df_train <- df_train[-((k*10)+(1:10)),]
    log_model <-  glm(as.formula(paste(c("Subtype ~ ", comps ), sep = "", collapse = "")), family = binomial("logit"), temp_df_train)
    pre_log <- predict(log_model, temp_df_test, type = "response")
    conf_df <- t(as.data.frame(confusionMatrix(round(pre_log, 0), temp_df_test$Subtype)[3]))
    pred_df <- rbind(pred_df, conf_df)
  }
  pred_tot <- rbind(pred_tot, t(as.data.frame(colMeans(pred_df, na.rm = FALSE))))
  #print(colMeans(pred_df, na.rm = TRUE))
}


log_model <-  glm(as.formula("Subtype ~ PC1 + PC2"), family = binomial("logit"), df_train) # Chose using two PC components
pre_log <- predict(log_model, df_test, type = "response")
confusionMatrix(round(pre_log, 0), df_test$Subtype) # The result

#------------------------ penalized logistic regression 
# alpha default is LASSO
lgp_model <- glmnet(x = as.matrix(df_train[,-(which(colnames(df_PC)=="Subtype"))]),
                      y = as.matrix(df_train$Subtype),
                      family = "binomial")                                      

pre_lgp <- predict(lg_model, as.matrix(df_test[,-(which(colnames(df_PC)=="Subtype"))]), s = 0.1, type = "class")
confusionMatrix(pre_lgp, df_test$Subtype)

#------------------------ Crossvalidation to fond best lambda
opt_lamb <- cv.glmnet(x = as.matrix(df_train[,-(which(colnames(df_PC)=="Subtype"))]),
                      y = as.matrix(df_train$Subtype),
                      family = "binomial",
                      type.measure = "class",
                      nfolds = 10)

pre_lg_opt <- predict(opt_lamb, as.matrix(df_test[,-(which(colnames(df_PC)=="Subtype"))]), s = "lambda.min", type = "class")
confusionMatrix(pre_lg_opt, df_test$Subtype)

#----------------------- Use best lambda to create optimal model
lgp_model_optlambda <- glmnet(x = as.matrix(df_train[,-(which(colnames(df_PC)=="Subtype"))]),
                              y = as.matrix(df_train$Subtype),
                              family = "binomial") 

pre_lgp_model_optlambda <- predict(lgp_model_optlambda, as.matrix(df_test[,-(which(colnames(df_PC)=="Subtype"))]), s = opt_lamb$lambda.min, type = "class")
confusionMatrix(pre_lgp, df_test$Subtype)


#------------------------------------------------------------
#---------------------Task 2.3
#-----------------------Create a training and a test-set
df_tot <- as.data.frame(t(GeneExpressionData))
df_tot$Subtype <- MetaData$SubType

sample(df_tot)
sample(df_tot)
sample(df_tot)
df_test <- sample(nrow(df_tot), nrow(df_tot) * 0.2)
df_train <- df_tot[-df_test,]
df_test <- df_tot[df_test,]


#----------------------- Use best lambda to create optimal model
opt_lamb <- cv.glmnet(x = as.matrix(df_train[,-(which(colnames(df_tot)=="Subtype"))]),
                      y = as.matrix(df_train$Subtype),
                      family = "binomial",
                      type.measure = "class",
                      nfolds = 10)

pre_lg_opt <- predict(opt_lamb, as.matrix(df_test[,-(which(colnames(df_tot)=="Subtype"))]), s = "lambda.min", type = "class")
confusionMatrix(pre_lg_opt, df_test$Subtype)
which (coef(opt_lamb, s="lambda.min") != 0 ) # fetch indecies of the non-zero genes
coef(opt_lamb, s="lambda.min")[which (coef(opt_lamb, s="lambda.min") != 0 )] # gets the values

#------------------ Måste få ut rätt gener

