#-- Som perceptroner, kan testa utan training
rm(list = ls())
library(caret); library(mlbench); library(adabag); library(plyr)

load("HWD_train_data.RData")
set.seed(2017)

df_numbers <- sample(data_train)
rm(data_train)
id_train <- sample(1:dim(df_numbers)[1], 0.8*dim(df_numbers)[1])
train_data <- df_numbers[id_train, ]; rownames(train_data) <- 1:dim(train_data)[1]
test_data <- df_numbers[-id_train, ]; rownames(test_data) <- 1:dim(test_data)[1]

rm(id_train)

#------------------ Train a model

ctrl <- trainControl(method = 'repeatedcv',
                     repeats = 5,
                     number = 10,
                     classProbs = T)

M_svm <- train(make.names(V1) ~ ., train_data,
               method = 'mlp',
               metric = "Accuracy", 
               trControl = ctrl)
