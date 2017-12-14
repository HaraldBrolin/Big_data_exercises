rm(list = ls())
library(caret); library(mlbench); library(RSNNS)

load("HWD_test_data.RData")
load("Model_mlp_25_final.RData")


train_error <- model_mlp_25
test_error <-  caret::confusionMatrix(predict(model_mlp_25, data_test), make.names(data_test$V1))