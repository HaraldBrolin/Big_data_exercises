rm(list = ls())
library(caret); library(mlbench); library(RSNNS)

load("HWD_train_data.RData")

colors <- c('white','black'); cus_col <- colorRampPalette(colors=colors)

set.seed(2017)
id_train <- sample(1:dim(data_train)[1], 0.8*dim(data_train)[1])
train_data <- data_train[id_train, ]; rownames(train_data) <- 1:dim(train_data)[1]
test_data <- data_train[-id_train, ]; rownames(test_data) <- 1:dim(test_data)[1]



grid <- expand.grid(size = c(25))
ctrl <- trainControl(method = "repeatedcv",
                     repeats = 2,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_mlp <- caret::train(y = make.names(data_train$V1),
                   x = data_train[,-1],
                   method = "mlp",
                   tuneGrid = grid,
                   trControl = ctrl)
