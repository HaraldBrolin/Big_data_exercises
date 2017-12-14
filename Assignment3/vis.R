# the way you visualize the images of digit

rm(list = ls())
library(caret); library(mlbench); library(RSNNS)

load("~/Documents/BD/Big_data_exercises/Assignment3/HWD_train_data.RData")

colors <- c('white','black'); cus_col <- colorRampPalette(colors=colors)

set.seed(2017)
id_train <- sample(1:dim(data_train)[1], 0.8*dim(data_train)[1])
train_data <- data_train[id_train, ]; rownames(train_data) <- 1:dim(train_data)[1]
test_data <- data_train[-id_train, ]; rownames(test_data) <- 1:dim(test_data)[1]

z <- matrix(as.numeric(train_data[10,257:2]),16,16,byrow=T)[,16:1] # Inverts the matrix

image(t(z),col=cus_col(256))

#---------svmLinear---------

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_svm <- train(make.names(V1) ~.,
               data = train_data,
               method = "svmLinear",
               metric = "Accuracy",
               trControl = ctrl)
model_svm

#----------Gaussian Process with polynomial Kernel---------

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_gaussprPoly <- train(make.names(V1) ~.,
               data = train_data,
               method = "gaussprPoly",
               metric = "Accuracy",
               trControl = ctrl)

model_gaussprPoly

#----------Grided Gaussian Process with polynomial Kernel---------

grid <- expand.grid(degree = c(1,3), scale = c(1.542958e-03, 2.5e-04, 3.4e-05))

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_gaussprPoly <- caret::train(make.names(V1) ~.,
                           data = train_data,
                           method = "gaussprPoly",
                           tuneGrid = grid,
                           trControl = ctrl)

model_gaussprPoly

#----------Multilayer Perceptron-----

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_mlp <- train(y = make.names(train_data$V1),
                 x = train_data[,-1],
                 method = "mlp",
                 metric = "Accuracy",
                 trControl = ctrl)
model_mlp


#----------Grided Multilayer Perceptron-------

grid <- expand.grid(size = c(18,20,21))

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_mlp_grid <- train(y = make.names(train_data$V1),
                   x = train_data[,-1],
                   method = "mlp",
                   tuneGrid = grid,
                   trControl = ctrl)
model_mlp_grid

#----------Grided boot Multilayer Perceptron-------

grid <- expand.grid(size = c(18,20,21))

ctrl <- trainControl(method = "boot",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_mlp_grid <- train(y = make.names(train_data$V1),
                        x = train_data[,-1],
                        method = "mlp",
                        tuneGrid = grid,
                        trControl = ctrl)
model_mlp_grid



#----------Oblique Random Forest-----BROKE!

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_ORFlog <- train(make.names(V1) ~.,
                 data = train_data,
                 method = "ORFlog",
                 metric = "Accuracy",
                 trControl = ctrl)
model_ORFlog

#----------Multilayer Perceptron with weight decay-----

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_mlpWeightDecay <- caret::train(y = make.names(train_data$V1),
                 x = train_data[,-1],
                 method = "mlpWeightDecay",
                 metric = "Accuracy",
                 trControl = ctrl)

model_mlpWeightDecay


#----------Grided Multilayer Perceptron with weight decay-----

grid <- expand.grid(size = c(16,17,18), decay = c(0.000002321, 0.000004643))

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1, 
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_mlpWeightDecay_grid <- caret::train(y = make.names(train_data$V1),
                      x = train_data[,-1],
                      method = "mlpWeightDecay",
                      tuneGrid = grid,
                      trControl = ctrl)

model_mlpWeightDecay_grid


#----------Naive Bayes-----

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

model_naive_bayes <- caret::train(y = make.names(train_data$V1),
                                     x = train_data[,-1],
                                     method = "naive_bayes",
                                     metric = "Accuracy",
                                     trControl = ctrl)

model_naive_bayes



