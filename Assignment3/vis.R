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

M_svm <- train(make.names(V1) ~.,
               data = train_data,
               method = "svmLinear",
               metric = "Accuracy",
               trControl = ctrl)
M_svm

#----------Gaussian Process with polynomial Kernel--------- DONE!

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

M_svm <- train(make.names(V1) ~.,
               data = train_data,
               method = "gaussprPoly",
               metric = "Accuracy",
               trControl = ctrl)
M_svm

#----------Multilayer Perceptron-----

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     number = 10,
                     search = "random", 
                     classProbs = T)

M_svm_1 <- train(y = make.names(data_train$V1),
                 x = data_train[,-1],
                 method = "mlp",
                 metric = "Accuracy",
                 trControl = ctrl)
M_svm_1


#----------Oblique Random Forest-----BROKE!

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

oblique_RF <- train(make.names(V1) ~.,
                 data = train_data,
                 method = "ORFlog",
                 metric = "Accuracy",
                 trControl = ctrl)
M_svm_1



library(tensorflow)
install_tensorflow()

sess = tf$Session()
hello <- tf$constant('Hello, TensorFlow!')
sess$run(hello)
