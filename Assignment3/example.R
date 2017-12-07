#-------------------------------#
#--- Demo of package 'caret' ---#
#------ Xijia Liu 2017 Nov -----#
#-------------------------------#

rm(list = ls())
library(caret); library(mlbench)

data("Sonar")

head(Sonar)
pairs(Sonar[,10:15], col = Sonar$Class)

set.seed(2017)
id_train <- sample(1:dim(Sonar)[1], 0.8*dim(Sonar)[1])
train_data <- Sonar[id_train, ]; rownames(train_data) <- 1:dim(train_data)[1]
test_data <- Sonar[-id_train, ]; rownames(test_data) <- 1:dim(test_data)[1]

# Caret -> svm

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     classProbs = T)

M_svm <- train(x = train_data[,-61],
               y = train_data[,61],
               method = "svmRadial",
               metric = "Accuracy", 
               trControl = ctrl)



# two training parameter

grid <- expand.grid(sigma = c(0.01, 3), C = c(0.75, 0.9, 1, 1.1, 1.25))

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 1,
                     number = 10,
                     search = "random", 
                     classProbs = T)

M_svm <- train(x = train_data[,-61],
               y = train_data[,61],
               method = "svmRadial",
               tuneGrid = grid,
               trControl = ctrl)

M_svm

predict(M_svm, test_data)
