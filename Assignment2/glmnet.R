#--------------------------------#
#--- Demo of package 'glmnet' ---#
#------ Xijia Liu 2017 Nov ------#
#--------------------------------#

library(glmnet); 

# the library for data and 'confusion matrix'
library(MASS); library(kernlab); library(caret)


# Main function in 'glment'

## ?glmnet
## ?cv.glmnet

# 1) regression + panelty

data <- Boston
dim(data)
head(data)

## try linear regression

m <- lm(medv~., data = data)
summary(m)

## Split data

set.seed(201711)

id_train <- sample(1:dim(data)[1], dim(data)[1]*0.8)
data_train <- data[id_train, ]
data_test <- data[-id_train, ]

## fit a lasso

model_0 <- glmnet(x = as.matrix(data_train[,-14]), 
                  y = data_train$medv,
                  family = "gaussian",
                  alpha = 1,
                  lambda = 1)

coef(model_0)

## cross validation 

model_1 <- cv.glmnet(x = as.matrix(data_train[,-14]), 
                     y = data_train$medv,
                     family = "gaussian",
                     alpha = 1,
                     #lambda = shrink_par,
                     #nfolds = dim(data_train)[1]
                     nfolds = 10)

plot(model_1)

model_1$lambda
model_1$lambda.min
coef(model_1, s = "lambda.min")

# 2) Classification, logistic regression with panelty

data(spam); data <- spam
head(data)
dim(data)

id <- sample(1:dim(data)[1], 100)
pairs(data[id,1:10], col=data$type[id])

## split data 

set.seed(2017)
id_train <- sample(1:dim(data)[1], 0.7*dim(data)[1])
data_train <- data[id_train, ]
data_test <- data[-id_train, ]

model_1 <- cv.glmnet(x = as.matrix(data_train[, 1:57]),
                     y = data_train$type,
                     family = "binomial",
                     type.measure = "class",
                     nfolds = 10)

plot(model_1)
pre_res <- predict(model_1, as.matrix(data_test[,-58]), s = "lambda.min", type = "class")
confusionMatrix(pre_res, data_test$type)
