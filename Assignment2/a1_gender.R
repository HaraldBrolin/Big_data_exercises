#All functions


#function for solving t-test problem
pVtest <- function(tA, tB){
  t <- try(t.test(tA,tB,silent = TRUE))
  if (is(t, "try-error")) return(NA) else return(t$p.value)
}

#Function for calculating mean c-value
cVal <- function(x, print = TRUE) {
  amean <- x[3,1]
  bmean <- x[3,2]
  cVal <- (amean + bmean)/2
  return(cVal)
}

#Calculates the cc-value for a principal component
cValCal <- function(pc, pcStat, mdAns, c){
  aC <- 0
  bC <- 0
  if (pcStat[3,1] < pcStat[3,2]){
    a <- 0
  } else{
    a <- 1
  }
  for (i in 1:length(pc)){
    if (pc[i] < c){
      if(mdAns$gender[i] == 'male' && a == 0){
        aC <- aC + 1
      }else if (mdAns$gender[i] == 'female' && a == 1){
        bC <- bC + 1
      }else{next}
    }else{ #pc[i] >= c
      if(mdAns$gender[i] == 'female' && a == 0){
        aC <- aC + 1
      }else if (mdAns$gender[i] == 'male' && a == 1){
        bC <- bC + 1
      }else{next}
    }
  }
  cc <- (aC + bC) / length(pc)
  return (cc)
}


#Function for finding the c-value with the highest correction rate
bestC <- function(pc, pcstat, mdtrain, cvalues){
  bestC <- -Inf
  rowIndex <- 0
  for (col in seq(1,16,1)){
    newC <- cValCal(pc, pcstat, mdtrain, cvalues[,col])
    if (newC > bestC){
      bestC <- newC
      rowIndex <- col
    }
  }
  return(c(bestC, rowIndex))
}

#calculates c-val using line regression
lineCval <- function(pc1, pc2){
  l1<- summary(lm(pc1[MD_train$gender == "male"] ~ pc2[MD_train$gender == "male"]))
  l2<- summary(lm(pc1[MD_train$gender == "female"] ~ pc2[MD_train$gender == "female"]))
  inter_A <- l1$coefficients[1]
  slope_A <- l1$coefficients[2]
  inter_B <- l2$coefficients[1]
  slope_B <- l2$coefficients[2]
  x_a <- pc1[MD_train$gender == "male"]
  x_b <- pc1[MD_train$gender == "female"]
  
  if(slope_A > slope_B){
    aS <- slope_A
    bS <- slope_B
  }else{
    aS <- slope_B
    bS <- slope_A
  }
  if(inter_A > inter_B){
    aI <- inter_A
    bI <- inter_B
  }else{
    aI <- inter_B
    bI <- inter_A
  }
  slopes <- seq(bS - pi/2, aS + pi/2, by = (aS - bS)/9)
  ints <- seq(bI, aI, by = (aI - bI)/9)
  ccVec <- c(0, 0, 0)
  names(ccVec) <- c("cc-val", "intercept", "slope")
  c <- 0
  
  for (i in slopes) {
    for (j in ints) {
      y4 <- i*x_a+j
      y5 <- i*x_b+j
      resi_a <- pc2[MD_train$gender == "male"] - y4
      resi_b <- pc2[MD_train$gender == "female"] - y5
      if (sum(resi_a>0)-sum(resi_a<0)/length(x_a) > sum(resi_b>0)-sum(resi_b<0)/length(x_b)) {
        cNew <- (sum(resi_a>0)+sum(resi_b<0))/(length(x_a)+length(x_b))
      }else { #((sum(resi_A>0)-sum(resi_A<0))/length(id_typeA_train) < (sum(resi_B>0)-sum(resi_B<0))/length(id_typeB_train)) {
        cNew <- (sum(resi_a<0)+sum(resi_b>0))/(length(x_a)+length(x_b))
      }
      
      if(cNew > c){
        ccVec[1] <- cNew
        ccVec[2] <- j
        ccVec[3] <- i
        c <- cNew
      }
      #print(c(i,j))
    } 
  }
  return(ccVec)
}

lineCEval <- function(pc1, pc2, inter, slope){
  x_a <- pc1[MD_test$gender == "male"]
  x_b <- pc1[MD_test$gender == "female"]
  
  ccVec <- c(0, 0, 0)
  names(ccVec) <- c("cc-val", "intercept", "slope")
  cNew <- 0
  
  y4 <- slope*x_a+inter
  y5 <- slope*x_b+inter
  resi_a <- pc2[MD_test$gender == "male"] - y4
  print(resi_a)
  resi_b <- pc2[MD_test$gender == "female"] - y5
  print(resi_b)
  if (sum(resi_a>0)-sum(resi_a<0)/length(x_a) > sum(resi_b>0)-sum(resi_b<0)/length(x_b)) {
    cNew <- (sum(resi_a>0)+sum(resi_b<0))/(length(x_a)+length(x_b))
  }else { #((sum(resi_A>0)-sum(resi_A<0))/length(id_typeA_train) < (sum(resi_B>0)-sum(resi_B<0))/length(id_typeB_train)) {
    cNew <- (sum(resi_a<0)+sum(resi_b>0))/(length(x_a)+length(x_b))
  }
  
  ccVec[1] <- cNew
  ccVec[2] <- inter
  ccVec[3] <- slope
  return(ccVec)
}

# Import data
library(readr)
GeneExpressionData <- read_delim("C:/Users/Rasmus/Desktop/Big Data/Assignment1/GeneExpressionData.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
View(GeneExpressionData)
MetaData <- read_delim("C:/Users/Rasmus/Desktop/Big Data/Assignment1/MetaData.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
View(MetaData)

GeneExpressionData <- GeneExpressionData[,-1]

# The script below creates the training sets (GED_train, MD_train) and the test sets 	
# (GED_test, MD_test), where GED and MD stands for Gene Expression Data and 
# Meta Data respectvely.
set.seed(2)
id_tyepA <- which(MetaData$gender == "male")
id_typeA_train <- sample(id_tyepA, 86)
id_tyepB <- which(MetaData$gender == "female")
id_typeB_train <- sample(id_tyepB, 64)
id_train <- c(id_typeA_train, id_typeB_train)
GED_train <- GeneExpressionData[, id_train]
GED_test <- GeneExpressionData[,-id_train]
MD_train <- MetaData[id_train, ]
MD_test <- MetaData[-id_train, ]


#Assignment of 10 PC
PC <- prcomp(t(GED_train))
for (int in c(1,2,3,4,5,6,7,8,9,10)){
  varName <- paste("PC_", int, sep = "")
  assign(varName, PC$x[,int])
}

#Values from the boxplot
PC1stat <- boxplot(PC_1[MD_train$gender=="male"], PC_1[MD_train$gender=="female"], names = c("male", "female"))
PC2stat <- boxplot(PC_2[MD_train$gender=="male"], PC_2[MD_train$gender=="female"], names = c("male", "female"))
PC3stat <- boxplot(PC_3[MD_train$gender=="male"], PC_3[MD_train$gender=="female"], names = c("male", "female"))
PC4stat <- boxplot(PC_4[MD_train$gender=="male"], PC_4[MD_train$gender=="female"], names = c("male", "female"))
PC5stat <- boxplot(PC_5[MD_train$gender=="male"], PC_5[MD_train$gender=="female"], names = c("male", "female"))
PC6stat <- boxplot(PC_6[MD_train$gender=="male"], PC_6[MD_train$gender=="female"], names = c("male", "female"))
PC7stat <- boxplot(PC_7[MD_train$gender=="male"], PC_7[MD_train$gender=="female"], names = c("male", "female"))
PC8stat <- boxplot(PC_8[MD_train$gender=="male"], PC_8[MD_train$gender=="female"], names = c("male", "female"))
PC9stat <- boxplot(PC_9[MD_train$gender=="male"], PC_9[MD_train$gender=="female"], names = c("male", "female"))
PC10stat <- boxplot(PC_10[MD_train$gender=="male"], PC_10[MD_train$gender=="female"], names = c("male", "female"))

#Plots the boxplots for different PC:s
boxplot(PC_1[MD_train$gender=="male"], PC_1[MD_train$gender=="female"], 
        names = c("male", "female"))
title("PC_1")
boxplot(PC_2[MD_train$gender=="male"], PC_2[MD_train$gender=="female"], 
        names = c("male", "female"))
title("PC_2")
boxplot(PC_3[MD_train$gender=="male"], PC_3[MD_train$gender=="female"], 
        names = c("male", "female"))
title("PC_3")
boxplot(PC_4[MD_train$gender=="male"], PC_4[MD_train$gender=="female"], 
        names = c("male", "female"))
title("PC_4")
boxplot(PC_5[MD_train$gender=="male"], PC_5[MD_train$gender=="female"], 
        names = c("male", "female"))
title("PC_5")
boxplot(PC_6[MD_train$gender=="male"], PC_6[MD_train$gender=="female"], names = c("male", "female"))
title("PC_6")
boxplot(PC_7[MD_train$gender=="male"], PC_7[MD_train$gender=="female"], names = c("male", "female"))
title("PC_7")
boxplot(PC_8[MD_train$gender=="male"], PC_8[MD_train$gender=="female"], names = c("male", "female"))
title("PC_8")
boxplot(PC_9[MD_train$gender=="male"], PC_9[MD_train$gender=="female"], names = c("male", "female"))
title("PC_9")
boxplot(PC_10[MD_train$gender=="male"], PC_10[MD_train$gender=="female"], names = c("male", "female"))
title("PC_10")


#Extract quantile, mean values from PCA boxplots
PC1stat <- cbind(PC1stat$stats)
PC2stat <- cbind(PC2stat$stats)
PC3stat <- cbind(PC3stat$stats)
PC4stat <- cbind(PC4stat$stats)
PC5stat <- cbind(PC5stat$stats)
PC6stat <- cbind(PC6stat$stats)
PC7stat <- cbind(PC7stat$stats)
PC8stat <- cbind(PC8stat$stats)
PC9stat <- cbind(PC9stat$stats)
PC10stat <- cbind(PC10stat$stats)



#Creates a dataframe with c-values
df <- data.frame(PC1stat[3,], PC2stat[3,], PC3stat[3,], PC4stat[3,], PC5stat[3,], PC6stat[3,], PC7stat[3,], PC8stat[3,], PC9stat[3,], PC10stat[3,])
c_values_box <- data.frame()
for (num in 1:10) {
  mini <- sort(c(df[1,num], df[2,num]), decreasing = FALSE)
  bo <- abs(mini[2]-mini[1])/15
  mini2 <- seq(mini[1],mini[2], bo)
  c_values_box <- rbind(c_values_box, mini2)
}


#Finding the best c-values
PC_1_c <- bestC(PC_1, PC1stat, MD_train, c_values_box)
PC_2_c <- bestC(PC_2, PC2stat, MD_train, c_values_box)
PC_3_c <- bestC(PC_3, PC3stat, MD_train, c_values_box)
PC_4_c <- bestC(PC_4, PC4stat, MD_train, c_values_box)
PC_5_c <- bestC(PC_5, PC5stat, MD_train, c_values_box)
PC_6_c <- bestC(PC_6, PC6stat, MD_train, c_values_box)
PC_7_c <- bestC(PC_7, PC7stat, MD_train, c_values_box)
PC_8_c <- bestC(PC_8, PC8stat, MD_train, c_values_box)
PC_9_c <- bestC(PC_9, PC9stat, MD_train, c_values_box)
PC_10_c <- bestC(PC_10, PC10stat, MD_train, c_values_box)

#Store the best c-values from each PC and rank them
collectedC <- data.frame(PC_1_c, PC_2_c, PC_3_c, PC_4_c, PC_5_c, PC_6_c, PC_7_c, PC_8_c, PC_9_c, PC_10_c)
colnames(collectedC) <- c('PC_1', 'PC_2', 'PC_3', 'PC_4', 'PC_5', 'PC_6', 'PC_7', 'PC_8', 'PC_9', 'PC_10')
rankBox <- collectedC[,order(collectedC[1,])]



#Chaos begins below


#ranks the regression plots

regDf <- data.frame(lineCval(PC_1, PC_2), lineCval(PC_1, PC_3), lineCval(PC_1, PC_4), lineCval(PC_1, PC_5), 
                    lineCval(PC_2,PC_3),lineCval(PC_2,PC_4), lineCval(PC_2, PC_5), 
                    lineCval(PC_3, PC_4), lineCval(PC_3, PC_5), 
                    lineCval(PC_4, PC_5))
colnames(regDf) <- c("PC_1 vs PC_2", "PC_1 vs PC_3","PC_1 vs PC_4","PC_1 vs PC_5",
                     "PC_2 vs PC_3","PC_2 vs PC_4","PC_2 vs PC_5",
                     "PC_3 vs PC_4","PC_3 vs PC_5",
                     "PC_4 vs PC_5")
rankReg <- regDf[,order(regDf[1,])]

#prints plots
labels <- (MD_train$gender == "male")+1
par(mfrow=c(2,5))
ans <- lineCval(PC_1, PC_2)
plot(PC_1, PC_2, col=labels, xlab="PC_1",ylab="PC_2", main= "PC_62 = a + b*PC_1")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))
ans <- lineCval(PC_1, PC_3)
plot(PC_1, PC_3, col=labels, xlab="PC_1",ylab="PC_3", main= "PC_3 = a + b*PC_1")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))
ans <- lineCval(PC_1, PC_4)
plot(PC_1, PC_4, col=labels, xlab="PC_1",ylab="PC_4", main= "PC_4 = a + b*PC_1")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))
ans <- lineCval(PC_1, PC_5)
plot(PC_1, PC_5, col=labels, xlab="PC_1",ylab="PC_5", main= "PC_5 = a + b*PC_1")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))

ans <- lineCval(PC_2, PC_3)
plot(PC_2, PC_4, col=labels, xlab="PC_2",ylab="PC_3", main= "PC_3 = a + b*PC_2")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))
ans <- lineCval(PC_2, PC_4)
plot(PC_2, PC_1, col=labels, xlab="PC_2",ylab="PC_4", main= "PC_4 = a + b*PC_2")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))
ans <- lineCval(PC_2, PC_5)
plot(PC_2, PC_4, col=labels, xlab="PC_2",ylab="PC_5", main= "PC_5 = a + b*PC_2")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))

ans <- lineCval(PC_3, PC_4)
plot(PC_3, PC_4, col=labels, xlab="PC_3",ylab="PC_4", main= "PC_4 = a + b*PC_3")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))
ans <- lineCval(PC_3, PC_5)
plot(PC_3, PC_5, col=labels, xlab="PC_3",ylab="PC_5", main= "PC_5 = a + b*PC_3")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))

ans <- lineCval(PC_4, PC_5)
plot(PC_4, PC_5, col=labels, xlab="PC_4",ylab="PC_5", main= "PC_5 = a + b*PC_4")
abline(ans[2], ans[3])
mtext(sprintf("cc = %f", round(ans[1], digits = 6), side = 3, line = 0.3))



#t-test

tA <- data.matrix(GeneExpressionData[1:1000, id_typeA_train])
tB <- data.matrix(GeneExpressionData[1:1000, id_typeB_train])

t <- t.test(tA[1,],tB[1,])

#matrix for storing result. Saves index and p-value
results <- matrix(data = NA, nrow = 1000, ncol = 2)
colnames(results) <- c("Index","p-value")
for (var in 1:1000){
  p <- pVtest(tA[var,], tB[var,])
  results[var, 1] <- var
  results[var, 2] <- p
}

#orders the matrix in increasing order
results <- results[order(results[,2], decreasing = FALSE),]
results[1:5,]

#Select five lowest. Perform Boxplot for each variable

par(mfrow = (c(1,1)))

var1 <- as.numeric(unlist(GED_train[499,]))
var2 <- as.numeric(unlist(GED_train[526,]))
var3 <- as.numeric(unlist(GED_train[763,]))
var4 <- as.numeric(unlist(GED_train[973,]))
var5 <- as.numeric(unlist(GED_train[417,]))

boxplot(var1[1:86], var1[87:150], names = c("male", "female"))
boxplot(var2[1:86], var2[87:150], names = c("male", "female"))
boxplot(var3[1:86], var3[87:150], names = c("male", "female"))
boxplot(var4[1:86], var4[87:150], names = c("male", "female"))
boxplot(var5[1:86], var5[87:150], names = c("male", "female"))

#get stats from boxplot for calculating c-values

var1stat <- (boxplot(var1[1:86], var1[87:150], names = c("male", "female")))
var2stat <- (boxplot(var2[1:86], var2[87:150], names = c("male", "female")))
var3stat <- (boxplot(var3[1:86], var3[87:150], names = c("male", "female")))
var4stat <- (boxplot(var4[1:86], var4[87:150], names = c("male", "female")))
var5stat <- (boxplot(var5[1:86], var5[87:150], names = c("male", "female")))

var1stat <- cbind(var1stat$stats)
var2stat <- cbind(var2stat$stats)
var3stat <- cbind(var3stat$stats)
var4stat <- cbind(var4stat$stats)
var5stat <- cbind(var5stat$stats)

#dataframe for storing c-values 

df <- data.frame(var1stat[3,],var2stat[3,],var3stat[3,],var4stat[3,],var5stat[3,])
c_values <- data.frame()
for (num in 1:5) {
  mini <- sort(c(df[1,num], df[2,num]), decreasing = FALSE)
  bo <- abs(mini[2]-mini[1])/15
  mini2 <- seq(mini[1],mini[2], bo)
  c_values <- rbind(c_values, mini2)
}

#calculates the best c-values for the variables

var1stat_c <- cValCal(var1, var1stat, MD_train, c_values[1,1])
var2stat_c <- cValCal(var2, var2stat, MD_train, c_values[2,1])
var3stat_c <- cValCal(var3, var3stat, MD_train, c_values[3,1])
var4stat_c <- cValCal(var4, var4stat, MD_train, c_values[4,1])
var5stat_c <- cValCal(var5, var5stat, MD_train, c_values[5,1])

#ranks the variables performance on the training data

collectedCVar <- data.frame(var1stat_c, var2stat_c, var3stat_c, var4stat_c, var5stat_c)
colnames(collectedCVar) <- c('var1', 'var2', 'var3', 'var4', 'var5')
rankBoxVar <- collectedCVar[,order(collectedCVar[1,])]#

#Test on new data

#prepares the principal components for evaluation

PCtest <- prcomp(t(GED_test))
for (int in c(1,2,3,4,5,6,7,8,9,10)){
  varName <- paste("PC_test_", int, sep = "")
  assign(varName, PCtest$x[,int])
}

#stores the boxplot stats from the principal components with test data

pc1ts <- cbind(boxplot(PC_test_1[MD_test$gender=="male"], PC_test_1[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc2ts <- cbind(boxplot(PC_test_2[MD_test$gender=="male"], PC_test_2[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc3ts <- cbind(boxplot(PC_test_3[MD_test$gender=="male"], PC_test_3[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc4ts <- cbind(boxplot(PC_test_4[MD_test$gender=="male"], PC_test_4[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc5ts <- cbind(boxplot(PC_test_5[MD_test$gender=="male"], PC_test_5[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc6ts <- cbind(boxplot(PC_test_6[MD_test$gender=="male"], PC_test_6[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc7ts <- cbind(boxplot(PC_test_7[MD_test$gender=="male"], PC_test_7[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc8ts <- cbind(boxplot(PC_test_8[MD_test$gender=="male"], PC_test_8[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc9ts <- cbind(boxplot(PC_test_9[MD_test$gender=="male"], PC_test_9[MD_test$gender=="female"], names = c("male", "female"))$stats)
pc10ts <- cbind(boxplot(PC_test_10[MD_test$gender=="male"], PC_test_10[MD_test$gender=="female"], names = c("male", "female"))$stats)


#Evaluates the c-values performance on the training data
evalbox1 <- cValCal(PC_test_1, pc1ts, MD_test, c_values_box[1,2])
evalbox2 <- cValCal(PC_test_2, pc2ts, MD_test, c_values_box[2,1])
evalbox3 <- cValCal(PC_test_3, pc3ts, MD_test, c_values_box[3,5])
evalbox4 <- cValCal(PC_test_4, pc4ts, MD_test, c_values_box[4,16])
evalbox5 <- cValCal(PC_test_5, pc5ts, MD_test, c_values_box[5,1])
evalbox6 <- cValCal(PC_test_6, pc6ts, MD_test, c_values_box[6,16])
evalbox7 <- cValCal(PC_test_7, pc7ts, MD_test, c_values_box[7,1])
evalbox8 <- cValCal(PC_test_8, pc8ts, MD_test, c_values_box[8,7])
evalbox9 <- cValCal(PC_test_9, pc9ts, MD_test, c_values_box[9,14])
evalbox10 <- cValCal(PC_test_10, pc10ts, MD_test, c_values_box[10,16])

#Stores the evaluations in a new data frame
evalbox <- data.frame(evalbox1, evalbox2, evalbox3, evalbox4, evalbox5, evalbox6, evalbox7, evalbox8, evalbox9, evalbox10)


regDfEval <- data.frame(lineCEval(PC_test_1, PC_test_2, regDf[2,1], regDf[3,1]), lineCEval(PC_test_1, PC_test_3, regDf[2,2], regDf[3,2]), lineCEval(PC_test_1, PC_test_4, regDf[2,3], regDf[3,3]), lineCEval(PC_test_1, PC_test_5, regDf[2,4], regDf[3,4]), 
                        lineCEval(PC_test_2,PC_test_3, regDf[2,5], regDf[3,5]),lineCEval(PC_test_2,PC_test_4, regDf[2,6], regDf[3,6]), lineCEval(PC_test_2, PC_test_5, regDf[2,7], regDf[3,7]), 
                        lineCEval(PC_test_3, PC_test_4, regDf[2,8], regDf[3,8]), lineCEval(PC_test_3, PC_test_5, regDf[2,9], regDf[3,9]), 
                        lineCEval(PC_test_4, PC_test_5, regDf[2,10], regDf[3,10]))
colnames(regDfEval) <- c("PC_1 vs PC_2", "PC_1 vs PC_3","PC_1 vs PC_4","PC_1 vs PC_5",
                         "PC_2 vs PC_3","PC_2 vs PC_4","PC_2 vs PC_5",
                         "PC_3 vs PC_4","PC_3 vs PC_5",
                         "PC_4 vs PC_5")

labels <- (MD_test$gender == "male")+1
par(mfrow=c(2,5))
plot(PC_test_1, PC_test_2, col=labels, xlab="PC_1",ylab="PC_2", main= "PC_2 = a + b*PC_1")
abline(regDf[2,1], regDf[3,1])
mtext(sprintf("cc = %f", round(regDfEval[1,1], digits = 6), side = 3, line = 0.3))

plot(PC_test_1, PC_test_3, col=labels, xlab="PC_1",ylab="PC_3", main= "PC_3 = a + b*PC_1")
abline(regDf[2,2], regDf[3,2])
mtext(sprintf("cc = %f", round(regDfEval[1,2], digits = 6), side = 3, line = 0.3))

plot(PC_test_1, PC_test_4, col=labels, xlab="PC_1",ylab="PC_4", main= "PC_4 = a + b*PC_1")
abline(regDf[2,3], regDf[3,3])
mtext(sprintf("cc = %f", round(regDfEval[1,3], digits = 6), side = 3, line = 0.3))

plot(PC_test_1, PC_test_5, col=labels, xlab="PC_1",ylab="PC_5", main= "PC_5 = a + b*PC_1")
abline(regDf[2,4], regDf[3,4])
mtext(sprintf("cc = %f", round(regDfEval[1,4], digits = 6), side = 3, line = 0.3))


plot(PC_test_2, PC_test_3, col=labels, xlab="PC_2",ylab="PC_3", main= "PC_3 = a + b*PC_2")
abline(regDf[2,5], regDf[3,5])
mtext(sprintf("cc = %f", round(regDfEval[1,5], digits = 6), side = 3, line = 0.3))

plot(PC_test_2, PC_test_4, col=labels, xlab="PC_2",ylab="PC_4", main= "PC_4 = a + b*PC_2")
abline(regDf[2,6], regDf[3,6])
mtext(sprintf("cc = %f", round(regDfEval[1,6], digits = 6), side = 3, line = 0.3))

plot(PC_test_2, PC_test_5, col=labels, xlab="PC_2",ylab="PC_5", main= "PC_5 = a + b*PC_2")
abline(regDf[2,7], regDf[3,7])
mtext(sprintf("cc = %f", round(regDfEval[1,7], digits = 6), side = 3, line = 0.3))


plot(PC_test_3, PC_test_4, col=labels, xlab="PC_3",ylab="PC_4", main= "PC_4 = a + b*PC_3")
abline(regDf[2,8], regDf[3,8])
mtext(sprintf("cc = %f", round(regDfEval[1,8], digits = 6), side = 3, line = 0.3))

plot(PC_test_3, PC_test_5, col=labels, xlab="PC_3",ylab="PC_5", main= "PC_5 = a + b*PC_3")
abline(regDf[2,9], regDf[3,9])
mtext(sprintf("cc = %f", round(regDfEval[1,9], digits = 6), side = 3, line = 0.3))


plot(PC_test_4, PC_test_5, col=labels, xlab="PC_4",ylab="PC_5", main= "PC_5 = a + b*PC_4")
abline(regDf[2,10], regDf[3,10])
mtext(sprintf("cc = %f", round(regDfEval[1,10], digits = 6), side = 3, line = 0.3))


#Gets the same variables from the test data that was investigated during training 
Evar1 <- as.numeric(unlist(GED_test[499,]))
Evar2 <- as.numeric(unlist(GED_test[526,]))
Evar3 <- as.numeric(unlist(GED_test[763,]))
Evar4 <- as.numeric(unlist(GED_test[973,]))
Evar5 <- as.numeric(unlist(GED_test[417,]))

#Stores the boxplotstats for each variable
evar1stat <- cbind(boxplot(Evar1[MD_test$gender=="male"], Evar1[MD_test$gender=="female"], names = c("male", "female"))$stats)
evar2stat <- cbind(boxplot(Evar2[MD_test$gender=="male"], Evar2[MD_test$gender=="female"], names = c("male", "female"))$stats)
evar3stat <- cbind(boxplot(Evar3[MD_test$gender=="male"], Evar3[MD_test$gender=="female"], names = c("male", "female"))$stats)
evar4stat <- cbind(boxplot(Evar4[MD_test$gender=="male"], Evar4[MD_test$gender=="female"], names = c("male", "female"))$stats)
evar5stat <- cbind(boxplot(Evar5[MD_test$gender=="male"], Evar5[MD_test$gender=="female"], names = c("male", "female"))$stats)


#Evaluates the c-values performance on the test data
evalVar1 <- cValCal(Evar1, evar1stat, MD_test, c_values[1,1])
evalVar2 <- cValCal(Evar2, evar2stat, MD_test, c_values[2,1])
evalVar3 <- cValCal(Evar3, evar3stat, MD_test, c_values[3,1])
evalVar4 <- cValCal(Evar4, evar4stat, MD_test, c_values[4,1])
evalVar5 <- cValCal(Evar5, evar5stat, MD_test, c_values[5,1])

#stores the evaluation in a new data frame
eValVar <- data.frame(evalVar1,evalVar2,evalVar3,evalVar4,evalVar5)

