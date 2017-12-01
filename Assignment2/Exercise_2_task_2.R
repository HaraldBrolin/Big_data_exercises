library(readr)
GeneExpressionData <- read_delim("~/Desktop/Big_data/Big_data_exercises/Assignment2/GeneExpressionData.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
MetaData <- read_delim("~/Desktop/Big_data/Big_data_exercises/Assignment2/MetaData.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
set.seed(2017)


df_test <- sample(nrow(GeneExpressionData), nrow(GeneExpressionData) * 0.2)
df_train <- GeneExpressionData[-df_test,]
df_test <- GeneExpressionData[df_test,]
