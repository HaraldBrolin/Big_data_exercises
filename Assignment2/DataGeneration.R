set.seed(201606)
N <- 20
x1 <- runif(N,-1,1)
x2 <- runif(N,-1,1)
X <- cbind(x1,x2)
y <- ifelse(x2>x1,-1,1)