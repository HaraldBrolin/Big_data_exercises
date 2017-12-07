

# Länk:
# http://www.cs.ukzn.ac.za/~hughm/dm/content/slides07.pdf



set.seed(201606)
N <- 20
x1 <- runif(N,-1,1)
x2 <- runif(N,-1,1)
x <- cbind(x1,x2)
Y <- ifelse(x2>x1,-1,1)

plot(x)
points(subset(x,Y==1),col="black",pch="+",cex=2)
points(subset(x,Y==-1),col="red",pch="-",cex=2)

distance.from.plane = function(z,w,b) {
  sum(z*w) + b # If False = 0, first iteration returns b (0), andra bara -1 b är -15
  }
classify.linear = function(x,w,b) {  # x, FALSE FALSE, 0
  distances = apply(x, 1, distance.from.plane, w, b)
  return(ifelse(distances < 0, -1, +1)) #first only 1111
  }

euclidean.norm = function(x) {sqrt(sum(x * x))} # multiply row-vector, the distance of the vector
perceptron = function(x, y, learning.rate=1) {
  w = vector(length = ncol(x)) # initialize w, with length = ncol, FALSE FALSE
  b = 0 # Initialize b
  k = 0 # count updates
  R = max(apply(x, 1, euclidean.norm)) # 1 refers to the rows, finds max distance 
  made.mistake = TRUE # to enter the while loop
  while (made.mistake) {
    made.mistake=FALSE # hopefully
    yc <- classify.linear(x,w,b) # predicts, first round 111111 50%
    for (i in 1:nrow(x)) {
      if (y[i] != yc[i]) { 
        w <- w + learning.rate * y[i]*x[i,] #w ökar för varje fel prediction
        b <- b + learning.rate * y[i]*R^2 # b ökar för varje fel, minskar r 2
        k <- k+1
        made.mistake=TRUE
        }
      } }
  s = euclidean.norm(w)
  return(list(w=w/s,b=b/s,updates=k))
  }

p <- perceptron(x, Y)
intercept <- - p$b / p$w[[2]]

slope <- - p$w[[1]] /p$ w[[2]]
abline(intercept,slope,col="green")
