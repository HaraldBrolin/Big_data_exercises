

# Länk:
# http://www.cs.ukzn.ac.za/~hughm/dm/content/slides07.pdf

set.seed(201606)
N <- 20
x1 <- runif(N,-1,1)
x2 <- runif(N,-1,1)
X <- cbind(x1,x2)
y <- ifelse(x2>x1,-1,1)

distance.from.plane = function(z,w,b) {
  sum(z*w) + b
  }
classify.linear = function(x,w,b) {
  distances = apply(x, 1, distance.from.plane, w, b)
  return(ifelse(distances < 0, -1, +1))
  }

euclidean.norm = function(x) {sqrt(sum(x * x))}
perceptron = function(x, y, learning.rate=1) {
  w = vector(length = ncol(x)) # initialize w
  b = 0 # Initialize b
  k = 0 # count updates
  R = max(apply(x, 1, euclidean.norm))
  made.mistake = TRUE # to enter the while loop
  while (made.mistake) {
    made.mistake=FALSE # hopefully
    yc <- classify.linear(x,w,b)
    for (i in 1:nrow(x)) {
      if (y[i] != yc[i]) {
        w <- w + learning.rate * y[i]*x[i,]
        b <- b + learning.rate * y[i]*R^2
        k <- k+1
        made.mistake=TRUE
        }
      } }
  s = euclidean.norm(w)
  return(list(w=w/s,b=b/s,updates=k))
  }

p <- perceptron(X, y)
