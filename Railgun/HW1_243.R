#============== Homework 1 R code ==================#

### Randomized matrix multiplication ###
# 1. Implement the algorithm presented in class for randomized matrix multiplication
RandomMatrixMultiplication = function(A,B,r){
  prob = c()
  l = dim(A)[2]
  for (i in 1:l){
    prob = c(prob, norm(matrix(A[,i]),type = c("2")) * norm(matrix(B[i,]),type = c("2") ))
  }
  prob = prob / sum(prob)
  
  S = 0
  for (i in 1:r){
    selected = sample(1:l,size = 1, replace = T, prob = prob)
    S = S + outer(A[,selected], B[selected,]) / prob[selected]
  }
  return(S / r)
}

# 2. Apply the algorithm to the provided matrices
library (readr)
A = read_csv("https://raw.githubusercontent.com/zhikuanquan/UCD-STA/master/STA243/HW1/STA243_homework_1_matrix_A.csv", col_names = FALSE)
B = read_csv("https://raw.githubusercontent.com/zhikuanquan/UCD-STA/master/STA243/HW1/STA243_homework_1_matrix_B.csv", col_names = FALSE)
A = as.matrix(A)
B = as.matrix(B)
Approx20 = RandomMatrixMultiplication(A,B,20)
Approx50 = RandomMatrixMultiplication(A,B,50)
Approx100 = RandomMatrixMultiplication(A,B,100)
Approx200 = RandomMatrixMultiplication(A,B,200)

# 3. Calculate the relative approximation error
ApproximateError = function(A,B,E){
  return(norm(E - A%*%B,c("F")) / (norm(A,c("F")) * norm(B,c("F"))))
}
error20 = ApproximateError(A,B,Approx20)
error50 = ApproximateError(A,B,Approx50)
error100 = ApproximateError(A,B,Approx100)
error200 = ApproximateError(A,B,Approx200)

# 4. Visualize the estimates
par(mfrow=c(2,2))
image(z = Approx20,main = "Estimates of r = 20")
image(z = Approx50,main = "Estimates of r = 50")
image(z = Approx100,main = "Estimates of r = 100")
image(z = Approx200,main = "Estimates of r = 200")



### Power MEthod ###
power_iteration = function(A, v0, eps = 1e-6, maxiter=100) {
  # Please implement the function power_iteration that takes 
  # in the matrix X and initial vector v0 and returns the eigenvector.
  v = v0
  for (i in 1:maxiter){
    c_t = A %*% v
    c_t = c_t / (norm(as.matrix(c_t),c("2")) + eps)
    v = c_t
  }
  return(c_t)
}

set.seed(5)
E = matrix(rnorm(100), 10, 10)
v = c(1, rep(0, 9))
lams = 1:10
prods = c()
for (lambda in lams) {
  X = lambda*outer(v, v) + E
  v0 = rep(1, nrow(E))
  v0 = v0/sqrt(sum(v0^2))
  vv = power_iteration(X, v0)
  prods = c(prods, abs(v %*% vv))
}
plot(lams, prods, "b")



###  Sketched-OLS  ###
library(pracma)
# 1. Algorithms 
sketched_OLS = function(X, y, error){
  # Preparation for X and y by resampling methods
  X = as.matrix(X)
  y = as.matrix(y)
  n = dim(X)[1]
  rN = 2^floor(log2(n))
  Reduce = sample(n,size = rN,replace = T)
  resamplingX = X[Reduce,]
  resamplingY = y[Reduce,]
  d = dim(X)[2]
  
  # r
  r = round(d * log(rN) / error)
  
  print("Building rd Matrix")
  rDMatrix = zeros(r,rN)
  sampleS = sample(c(1:rN),size = r,replace = T,prob = rep(1/rN,rN))
  sampleD = sample(c(1,-1),size = rN,replace = T,prob = c(0.5,0.5))
  for (i in 1:r){
    rDMatrix[i,] = fhm(generate_canonical(sampleS[i],rN) * sqrt(rN / r)) 
  }
  for (j in 1:rN){
    if (sampleD[j] == -1){
      rDMatrix[,j] = rDMatrix[,j] * -1
    }
  }
  
  #
  diaX = rDMatrix %*% resamplingX
  diaY = rDMatrix %*% resamplingY
  beta = solve(t(diaX) %*% diaX) %*% t(diaX) %*% diaY
  
  # Get the result
  result = list()
  result$beta = beta
  result$diaX = diaX
  result$diaY = diaY
  return(result)
}

# 2. Test
DesignMatrix = rand(2048,20)
Y = rand(2048,1)
error = 0.1
testResult = sketched_OLS(DesignMatrix,Y,error)



