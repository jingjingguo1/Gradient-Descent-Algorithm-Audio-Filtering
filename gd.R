###Load original audio files:
library(audio)
s1 <- load.wave("mike1.wav")
s2 <- load.wave("mike2.wav")
s3 <- load.wave("mike3.wav")

x1 <- matrix(s1, nrow=1)
x2 <- matrix(s2, nrow=1)
x3 <- matrix(s3, nrow=1)

X <- rbind(x1,x2,x3)
COV <- diag(3)
for(i in 1:3){
  for(j in 1:3){
    COV[i,j] <- cov(X[i,],X[j,])
  }
}
print(COV)

#scatter plot
# d12 <- data.frame(x1=head(t(x1)), x2=head(t(x2)))
# d23 <- data.frame(x2=head(t(x2)), x3=head(t(x3)))
# d31 <- data.frame(x3=head(t(x3)), x1=head(t(x1)))

par(mfrow=c(1,3)) 
# plot(x1[1:100],x2[1:100])
# plot(x2[1:100],x3[1:100])
# plot(x3[1:100],x1[1:100])
plot(x1,x2)
plot(x2,x3)
plot(x3,x1)

### Calculate Square root of the covariance matrix, sqCov:
library(expm)
options( warn = -1 )
sqCov <- sqrtm(COV)
X_white <- solve(sqCov) %*% X
#cov of X_white
sCov <- diag(3)
for(i in 1:3){
  for(j in 1:3){
    sCov[i,j] <- cov(X_white[i,],X_white[j,])
  }
}
sCov <- round(sCov, digits = 6)
print(sCov)


# Gradient Descent algorithm
gd <- function(X) {
  T <- ncol(X)
  print(T)
  # initialization
  eta <- 1
  tol <- 1e-2 # stopping criterion
  
  n_iter <- 1 
  
  set.seed(2)
  W_new <- matrix(runif(9), nrow=3)
  # W_new[lower.tri(W_new)] = t(W_new)[lower.tri(W_new)]
  # W_new <- matrix(c(1,6,3,7,2,1,9,8,2), nrow=3)
  W_old <- matrix(0,nrow=3 , ncol=3)
  
  logP <- list()
  while( norm(W_new - W_old, type="2") > tol ){
    W_old <- W_new
    W_new <- W_old + eta *( t(solve(W_old)) + 1/T * (-tanh(W_old %*% X) ) %*% t(X) )
    #logP[n_iter] <- T * log( abs( det(W_new) ) ) + sum( log( y(W_new%*%X) ))
    logP[n_iter] <- T*log(abs(det(W_new)))+sum( log( 1/(pi*cosh( W_new%*%X )) ) )
    # logP[n_iter] <- sum( log( 1/(pi*cosh( W_new%*%X )) ) )
    n_iter <- n_iter + 1
  }
  dM <-W_new %*% solve(sqCov)
  return(list(n=n_iter, W_h =W_new, l=logP, demixM=dM))
}

outs <- gd(X_white)
print(outs[c(1,2,4)])

maxL <- length(outs$l)
plot(x=1:maxL, y=outs$l, xlab="N", ylab="log Likelihood")


### Recover Source:

yhat <- recoverX
COV <- diag(3)
for(i in 1:3){
  for(j in 1:3){
    COV[i,j] <- cov(X[i,],X[j,])
  }
}
