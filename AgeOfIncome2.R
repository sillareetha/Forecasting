library(SemiPar)
library(MASS)
par(mar=c(3.1,3.1,1.1,1.1)) 

data(age.income)
input <- age.income$age
output <- age.income$log.income

ind <- sample(2,nrow(age.income), replace = T, prob = c(0.7,0.3))
ageIncome.train <- age.income[ind==1, ]
ageIncome.test <- age.income[ind==2, ]
Sum<-0

N <- 205
n <- 148# length(ageIncome.train) # nb of observations input( X <- design point)
nbsim <- 100 #nb of simulations
sigN <- seq(0.1,5,by=0.1) #sd of noise
sd<- diag(sigN)
theta <- 1:50
t<- diag(theta)
u <- seq(20,70, length = N) #discretized input set
result<- c()
##Matern 5/2 covariance function
for (p in 1:50) {
  for (q in 1:50) {
    
    k <- function(x, xp){
      (1+sqrt(5)*(abs(x-xp))/theta[p]+
         (5*(x-xp)^2)/(3*theta[p]^2))*exp(-sqrt(5)*(abs(x-xp))/theta[p])
    }
    kx <- matrix(NA,nrow=1,ncol=n)
    for(j in 1:n){
      for(i in 1 :1){
        kx[i,j]=k(u[i],ageIncome.train$age[j])
      }
    }
    Kxx <- matrix(NA,n,n)
    for(i in 1:n){
      for(j in 1 :n){
        Kxx[i,j]=k(ageIncome.train$age[i],ageIncome.train$age[j])
      }
    }
    kxx <- matrix(NA,1,1)
    for(j in 1:1){
      for(i in 1 : 1){
        kxx[i,j]=k(u[i],u[j])
      }
    }
    C <- kxx-kx%*%chol2inv(chol(Kxx+sigN[q]^2*diag(n)))%*%t(kx)
    zeta <- kx%*%solve(Kxx+sigN[q]^2*diag(n))%*%ageIncome.train$log.income ##zeta vector of dimantion
    result<-append(result,zeta)
    print(zeta)
    M <- (ageIncome.train$log.income - zeta) ^ 2
    Sum <- Sum + M
    Mean1 <- Sum / n
       }
}
r<-matrix(result,nrow = 50)
r
min(r)

id <- which(r == min(r), arr.ind = TRUE)
id
sd[1,1]

##Matern 5/2 covariance function
sigN <- 0.1 
theta <- 50
n<- 57
u <- seq(20,70, length = N) #discretized input set

k <- function(x, xp){
  (1+sqrt(5)*(abs(x-xp))/theta+
     (5*(x-xp)^2)/(3*theta^2))*exp(-sqrt(5)*(abs(x-xp))/theta)
}
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.test$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.test$age[i],ageIncome.test$age[j])
  }
}
##zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*dig))
kxx <- matrix(NA,N,N)
for(j in 1:N){
  for(i in 1 : N){
    kxx[i,j]=k(u[i],u[j])
  }
}
C <- kxx-kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%t(kx)
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.test$log.income ##zeta vector of dimantion
Sum1<-0
for (l in 1:57) {
  M1 <- (ageIncome.test$log.income[l] - zeta[l]) ^ 2
 print(M1) 
 Sum1 <- Sum1 + M1 
}
res<-Sum1
Mean2 <- res / n

Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.test$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
#lines(u,col='red',lwd=2,lty=2)
points(ageIncome.test$age,ageIncome.test$log.income,pch=19)

#3/2 mattern
k <- function(x, xp){
  (1+(sqrt(3)*abs(x-xp)/theta))
  exp(-sqrt(3)*abs(x-xp)/theta)
}
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.test$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.test$age[i],ageIncome.test$age[j])
  }
}
kxx <- matrix(NA,N,N)

for(j in 1:N){
  for(i in 1 : N){
    kxx[i,j]=k(u[i],u[j])
  }
}
C <- kxx-kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%t(kx)
#C <- kxx-kx%*%solve(Kxx)%*%t(kx) ##conditional covariance matrix
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.test$log.income ##zeta vector of dimantion
Sum2<-0
for (l in 1:n) {
  M1 <- (ageIncome.test$log.income[l] - zeta[l]) ^ 2
  print(M1) 
  Sum2 <- Sum2 + M1 
}
res1<-Sum2
Mean3 <- res1 / n

Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.test$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
points(ageIncome.test$age,ageIncome.test$log.income,pch=19)

#expo function
k <- function(x, xp){
  exp(-(abs(x-xp))/theta)
}
  kx <- matrix(NA,nrow=N,ncol=n)
  for(j in 1:n){
    for(i in 1 :N){
      kx[i,j]=k(u[i],ageIncome.test$age[j])
    }
  }
  Kxx <- matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1 :n){
      Kxx[i,j]=k(ageIncome.test$age[i],ageIncome.test$age[j])
    }
  }
  kxx <- matrix(NA,N,N)
  for(j in 1:N){
    for(i in 1 : N){
      kxx[i,j]=k(u[i],u[j])
    }
  }
  C <- kxx-kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%t(kx)
  #C <- kxx-kx%*%solve(Kxx)%*%t(kx) ##conditional covariance matrix
  zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.test$log.income ##zeta vector of dimantion
  Sum3<-0
  for (l in 1:57) {
    M1 <- (ageIncome.test$log.income[l] - zeta[l]) ^ 2
    print(M1) 
    Sum3 <- Sum3 + M1 
  }
  res2<-Sum3
  Mean4 <- res2 / n
  
  Y <- t(mvrnorm(nbsim,zeta,C))
  matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.test$log.income))
  title(xlab = 'x',ylab = 'Y(x)',line = 2)
  lines(u,zeta,col='black',lwd=2)
  #lines(u,col='red',lwd=2,lty=2)
  points(ageIncome.test$age,ageIncome.test$log.income,pch=19)

#squared expo function
k <- function(x, xp){
  exp(-(x-xp)^2/(2*(theta^2)))
}
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.test$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.test$age[i],ageIncome.test$age[j])
  }
}
##zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*dig))
kxx <- matrix(NA,N,N)
for(j in 1:N){
  for(i in 1 : N){
    kxx[i,j]=k(u[i],u[j])
  }
}
C <- kxx-kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%t(kx)
#C <- kxx-kx%*%solve(Kxx)%*%t(kx) ##conditional covariance matrix
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.test$log.income ##zeta vector of dimantion
Sum4<-0
for (l in 1:57) {
  M1 <- (ageIncome.test$log.income[l] - zeta[l]) ^ 2
  print(M1) 
  Sum4 <- Sum4 + M1 
}
res4<-Sum4
Mean5 <- res4 / n

Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.test$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
#lines(u,col='red',lwd=2,lty=2)
points(ageIncome.test$age,ageIncome.test$log.income,pch=19)

