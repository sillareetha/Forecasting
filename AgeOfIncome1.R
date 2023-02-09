library(SemiPar)
par(mar=c(3.1,3.1,1.1,1.1)) # adapt margins
# data <- read.table('ageIncome.txt',header=T)
data(age.income)
input <- age.income$age
output <- age.income$log.income
 plot(input,output,pch='*',ylab='',xlab='')
 title(xlab='age',ylab='log income',line=2)

ind <- sample(2,nrow(age.income), replace = T, prob = c(0.8,0.2))
ageIncome.train <- age.income[ind==1, ]
ageIncome.test <- age.income[ind==2, ]

library(MASS)
N <- 500
n <- 168 # nb of observations input( X <- design point)
nbsim <- 100 #nb of simulations
sigN <- 0.5 

theta <- 20
u <- seq(20,70, length = N) #discretized input set
##Matern 5/2 covariance function
k <- function(x, xp){
  (1+sqrt(5)*(abs(x-xp))/theta+
     (5*(x-xp)^2)/(3*theta^2))*exp(-sqrt(5)*(abs(x-xp))/theta)
}
##inbult function for matern
#k <- function(h,theta){
# matern.covariance()
#}



kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.train$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.train$age[i],ageIncome.train$age[j])
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
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.train$log.income ##zeta vector of dimantion
Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.train$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
#lines(u,col='red',lwd=2,lty=2)
points(ageIncome.train$age,ageIncome.train$log.income,pch=19)

#3/2 mattern
k <- function(x, xp){
  (1+(sqrt(3)*abs(x-xp)/theta))
  exp(-sqrt(3)*abs(x-xp)/theta)
}
##inbult function for matern
#k <- function(h,theta){
# matern.covariance()
#}



kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.train$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.train$age[i],ageIncome.train$age[j])
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
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.train$log.income ##zeta vector of dimantion
Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.train$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
#lines(u,col='red',lwd=2,lty=2)
points(ageIncome.train$age,ageIncome.train$log.income,pch=19)

#expo function
k <- function(x, xp){
  exp(-(abs(x-xp))/theta)
}
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.train$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.train$age[i],ageIncome.train$age[j])
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
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.train$log.income ##zeta vector of dimantion
Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.train$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
#lines(u,col='red',lwd=2,lty=2)
points(ageIncome.train$age,ageIncome.train$log.income,pch=19)

#squared expo function
k <- function(x, xp){
  exp(-(x-xp)^2/(2*(theta^2)))
}
kx <- matrix(NA,nrow=N,ncol=n)
for(j in 1:n){
  for(i in 1 :N){
    kx[i,j]=k(u[i],ageIncome.train$age[j])
  }
}
Kxx <- matrix(NA,n,n)
for(i in 1:n){
  for(j in 1 :n){
    Kxx[i,j]=k(ageIncome.train$age[i],ageIncome.train$age[j])
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
zeta <- kx%*%chol2inv(chol(Kxx+sigN^2*diag(n)))%*%ageIncome.train$log.income ##zeta vector of dimantion
Y <- t(mvrnorm(nbsim,zeta,C))
matplot(u,Y,type='l',lwd = 1,col ='gray',lty = 1,ylab = '',xlab = '',ylim = range(ageIncome.train$log.income))
title(xlab = 'x',ylab = 'Y(x)',line = 2)
lines(u,zeta,col='black',lwd=2)
#lines(u,col='red',lwd=2,lty=2)
points(ageIncome.train$age,ageIncome.train$log.income,pch=19)
