### Libraries
library('xtable')
library('stats4')
library('locfit')
library('SDMTools')
library('plotrix')
library('reshape2')
library('nlme')
library('parallel')
library('faraway')

### Problem 1
source('hw4/iris.R')
data.1 <- data.frame(data,variety)

## (i)
mu   <- apply(data,2,mean)
mu.b <- apply(data,2,function(x) tapply(x,variety,mean))
data.b <- matrix(rep(mu.b - matrix(rep(mu,each=3),ncol=4,byrow=F),each=50),ncol=4)
S.b  <- t(data.b) %*% data.b / 2
data.w <- data - matrix(rep(mu.b,each=50),ncol=4,byrow=F)
S.w  <- t(data.w) %*% data.w / 147 

## (ii)
eig <- eigen(solve(S.w) %*% S.b)

H <- eig$vectors
