library('alr3')
library('faraway')
library('car')
library('MASS')

#############################
# Problem 1
#############################

### Part d)
# Load the data and create the matrices 
data('savings')
attach(savings)
n<-nrow(savings)
x<-as.matrix(data.frame(rep(1,n),savings[c(2,3,4,3,5)]))
y<-as.matrix(savings[1])
qr<-qr(x)

R<-qr.R(qr)[1:qr$rank,]
Q<-qr.Q(qr)[,1:qr$rank]
W<-R[,qr$pivot]
beta<-t(W) %*% solve(W%*%t(W)) %*% t(Q) %*% y 
beta 

### part e)

model<-lm(y~x)
summary(model)

### part f)
qr.R(qr)

#############################
# Problem 3 
#############################

### part b
# Linear model and output
model<-lm(sr ~ pop15 + pop75 + dpi + ddpi)
X <- model.matrix(model)
beta<-model$coef
sigma2 <- summary(model)$sigma^2
df <- summary(model)$df[2]

# Set up for ellipse
R<-matrix(c(0,0,1,-1,-1,0,0,0,0,1),nrow=2,ncol=5)
variance.pred <- R %*% solve(t(X) %*% X) %*% t(R) + diag(2) 
variance.mean <- R %*% solve(t(X) %*% X) %*% t(R)

mycenter <- R%*%beta
radius9<- sqrt(2* sigma2 * qf(.1,2,df, lower.tail=FALSE))
radius95<- sqrt(2* sigma2 * qf(.05,2,df, lower.tail=FALSE))  

# Draw the ellipses
# Ellipse for 90% CI
pdf("hw4_3_b_ci.pdf")
plot.new()
plot.window(xlim=c(-10,15),ylim=c(-10,15))
title(main='Prediction Interval of z', xlab='z1', ylab='z2')
car::ellipse(center=c(mycenter), shape=variance.pred, radius=radius9, lty=2, main='Confidence Interval for Prediction', col=palette()[1])
car::ellipse(center=c(mycenter), shape=variance.pred, radius=radius95,lty=3, add=TRUE,col=palette()[1])
legend('topright', legend=c('90% CI','95% CI'), lty=c(2,3))
axis(1)
axis(2)
dev.off()

### part c
n<-1000
predictions <- X %*% beta 
inv.variance.pred <- solve(variance.pred)
in.range <- function(u) {
    m<-lm(predictions+rnorm(50,0,sqrt(sigma2)) ~ 0 + X)
    p<-mvrnorm(1,R%*%m$coef,summary(m)$sigma^2 * diag(2))
    e<-mycenter - p
    return(t(e) %*% inv.variance.pred %*% e <= radius9^2)
}
boots<-sapply(1:n,in.range)
boot.prob<-sum(boots)/n
boot.prob

###########################
# Problem 4
###########################
data(wm1)
attach(wm1)

###2.13
#1
pdf("hw4_4_2_13_1_plot.pdf")
plot(RSpd,CSpd)
dev.off()
# 2
reg<-lm(CSpd~RSpd)
summary(reg)
# 3
pred<-predict(reg,data.frame(CSpd=NA,RSpd=7.4285), interval="prediction")
# 5 
m<-62039
sig2<-summary(reg)$sigma^2
bxa<-7.4285
bx<-mean(RSpd)
sxx<-sum((RSpd-mean(RSpd))^2)
n<-1116
se<-sqrt(sig2/m + sig2/n + (bxa-bx)^2/sxx)

ci = pred[1]+qt(.025,1116)*c(-se,se)

### 4.6
sample.mean<-function(n) mean(sample(CSpd,n,replace=TRUE))
means<-sapply(rep(1116,999),sample.mean)
quantile(means,c(.025,.975))
pdf('hw4_4_4pt6_qq.pdf')
qqnorm(reg$residuals,ylab="Residuals",xlab="Normal Score", main="Wind Speed Residuals")
qqline(reg$residuals)
dev.off()
