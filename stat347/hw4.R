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
library('mvtnorm')

### Problem 1
source('hw4/iris.R')
data.1 <- data.frame(data,variety)

## (i)
calc <- function(data,variety) {
    mu   <- apply(data,2,mean)
    mu.b <- apply(data,2,function(x) tapply(x,variety,mean))
    data.b <- matrix(rep(mu.b - matrix(rep(mu,each=3),ncol=4,byrow=F),each=50),ncol=4)
    S.b  <- t(data.b) %*% data.b / 2
    data.w <- data - matrix(rep(mu.b,each=50),ncol=4,byrow=F)
    S.w  <- t(data.w) %*% data.w / 147 
    return(list(mu,mu.b,S.b,S.w))
}
calc.reg <- calc(data     ,variety)
calc.log <- calc(log(data),variety)

## (ii)
eig <- eigen(solve(calc.reg[[4]]) %*% calc.reg[[3]])
H   <- eig$vectors
lam <- eig$values

## (iii)
trace <- function(X) return(sum(diag(X)))

est.theta <- function(calc.list) {
    return((trace(calc.list[[3]])/trace(calc.list[[4]])-1)*150*2/(150^2-3*50^2))
}
theta.reg <- est.theta(calc.reg)
theta.log <- est.theta(calc.log) 

## (iv)
pred.data <- rbind(c(6.3,2.9,4.9,1.7),
                   c(5.5,3.1,2.9,0.8),
                   c(5.8,3.2,3.5,1.1),
                   c(5.6,3.2,3.2,1.0)
                   )
pred <- function(new.obs,theta,lam,mu,var.means,Sigma) {
    nb  <- 50
    mub <- apply(var.means,1,function(ybar) (mu +nb*theta*ybar)/(1+nb*theta))
    sigb <- Sigma * (1 + theta/(1+nb*theta))
    ps   <- nb*dmvnorm(t(new.obs-mub),sigma=sigb)
    ps  <- c(ps,lam*dmvnorm(t(new.obs - mu), sigma=Sigma*(1+theta)))
    ps  <- ps/sum(ps)
    names(ps)[4] <- "other"
    return(ps)
}

preds <- function(new.data,theta,lam,mu,var.means,Sigma) {
    out <- t(apply(new.data,1,function(new.obs) pred(new.obs,theta,lam,mu,var.means,Sigma)))
    rownames(out) <- apply(pred.data,1,function(x) paste("(",paste(x,collapse=', '),")",sep=""))
    return(out) 
}

pred.prob.10 <- preds(new.data=pred.data,theta=10,lam=1,mu=calc.reg[[1]],var.means=calc.reg[[2]],Sigma=calc.reg[[4]])
print(xtable(pred.prob.10,digits=3),type="latex",file="hw4/A1_iv.tex")
    
## (v)
pred.prob.5 <- preds(new.data=pred.data,theta=5,lam=1,mu=calc.reg[[1]],var.means=calc.reg[[2]],Sigma=calc.reg[[4]])
print(xtable(pred.prob.5,digits=3),type="latex",file="hw4/A1_v_5.tex")
pred.prob.50 <- preds(new.data=pred.data,theta=50,lam=1,mu=calc.reg[[1]],var.means=calc.reg[[2]],Sigma=calc.reg[[4]])
print(xtable(pred.prob.50,digits=3),type="latex",file="hw4/A1_v_50.tex")

## (vii)
pdf('hw4/A1_vii.pdf')
plot(data %*% H[,1:2],xlab='First Column YH',ylab='Second Column YH',main="Projection of YH",col=rep(2:4,each=50))
points(rbind(calc.reg[[1]],calc.reg[[2]])%*%H[,1:2],col=1:4,pch=17)
legend('bottom',c('Observation','Mean','Overall','Setosa','Versicolor','Virginica'),col=c(rep('gray50',2),1:4),pch=c(1,17,rep(15,4)),ncol=3)
dev.off()

center.prob <-  pred(new.obs=calc.reg[[1]],theta=10,lam=1,mu=calc.reg[[1]],var.means=calc.reg[[2]],Sigma=calc.reg[[4]])
center.prob <- as.matrix(t(center.prob))
rownames(center.prob) <- 'Mean y'
print(xtable(center.prob,digits=3),type="latex",file="hw4/A1_vii.tex")

# (viii)
pdf('hw4/A1_viii.pdf')
plot(svd(data)$u, col=rep(2:4,each=50),main='Projection onto First Two Principal Components',ylab='Second Component',xlab='First Component')
legend('topleft',c('Setosa','Versicolor','Virginica'),col=2:4,pch=c(1,1,1))
dev.off()

# (ix)
pred.prob.log <- preds(new.data=log(pred.data),theta=theta.log,lam=1,mu=calc.log[[1]],var.means=calc.log[[2]],Sigma=calc.log[[4]])
print(xtable(pred.prob.log,digits=3),type="latex",file="hw4/A1_ix.tex")
 
### problem 5
source('hw4/birth-death.R')
devst <- function(X) {
    n <- sum(X)
    rowp <- apply(X,1,sum)/n
    colp <- apply(X,2,sum)/n
    Xp   <- outer(rowp,colp)*n
    return(2*sum(X*ifelse(X==0,0,(log(X)-log(Xp)))))
}
nomchi <- chisq.test(PFtable)$statistic
nomdev <- devst(PFtable)

# simulate tables
b <- 10000
tables <- r2dtable(b,apply(PFtable,1,sum),apply(PFtable,2,sum))
pear.boots <- unlist(mclapply(tables,function(X) chisq.test(X)$statistic,mc.cores=4))
dev.boots  <- unlist(mclapply(tables,devst,mc.cores=4))
pear.ex <- sum(pear.boots>nomchi)
dev.ex  <- sum(dev.boots>nomdev)

# plot 
mi = min(c(pear.boots,dev.boots))
ma = max(c(pear.boots,dev.boots))
xseq <- seq(from=mi,to=ma,length.out=500)
pdf('hw4/A5_chi.pdf')
hist(pear.boots,freq=F,breaks=30,main='Bootstrapped Pearson Chi-Sq vs. Expected',ylab='Density',xlab='Statistic',xlim=c(mi,ma),ylim=c(0,.03))
points(xseq,dchisq(xseq,121),type='l')
dev.off()

pdf('hw4/A5_dev.pdf')
hist(dev.boots,freq=F,breaks=30,main='Bootstrapped Deviance vs. Expected',ylab='Density',xlab='Statistic',xlim=c(mi,ma),ylim=c(0,.03))
points(xseq,dchisq(xseq,121),type='l')
dev.off()

# aggregate by month difference
douX <- cbind(PFtable,PFtable)
n <-sum(PFtable)
rowp <- apply(PFtable,1,sum)/n
colp <- apply(PFtable,2,sum)/n
Xp   <- outer(rowp,colp)*n
douXp <-cbind(Xp,Xp)

diff.tot <- rep(NULL,12)
Xp.tot <- rep(NULL,12)
for (i in 1:12) {
    diff.tot[i]<-sum(diag(douX[,i:(i+11)]))
    Xp.tot[i]<-sum(diag(douXp[,i:(i+11)]))
}
test <- diff.tot
diff.tot = c(diff.tot[8:12],diff.tot[1:7])
names(diff.tot) <- -5:6
Xp.tot = c(Xp.tot[8:12],Xp.tot[1:7])
names(Xp.tot) <- -5:6
chisq.tot <- chisq.test(diff.tot,p=Xp.tot/n)

table <- rbind(diff.tot,Xp.tot)
rownames(table) <- c('Actual','Expected')
print(xtable(table,digits=3),type="latex",file="hw4/A5_sum_diff.tex")


