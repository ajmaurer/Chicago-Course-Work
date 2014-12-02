library('parallel')
library('locfit')
library('gam')
library('MASS')

############################################
# Problem 3 
############################################

fiji<-read.table('fijiquakes.dat',header=TRUE)
n<-nrow(fiji)

### a)
mag.den<-density(fiji$mag,bw="nrd",kernel="gaussian")
h<-mag.den$bw

pdf("hw3_3_a.pdf")
plot(mag.den,main=paste("Fiji Earthquake Density Estimates, Bandwidth of",format(h,digits=3)),xlab="Magnitude", ylab="Density",lty=1)
dev.off()

### b)
mag.ecdf<-ecdf(fiji$mag)
conf<-sqrt((1/(2*n))*log(2/.05))
pdf("hw3_3_b.pdf")
curve(mag.ecdf,xlim=c(3,7),n=1000,main="Empirical Cumulative Density of Magnitudes",xlab="Magnitude", ylab="Cumulative Density")
curve(pmin(mag.ecdf(x)+conf,1),xlim=c(3,7),n=1000,add=TRUE,lty=2)
curve(pmax(mag.ecdf(x)-conf,0),xlim=c(3,7),n=1000,add=TRUE,lty=2)
dev.off()

### c)
int.mag.den <- function(x) {
    mat<-matrix(rep(fiji$mag,length(x)),byrow=TRUE,nrow=length(x))
    out<-rowMeans(pnorm(x,mat,h))
    return(out)
}

int.est <- int.mag.den(4.9) - int.mag.den(4.3)
ecdf.est <- mag.ecdf(4.9) - mag.ecdf(4.3)

ub<-ecdf.est-pnorm(.025)*ecdf.est*(1-ecdf.est)/sqrt(n)
lb<-ecdf.est+pnorm(.025)*ecdf.est*(1-ecdf.est)/sqrt(n)

############################################
# Problem 4 
############################################

ru<-runif(50)
th<-max(ru)
n<-50
b<-1000

### a)
p.max.unif <- function(x) n*x^{n-1}

np.boots<-replicate(b,sample(ru,50,replace=TRUE))
np.th <- apply(np.boots,2,max)

p.boots<-replicate(b,runif(50,0,th))
p.th <- apply(p.boots,2,max)

lb<-.75
pdf("hw3_4_a_np.pdf")
hist(np.th,main=paste("Non Parametric Bootstrap, Theta Hat =",format(th,digits=3)),xlab="Theta Hat Star",ylab="Density",xlim=c(lb,1),probability=TRUE,breaks=seq(lb,1,.005),ylim=c(0,100))
curve(p.max.unif,xlim=c(lb,1),lty=2,add=TRUE)
legend('topleft', legend="Analytic Distribution", lty=2)
dev.off()

pdf("hw3_4_a_p.pdf")
hist(p.th,main=paste("Non Parametric Bootstrap, Theta Hat =",format(th,digits=3)),xlab="Theta Hat Star",ylab="Density",xlim=c(lb,1),probability=TRUE,breaks=seq(lb,1,.005),ylim=c(0,100))
curve(p.max.unif,xlim=c(lb,1),lty=2,add=TRUE)
legend('topleft', legend="Analytic Distribution", lty=2)
dev.off()

############################################
# Problem 5 
############################################

B<-1000
exp<- 500 
n<- 100

### a)
b0 <- -1
b1 <- 2
b2 <- -1

g<-function(x) b0 + b1*x + b2*x^2
X<-runif(n,0,2)
Y<-g(X)+rnorm(n,0,.2)

a.boot.theta<-function(c,nf,Xf,Yf) {
    sample<-sample(1:nf,nf,replace=TRUE)
    bh<-lm(I(Yf[sample])~I(Xf[sample])+I(Xf[sample]^2))$coef
    return(-bh[1]/(2*bh[2]))
}
a.boots<-mcmapply(a.boot.theta,c=1:B,MoreArgs=list(nf=n,Xf=X,Yf=Y),mc.cores=4)
a.ci<-quantile(a.boots,c(.025,.975))

Xt<-runif(exp,0,2)
Yt<-g(Xt)+rnorm(n,0,.2)
bht<-lm(Yt~Xt+I(Xt^2))$coef
a.true<- -bht[1]/(2*bht[2])

### b)
b.Z<-rnorm(n)
b.X<- 10*b.Z+rnorm(n)
b.Y<- 10*b.Z+rnorm(n)

b.matrix<-matrix(c(b.X,b.Y,b.Z),ncol=3)

b.theta <- function(M) {
    siginv<-solve(cov(M))
    return(-siginv[1,2]/sqrt(siginv[1,1]*siginv[2,2]))
}

b.boot.theta<-function(c,nf,M) {
    sample<-sample(1:nf,nf,replace=TRUE)
    return(b.theta(M[sample,]))
}
b.boots<-mcmapply(b.boot.theta,c=1:B,MoreArgs=list(nf=n,M=b.matrix),mc.cores=4)
b.ci<-quantile(b.boots,c(.025,.975))

Zt<-rnorm(exp)
b.true<-b.theta(matrix(c(10*Zt+rnorm(exp),10*Zt+rnorm(exp),Zt),ncol=3))

### c)
c.X<-mvrnorm(n=n,rep(0,10),diag(10))
c.theta <- function(M) min(eigen(cov(M))$values)
c.boot.theta <- function(c,nf,M) {
    sample<-sample(1:nf,nf,replace=TRUE)
    return(c.theta(M[sample,]))
}
c.boots<-mcmapply(c.boot.theta,c=1:B,MoreArgs=list(nf=n,M=c.X),mc.cores=4)
c.ci<-quantile(c.boots,c(.025,.975))
c.X.t<-mvrnorm(n=exp,rep(0,10),diag(10))
c.true<-c.theta(c.X.t)












    
