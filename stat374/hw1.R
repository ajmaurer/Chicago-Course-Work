library('parallel','locfit')

# Option on which section of code to run

############################################
# Problem 1
############################################

#Part a
m<-1
s<-3
b<-100
nseq<-1:200
uhat<-function(n) mean(rnorm(n,m,s))
meansqerr<-function(n) mean((m-sapply(rep(n,times=b),uhat))^2)
empmsqerr<-sapply(nseq, meansqerr)
anamsqerr<-sapply(nseq, function(n) (s^2)/n)

pdf('hw1_1_a_reg.pdf')
plot(1:200,empmsqerr,main="Empirical vs Analytic Mean Squared Error", pch=1, cex=.5, xlab="N", ylab="Mean Squarred Error")
points(1:200,anamsqerr,cex=.5, pch=18)
legend('topright', legend=c('Empirical','Analytic'), pch=c(1,18))
dev.off()

pdf('hw1_1_a_log.pdf')
plot(1:200,empmsqerr,main="Empirical vs Analytic Mean Squared Error", pch=1, cex=.5, xlab="Log N", ylab="Log Mean Squarred Error", log='xy')
points(1:200,anamsqerr,cex=.5, pch=18)
legend('topright', legend=c('Empirical','Analytic'), pch=c(1,18))
dev.off()

#Part b 
nseq<-c(10,50,100)
genZ<-function(n) sqrt(n)*(m-sapply(rep(n,times=b),uhat))
for (n in nseq) {
    pdf(paste('hw1_1_b_',toString(n),'.pdf', sep=""))
    plot(density(genZ(n)), main="Empirical vs Analytic Density of Z", sub=paste('N=',toString(n),sep=""), lty=2 , xlab="Z", ylab="Density", xlim=c(-10,10), ylim=c(0,.25))
    curve(dnorm(x,0,s), add=TRUE, lty=1)
    legend('topright', legend=c('Empirical','Analytic'), lty=c(2,1))
    dev.off()
}


############################################
# Problem 2
############################################

### Sigma^2 = 1

#plot the risk in terms of b
pdf('hw1_2_plot_1.pdf')
curve(1000*x^2+(1-x)^2*1.08232,xlim=c(-.2,1.2),, ylim=c(0,2), xlab='b', ylab='Risk', main='Risk As a Function of b')
dev.off()

#simulate the risk
n<-1000
b<-100
s<-1
th.i<- function(i) rnorm(1,1/i^2,s)
JS<- function(n) max(1- n*s^2/sum(sapply(1:n,th.i)^2),0)
JSrisk<- function(n) {
    norm<-rnorm(n,0,s)
    theta<-1/((1:n)^2)
    z<-norm+theta
    js<-1- n*s^2/sum(z^2)
    return((js*z-theta)^2)
}
simrisk<-sapply(rep(n,times=b),JSrisk)
pdf('hw1_2_risk_1.pdf')
hist(simrisk, xlab='Simulated James-Stein Estimator Mean Squarred Error', ylab='Frequency', main=paste('Performance James-Stein Estimator, Mean=',toString(format(mean(simrisk),trim=TRUE))))
dev.off()
sim<-sapply(rep(n,times=b),JS)
mean<-format(mean(sim),trim=TRUE)
pdf('hw1_2_sim_1.pdf')
hist(sim, xlab='Simulated James-Stein Estimator', ylab='Frequency', main=paste('Performance James-Stein Estimator, Mean=',toString(mean)))
dev.off()


### sigma^2 = 1/n

#plot the risk in terms of b
pdf('hw1_2_plot_n.pdf')
curve(x^2+(1-x)^2*1.08232,xlim=c(-.2,1.2),, ylim=c(0,2), xlab='b', ylab='Risk', main='Risk As a Function of b')
dev.off()

#simulate the risk
n<-1000
b<-100
s<-sqrt(1/n)
th.i<- function(i) rnorm(1,1/i^2,s)
JS<- function(n) 1- n*s^2/sum(sapply(1:n,th.i)^2)
JSrisk<- function(n) {
    norm<-rnorm(n,0,s)
    theta<-1/((1:n)^2)
    z<-norm+theta
    js<-1- n*s^2/sum(z^2)
    return((js*z-theta)^2)
}
simrisk<-sapply(rep(n,times=b),JSrisk)
pdf('hw1_2_risk_n.pdf')
hist(simrisk, xlab='Simulated James-Stein Estimator Mean Squarred Error', ylab='Frequency', main=paste('Performance James-Stein Estimator, Mean=',toString(format(mean(simrisk),trim=TRUE))))
dev.off()
sim<-sapply(rep(n,times=b),JS)
mean<-format(mean(sim),trim=TRUE)
pdf('hw1_2_sim_n.pdf')
hist(sim, xlab='Simulated James-Stein Estimator', ylab='Frequency', main=paste('Performance James-Stein Estimator, Mean=',toString(mean)))
dev.off()


############################################
# Problem 3
############################################

z<-read.csv(file="assn1-prob3-data.txt",head=FALSE)$V1

# This data set is known to be 1000 draws of N(0,1) plus value.
# My suspicion is that the means are draws from a pareto distribution with min ~2.5 and shape ~2.25.
# define function
rpareto<-function(n,a,m) m/(runif(n)^(1/a)) #random pareto vector
dpareto<-function(x,a,m) a*m^a/(x^(a+1))            #pareto density
den<-function(x,z,a,m) dpareto(x,a,m)*dnorm(z-x) #conditional density given z

# parameters
n<-1000
b<-25
al<-2
min<-2.4

# Simulate
p<-rnorm(1000)+rpareto(1000,al,min)
summary(p)
summary(z)

# QQ plot two simulated
pdf('hw1_3_qq.pdf')
qqplot(p,z, xlim=c(0,max(p,z)), main= 'Z versus Simulated', ylim=c(0,max(p,z)), xlab='Simulated Pareto + N(0,1)', ylab='Z')
curve(x+0, add=TRUE, lty=1)
dev.off()

# Emperical versus simulated density
pdf('hw1_3_den.pdf')
plot(density(z), lty=2, ylab='Density', xlab='Z', main='Density of Z', xlim=c(-.5,75))
condden<-function(zp) integrate(function(p) den(p,zp,al,min), lower=min, upper=Inf, rel.tol=.Machine$double.eps^.85)$value
curve(sapply(x, condden), add=TRUE)
legend('topright', legend=c('Kernel Density Estimate of Z','Analytic Pareto + N(0,1)'), lty=c(2,1))
dev.off()

# Maximize each estimate th
th<-rep(0, times=1000)
for (i in 1:1000) {
    scale=1/integrate(function(x) den(x,z[i],al,min),lower=min,upper=Inf, rel.tol=.Machine$double.eps^.85)$value
    minfunc<-function(b) scale*integrate(function(x) (x-b)^2 * den(x,z[i],al,min),lower=min,upper=Inf, rel.tol=.Machine$double.eps^.85)$value
    th[i]<-optimize(f=minfunc,lower=0,c(z[i]/2,100))$minimum
}
th[z>30]<-(z-1/z)[z>30] # We run into numerical issues with large y, so we use stein

#Plot predicted theta versus Z
pdf('hw1_3_prediction.pdf')
plot(z,th, xlab='Original Z', ylab='Predicted Theta', main='Z vs. Predicted Thetas', xlim=c(-.5,75), ylim=c(-.5,75))
curve(x+0, add=TRUE, lty=1)
dev.off()

# Output 
write.table(th,file="assn1-prob3-ajmaurer.txt", row.names=FALSE, col.names=FALSE)


############################################
# Problem 5
############################################

doppler<-function(x) sqrt(x*(1-x))*sin(2.1*pi/(x+.05))
n<-1000
x<-(1:n)/n

# sigma = .1
y.1 <-sapply(x,doppler)+.1*rnorm(n)
gcv.1<-gcvplot(y.1~x,alpha=as.matrix((5:250)/1000))
pdf('hw1_5_p1_bandwith.pdf')
plot(gcv.1$alpha,gcv.1$values, main='CV Score, Sigma=.1', xlab='Bandwith', ylab='CV Score')
dev.off()

band.1 <- (gcv.1$alpha[order(gcv.1$values)])[1]
loc.1<-locfit(y.1~x,alpha=band.1)
crit(loc.1)<-crit(loc.1,cov=.95)
pdf('hw1_5_p1_fit.pdf')
plot(loc.1,band='local',main=paste('Doppler With Sigma=.1, Fit With Bandwith=',toString(band.1)),ylab='Y',xlab='x')
points(x,y.1, cex=.5)
dev.off()

# sigma = 1
y1 <-sapply(x,doppler)+rnorm(n)
gcv1<-gcvplot(y1~x,alpha=as.matrix((50:500)/1000))
pdf('hw1_5_1_bandwith.pdf')
plot(gcv1$alpha,gcv1$values, main='CV Score, Sigma=1', xlab='Bandwith', ylab='CV Score')
dev.off()

band1 <- (gcv1$alpha[order(gcv1$values)])[1]
loc1<-locfit(y1~x,alpha=band1)
crit(loc1)<-crit(loc1,cov=.95)
pdf('hw1_5_1_fit.pdf')
plot(loc1,band='local',main=paste('Doppler With Sigma=1, Fit With Bandwith=',toString(band1)),ylab='Y',xlab='x')
points(x,y1, cex=.5)
dev.off()







    

