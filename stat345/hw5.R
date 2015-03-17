library('xtable')
library('stats4')
library('locfit')
library('coefplot')
library('SDMTools')
library('plotrix')
library('combinat')
library('reshape2')
library('nlme')
library('parallel')
library('faraway')

######## Problem 1
# 2^4 factorial design, 3,4x interactions 0, $Z$ is distributed the same as the estimate of variance gotten after
# testing for the presence of 2x interactions (which are the betas)
Zgen<-function(i,betas) {
    Sy<-sum(rnorm(5,0,1)^2)/5
    Sx<-sum(rnorm(6,2*betas,1)^2)/6
    return(ifelse(Sx/Sy>4.95,Sy,(5*Sy+6*Sx)/11)) 
}

### b)
boots.1.b<-1000
betas.1.b<-rep(0,6)
Z.1.b<-mcmapply(Zgen,i=1:boots.1.b,MoreArgs=list(betas<-betas.1.b),mc.cores=4)
meanZ.1.b<-mean(Z.1.b)
mseZ.1.b<-mean((Z.1.b-1)^2)

### c)
boots.1.c<-1000
meanZ.1.c<-NULL
mseZ.1.c<-NULL
for (b in 1:5) {
    betas<-c(2^(b-2),rep(0,5))
    Z<-mcmapply(Zgen,i=1:boots.1.c,MoreArgs=list(betas<-betas),mc.cores=4)
    meanZ.1.c[b]<-mean(Z)
    mseZ.1.c[b]<-mean((Z-1)^2)
}

table.1.c<-rbind(meanZ.1.c,mseZ.1.c)
rownames(table.1.c)<-c("Z Mean","MSE")
colnames(table.1.c)<-2^(-1:3)

print(xtable(table.1.c,digits=3),type="latex",file="hw5/hw5_1c.tex")

######## Problem 2
shr<-read.table("hw5/shrinkage.dat",sep="",header=T)

shr2<-shr[1:4]
for (i in 1:4) shr2[i]<-ifelse(shr[i]=="+",1,-1)
shr2["E"]<-shr2[1]*shr2[2]*shr2[3]
shr2["F"]<-shr2[4]*shr2[2]*shr2[3]
shr2["shr"]<-shr$Shrinkage

options(contrasts=c("contr.sum","contr.poly"))

### b)
coef.2.b.s <- coef(lm(Shrinkage~A*B*C*D,data=shr))
pdf("hw5/hw5_2b_halfnorm.pdf")
halfnorm(coef.2.b.s,5,names(coef.2.b.s),xlim=c(0,2),ylab="Coefficients",main="Coefficient Half-Normal Plot")
dev.off()

print(xtable(summary(lm(Shrinkage~A*B*C*D,data=shr))$coef,digits=3),type="latex",file="hw5/hw5_2b_rtable.tex")

res.2.b<-residuals(lm(shr~A+B+C+D+E+F,data=shr2))
print(xtable(summary(lm(log(abs(res.2.b))~A+B+C+D+E+F,data=shr2))$coef,digits=3),type="latex",file="hw5/hw5_2b_lrtable.tex")

### c)
coef.2.c.s <- coef(lm(Shrinkage~A*C*D,data=shr,subset=(B=="-")))
pdf("hw5/hw5_2c_halfnorm.pdf")
halfnorm(coef.2.c.s,5,names(coef.2.c.s),xlim=c(0,2),ylab="Coefficients",main="Coefficient Half-Normal Plot")
dev.off()

print(xtable(summary(lm(Shrinkage~A*C*D,data=shr,subset=(B=="-")))$coef,digits=3),type="latex",file="hw5/hw5_2c_rtable.tex")









