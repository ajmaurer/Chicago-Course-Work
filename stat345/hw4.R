library('xtable')
library('stats4')
library('locfit')
library('coefplot')
library('SDMTools')
library('plotrix')
library('combinat')
library('reshape2')
library('nlme')

### Problem 2
rye<-read.table("hw4/ryegrass.dat",sep="",header=T)
rye$Block<-as.factor(rye$Block)
rye$str.man<-rep(1:8,each=4)
rye.rs<-dcast(melt(rye,id=c("Strain","Manure","Block","str.man"),variable.name=c("Yield"),measure.vars="Yield"),Strain+Manure+str.man~Block+Yield)

# a)
pdf('hw4/hw4_2a_plot.pdf')
matplot(x=rye.rs$str.man[rye.rs$Manure=="A"]-.2,y=rye.rs[rye.rs$Manure=="A",4],axes=F,pch="1",cex=.8,ylim=c(0,450),xlim=c(.9,8.1),col=1,main="Yield by Strain/Manure/Block",xlab="Strain & Manure Level",ylab="Yield")
matplot(x=rye.rs$str.man[rye.rs$Manure=="H"]-.2,y=rye.rs[rye.rs$Manure=="H",4],cex=.8,col=2,add=T,pch="1")
for (i in 1:3) {
    for (l in c("A","H")){
        c<- ifelse(l=="A",1,2)
        matplot(x=rye.rs$str.man[rye.rs$Manure==l]-.2+.08*i,y=rye.rs[rye.rs$Manure==l,4+i],cex=.8,col=c,add=T,pch=paste(i+1))
    }
}
mtext("Digit Indicates Block; Red is High Manure, Black is Average Manure")
axis(1,labels=paste(rep("",8)),at=1:8,pos=0,lwd=0,lwd.ticks=1)
axis(1,labels=paste(c("S23","","NZ","","X","","Kent",""),rep(c("- High",""),4)),at=1:8,pos=5,lwd=0,lwd.ticks=0)
axis(1,labels=paste(c("","S23","","NZ","","X","","Kent"),rep(c("","- Avg"),4)),at=1:8,pos=-12.5,lwd=0,lwd.ticks=0)
axis(1,labels=c("",""),at=c(0,8.5),pos=0,lwd=1,lwd.ticks=0)
axis(2,at=seq(0,450,by=50))
dev.off()

# b)
rye.lme <- lme(Yield~ Strain*Manure, random= ~1 | Block, data=rye)
f.coef<-as.matrix(summary(rye.lme)$coefficients$fixed)
colnames(f.coef)<-"Fixed Effect"
r.coef<-as.matrix(summary(rye.lme)$coefficients$random$Block)
colnames(r.coef)<-"Random Effect"
print(xtable(f.coef,digits=2),type="latex",file="hw4/hw4_2b_lme_fixed.tex")
print(xtable(r.coef,digits=2),type="latex",file="hw4/hw4_2b_lme_random.tex")

# c)

c<-cbind(c(0,0,0,0,1,0,0,0),c(0,0,0,0,1,1,0,0),c(0,0,0,0,1,0,1,0),c(0,0,0,0,1,0,0,1))

coefficients<-NULL
resid<-fitted(rye.lme)-rye$Yield
res.blk <- summary(rye.lme)$coefficients$random$Block
for (i in 1:1000) {
    sim.yield<- summary(rye.lme)$fitted[,"fixed"] + sample(resid,32,replace=T) + rep(sample(res.blk,4,replace=T),8)
    #sim.yield<- summary(rye.lme)$fitted[,"fixed"] + rnorm(32,0,28.52) + rep(rnorm(4,0,21.8314),8)
    coefficients<-cbind(coefficients,t(c)%*%as.matrix(summary(lme(sim.yield~ Strain*Manure, random= ~1 | Block, data=rye))$coefficients$fixed))
}

pred<-t(c)%*%f.coef

# Parametric
cov<-summary(rye.lme)$varFix
par.std <- sqrt(diag(t(c)%*% cov %*% c))
par.p<- 2*pt(-abs(pred/par.std),21)
par<-cbind(pred,par.std,par.p)
colnames(par)<-c("Coefficient","SE","P-Value")
rownames(par)<-c("Kent","NZ","S23","X")

f.mat <- matrix(rep(pred,1000),ncol=1000)
bs.pval<-apply(abs(coefficients-f.mat)>f.mat,1,sum)/1000
bs.std<-apply(coefficients,1,sd)
bs<-cbind(pred,bs.std,bs.pval)
colnames(bs)<-c("Coefficient","SE","P-Value")
rownames(bs)<-c("Kent","NZ","S23","X")

print(xtable(par,digits=2),type="latex",file="hw4/hw4_2c_lme_par_ttable.tex")
print(xtable(bs,digits=2),type="latex",file="hw4/hw4_2c_lme_bs_ttable.tex")

# d) 
d<-as.matrix(c(0,-1,1,0,0,-1,1,0))
est <- t(d) %*% f.coef
var <- t(d) %*% cov %*% d
ci = est + c(-1,1)*sqrt(var)*pt(.025,21)



