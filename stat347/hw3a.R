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
library('lmtest')
library('Hmisc')

### Problem 0
comm2   <- read.table('hw3/commensal2.txt',header=T)[1:15,]
comm2.col <- data.frame(Generation=as.numeric(comm2[,1]),homo=comm2[,2]+comm2[,5],hetero=comm2[,3]+comm2[,4])

comm.glm <- glm(cbind(homo,hetero)~poly(Generation,2),family=binomial,data=comm2.col)
comm.glm.sum <-summary(comm.glm)
print(xtable(comm.glm.sum),type="latex",file="hw3/A0_ix_logit.tex")

comm.cs <- chisq.test(comm2[1:15,2:5])
comm.col.cs <- chisq.test(comm2.col[1:15,2:3])

### Problem 3 
byss <- read.table('hw3/byss.txt',header=T)
byss.col <- colnames(byss)
for (c in 3:7) {
    byss[,byss.col[c]] <- as.factor(byss[,byss.col[c]])
}
## (a)
byss.glm <- glm(cbind(affected,not_affected)~race+sex+smok+empl+dust,family=binomial,data=byss)
print(xtable(summary(byss.glm)),type="latex",file="hw3/A3_a_logit.tex")

## (b) estimate odds ratio CI for sex2
sex2.ci <- .1239+c(-1,1)*qnorm(.05)*.2288
sex2.od.ci <- exp(sex2.ci)

## (c) drop sex then race, in that order
byss.glm2 <- glm(cbind(affected,not_affected)~smok+empl+dust,family=binomial,data=byss)
print(xtable(summary(byss.glm2)),type="latex",file="hw3/A3_c_logit.tex")

## (d) largest interaction - turns out to be sex*dust
byss.glm.int <- glm(cbind(affected,not_affected)~race+sex+smok+empl+dust+sex:dust,family=binomial,data=byss)
byss.lrtest <- lrtest(byss.glm.int,byss.glm)
byss.glm.int2 <- glm(cbind(affected,not_affected)~sex+smok+empl+dust+sex:dust,family=binomial,data=byss)
print(xtable(summary(byss.glm.int2)),type="latex",file="hw3/A3_d_logit.tex")

## (e) summary
beta <- coef(byss.glm.int2)
sig <- vcov(byss.glm.int2)
lincombo <- matrix(,nrow=5,ncol=9)
lincombo[1,] <- c(0,0,0,0,0,1,0,0,0) # man, dust 2
lincombo[2,] <- c(0,0,0,0,0,0,1,0,0) # man, dust 3 
lincombo[3,] <- c(0,1,0,0,0,0,0,0,0) # woman, dust 1
lincombo[4,] <- c(0,1,0,0,0,1,0,1,0) # woman, dust 2
lincombo[5,] <- c(0,1,0,0,0,0,1,0,1) # woman, dust 3 

# estimates and standard errors for given subgroups
est <- lincombo %*% beta
std <- sqrt(diag(lincombo %*% sig %*% t(lincombo) ))
odds_est <-exp(cbind(est,est-1.96*std,est+1.96*std))
odds_est <-cbind(rbind(c(1,1,1),odds_est[1:2,]),odds_est[3:5,])

pdf('hw3/A3_e_plot.pdf')
matplot(1:3,odds_est,type='bbbbbb',pch=c(20,2,6,20,2,6),col=c(4,4,4,2,2,2),xaxt="n",yaxt="n",main="Odds Ratio of Byssinosis Compared to Low Dust Men",ylab="Odds Ratio",xlab="Dust Level",lty=c(1,3,3,1,3,3))
abline(a=1,b=0,lty=2)
axis(1,at=1:3,labels=c('Low','Medium','High'))
axis(2)
minor.tick(ny=5,nx=0)
legend('topleft',legend=c('Men','Women','Estimate','95% Upper Bound','95% Lower Bound'),lty=c(1,1,1,3,3),col=c(4,2,1,1,1),pch=c(NA,NA,20,6,2))
dev.off()


