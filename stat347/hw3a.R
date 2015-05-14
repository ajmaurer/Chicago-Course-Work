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

