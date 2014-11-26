library('parallel')
library('alr3')
library('faraway')
library('car')
library('MASS')
library('locfit')
library('xtable')
library('nlme')
library('quantreg')
library('splines')
library('rpart')
library('randomForest')

#digits of precision
d<- 3

#####################
# Problem 1
#####################

data(aatemp)

### a)
temp.lm <- lm(temp~year,data=aatemp)
print(xtable(summary(temp.lm),digits=d),type="latex",file="hw7/hw7_1_a.tex")

pdf("hw7/hw7_1_a.pdf")
plot(aatemp,xlab="Year",ylab="Temperature",main="Linear Fit",xlim=c(1850,2000),pch=20)
abline(temp.lm,xlim=c(1850,2000))
dev.off()

### b)
temp.alm<-gls(temp~year,correlation=corAR1(form=~year),method="ML",data=aatemp)
print(xtable(summary(temp.alm)$tTable,digits=d),type="latex",file="hw7/hw7_1_b_ttable.tex")
print(xtable(intervals(temp.alm)$corStruct,digits=d),type="latex",file="hw7/hw7_1_b_acor.tex")

pdf("hw7/hw7_1_b.pdf")
matplot(aatemp$year, cbind(aatemp$temp,temp.alm$fit),type="pl",pch=20,lty=1,col=1,xlim=c(1850,2000),ylab="Temperature",xlab="Year",main="Autoregressive Fit")
dev.off()

### c)
crit<-.2
p<-1
deg<-11
while (p>crit) {
    deg <- deg-1
    temp.poly<-poly(aatemp$year,deg)
    poly.data<-data.frame(temp=aatemp$temp,temp.poly=temp.poly)
    temp.poly.lm <- lm(temp ~ .,data=poly.data)
    p<-coef(summary(temp.poly.lm))[deg+1,4]
}
print(xtable(summary(temp.poly.lm),digits=d),type="latex",file="hw7/hw7_1_c.tex")

pdf("hw7/hw7_1_c.pdf")
plot(aatemp,xlab="Year",ylab="Temperature",main="Polynomial Fit",xlim=c(1850,2000),pch=20)
curve(predict(temp.poly,x) %*% coef(temp.poly.lm)[-1] + coef(temp.poly.lm)[1],xlim=c(1850,2000),add=TRUE)
dev.off()

pred.data<-as.data.frame(predict(temp.poly,2020))
names(pred.data)<-paste("temp.poly.",1:d,sep="")
temp.pred<-predict(temp.poly.lm,newdata=data.frame(temp.poly=predict(temp.poly,2020)),interval="prediction")
print(xtable(temp.pred,digits=d),type="latex",file="hw7/hw7_1_c_pred.tex")

### d)
temp.ls.lm <- lm(temp~ifelse(year<1930,0,year-1930),data=aatemp)
print(xtable(summary(temp.ls.lm),digits=d),type="latex",file="hw7/hw7_1_d.tex")

pdf("hw7/hw7_1_d.pdf")
matplot(aatemp$year, cbind(aatemp$temp,temp.ls.lm$fit),type="pl",pch=20,lty=1,col=1,xlim=c(1850,2000),ylab="Temperature",xlab="Year",main="First Order Spline Fit")
dev.off()

### e)
y.ma <- max(aatemp$year)
y.mi <- min(aatemp$year)
spl <- splineDesign(c(rep(y.mi,3),seq(y.mi,y.ma,length.out=5),rep(y.ma,3)),aatemp$year)
temp.cs.lm <- lm(aatemp$temp ~ spl)
print(xtable(summary(temp.cs.lm),digits=d),type="latex",file="hw7/hw7_1_e.tex")

pdf("hw7/hw7_1_e.pdf")
matplot(aatemp$year, cbind(aatemp$temp,temp.cs.lm$fit),type="pl",pch=20,lty=1,col=1,xlim=c(1850,2000),ylab="Temperature",xlab="Year",main="Cubic Spline Fit")
dev.off()

#####################
# Problem 2 
#####################
data(ozone)

### a)
oz.lm<-lm(O3~temp+humidity+ibh, data=ozone)
print(xtable(summary(oz.lm),digits=d),type="latex",file="hw7/hw7_2_a.tex")

pdf("hw7/hw7_2_a.pdf")
boxcox(oz.lm,plotit=TRUE, lambda=seq(0,.5,.1))
title(main="Box-Cox")
dev.off()

lam<-.25
oz.t.lm<-lm(O3^lam ~temp+humidity+ibh, data=ozone)
print(xtable(summary(oz.t.lm),digits=d),type="latex",file="hw7/hw7_2_a_trans.tex")

### b)
oz.t.rt<-rpart(O3^lam ~temp+humidity+ibh, data=ozone)
oz.t.rf<-randomForest(O3^lam ~temp+humidity+ibh, data=ozone)

pdf("hw7/hw7_2_b_tree.pdf")
plot(oz.t.rt,main="Regression Tree",uniform=TRUE)
text(oz.t.rt)
dev.off()

#Linear diagnoistics
pre<-predict(oz.t.lm)
res<-residuals(oz.t.lm)
pdf("hw7/hw7_2_b_lm_res.pdf")
plot(pre,res,xlab="Predicted",ylab="Residual",main="Linear Model Residual Plot")
lines(predict(loess(res~pre)))
legend('topright',"LOESS Fit",lty=1)
dev.off()
pdf("hw7/hw7_2_b_lm_qq.pdf")
qqnorm(res,ylab="Linear Model Residual Quantiles",xlab="Normal Theoretical Quantiles",main="Linear Model Residual QQ Plot")
qqline(res)
dev.off()

#regression tree diagnoistics
pre<-predict(oz.t.rt)
res<-residuals(oz.t.rt)
pdf("hw7/hw7_2_b_rt_res.pdf")
plot(pre,res,xlab="Predicted",ylab="Residual",main="Regression Tree Residual Plot")
lines(predict(loess(res~pre)))
legend('topright',"LOESS Fit",lty=1)
dev.off()
pdf("hw7/hw7_2_b_rt_qq.pdf")
qqnorm(res,ylab="Regression Tree Residual Quantiles",xlab="Normal Theoretical Quantiles",main="Regression Tree Residual QQ Plot")
qqline(res)
dev.off()

# Random Forrest diagnoistics
pre<-predict(oz.t.rf)
res<-ozone$O3^lam-predict(oz.t.rf)
pdf("hw7/hw7_2_b_rf_res.pdf")
plot(pre,res,xlab="Predicted",ylab="Residual",main="Random Forrest Residual Plot")
lines(predict(loess(res~pre)))
legend('topright',"LOESS Fit",lty=1)
dev.off()
pdf("hw7/hw7_2_b_rf_qq.pdf")
qqnorm(res,ylab="Random Forrest Residual Quantiles",xlab="Normal Theoretical Quantiles",main="Random Forrest Residual QQ Plot")
qqline(res)
dev.off()

#####################
# Problem 3 
#####################
data(prostate)
crit<-.2

### a)
exclude<-NULL
excl.pval<-NULL
go<-TRUE
p.lm<-lm(lpsa~.,data=prostate)

while (go) {
    p.lm<-lm(lpsa~.,data=prostate[!(names(prostate) %in% exclude)])
    pvals<-summary(p.lm)$coefficients[-1,4]
    max.ind <- which.max(pvals)
    max.pval <- pvals[max.ind]
    if (max.pval>crit) {
        excl.pval<-c(excl.pval,max.pval)
        exclude<-c(exclude,names(max.ind))
    }
    else { 
        go<-FALSE
    }
}
ecl.ma<-cbind(1:length(excl.pval),excl.pval)
colnames(ecl.ma)<-c("Order","p-value")
print(xtable(summary(p.lm),digits=d),type="latex",file="hw7/hw7_3_a_lm.tex")
print(xtable(ecl.ma,digits=d),type="latex",file="hw7/hw7_3_a_excl.tex")

### b)

crit<-.2
include<-NULL
incl.pvals<-NULL
remaining<-1:8
go<-TRUE

while (go) {
    mpval<-1
    for (nv in remaining) {
        p.lm<-lm(lpsa~.,data=prostate[c(include,nv,9)])
        pval<-summary(p.lm)$coefficients[length(include)+2,4]
        if (mpval>pval) {
            bv   <-nv
            mpval<-pval
        }
    }
    if (mpval<crit) {
        include<-c(include,bv)
        incl.pvals<-c(incl.pvals,mpval)
        remaining<-remaining[-which(remaining==bv)]
    }
    else {
        go<-FALSE
    }
}
names(pvals)<-names(prostate)[include]
incl.ma <- cbind(1:length(incl.pvals),incl.pvals)
rownames(incl.ma)<-names(prostate)[include]
colnames(incl.ma)<-c("Order","p-value")

p.lm<-lm(lpsa~.,data=prostate[c(include,9)])

print(xtable(summary(p.lm),digits=d),type="latex",file="hw7/hw7_3_b_lm.tex")
print(xtable(incl.ma,digits=d),type="latex",file="hw7/hw7_3_b_incl.tex")

       












