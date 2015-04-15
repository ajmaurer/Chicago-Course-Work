library('parallel')
library('alr3')
library('faraway')
library('car')
library('MASS')
library('locfit')
library('xtable')

#############################
# Problem 1
#############################

data('sat')
attach(sat)
lm.sat<-lm(total~expend + salary + ratio + takers)
n<-nrow(sat)
names<-row.names(sat)

### a)
pdf("hw5_1_a_res.pdf")
plot(fitted(lm.sat),residuals(lm.sat),xlab="Fitted SAT Score",ylab="Residuals",main="Residuals Plot")
abline(h=0)
legend('topleft', legend=c('y=0'), lty=c(1))
dev.off()

lm.res<-lm(sqrt(abs(residuals(lm.sat)))~predict(lm.sat))
pdf("hw5_1_a_sqrtres.pdf")
plot(fitted(lm.sat),sqrt(abs(residuals(lm.sat))),xlab="Fitted SAT Score",ylab="Square Root of Absolute Value of Residuals",main="Square Root Residuals Plot")
abline(lm.res)
legend('topleft', legend=c("Linear Fit"), lty=c(1))
dev.off()

### b)
pdf("hw5_1_b_qq.pdf")
qqnorm(residuals(lm.sat),ylab="Residuals",xlab="Theoretical Normal Quantiles", main="Residuals vs. Normal Qunatile Plot")
qqline(residuals(lm.sat))
dev.off()
shap<-shapiro.test(residuals(lm.sat))

### c)
pdf("hw5_1_c_halfnorm.pdf")
halfnorm(hatvalues(lm.sat),labs=names,ylab="Leverages",main="Leverages vs. Half-Normal Quantiles")
dev.off()

### d)
pdf("hw5_1_d_tqq.pdf")
qqplot(rt(1000,df=n-5),rstudent(lm.sat),xlab="T Distribution With 45 Degrees of Freedom", ylab="Studentized Residuals", main="QQ Plot Studentized Resisduals vs. T Distribution")
abline(0,1)
dev.off()
stud.max<-rstudent(lm.sat)[which.max(abs(rstudent(lm.sat)))]
qtpt10<-qt(.1/100,45)
qtpt05<-qt(.05/100,45)

### e)
pdf("hw5_1_e_halfnorm.pdf")
halfnorm(cooks.distance(lm.sat),3,labs=names,ylab="Cook's Distance",main="Cook's Distance vs. Half-Normal Quantiles")
dev.off()
vars<-c("Intercept",attr(terms(lm.sat),"term.labels"))
for (i in 1:5) {
    minmax<-c(min(dfbeta(lm.sat)[,i]),max(dfbeta(lm.sat)[,i]))
    cutoff<-sort(abs(dfbeta(lm.sat)[,i]))[n-2]
    extreme<-abs(dfbeta(lm.sat)[,i])>=cutoff
    pdf(paste("hw5_1_e_",vars[i],".pdf",sep=""))
    plot(dfbeta(lm.sat)[!extreme,i],ylab=paste("Change in",vars[i],"coef"),xlab="Index",main=paste("Influence of Observations on",vars[i],"coef"),xlim=c(0,n),ylim=minmax)
    text((1:n)[extreme],dfbeta(lm.sat)[extreme,i],names[extreme])
    dev.off()
}

### d)
lf.res<-locfit(residuals(lm.sat)~lp(predict(lm.sat),deg=1))
lm.res<-lm(residuals(lm.sat)~predict(lm.sat))
pdf("hw5_1_f_locfitres.pdf")
plot(predict(lm.sat),residuals(lm.sat),main="Residuals vs. Predicted",ylab="Residuals",xlab="Predicted Values")
abline(lm.res,lty=1)
lines(lf.res,lty=2)
legend('bottomleft', legend=c('Linear Fit','Local Linear Fit'), lty=c(1,2))
dev.off()
for (i in 2:5) {
    x<-model.matrix(lm.sat)[,i]
    partial.res<-residuals(lm.sat)+lm.sat$coef[i]*(x-mean(x))
    pres.lf<-locfit(partial.res~lp(x,deg=1))
    pres.lm<-lm(partial.res~x)
    pdf(paste("hw5_1_f_",vars[i],".pdf",sep=""))
    plot(x,partial.res,xlab=vars[i],ylab=paste(vars[i],"Partial Residual"),main=paste(vars[i],"Partial Residual Plot"))
    abline(pres.lm,lty=1)
    lines(pres.lf,lty=2)
    legend('topright', legend=c('Linear Fit','Local Linear Fit'), lty=c(1,2))
    dev.off()
}

#############################
# Problem 2 
#############################

data('faithful')
attach(faithful)
n<-nrow(faithful)

### a)
pdf("hw5_2_a.pdf")
lm.eru<-lm(eruptions~waiting)
plot(waiting,eruptions,xlab="Waiting Time, Min",ylab="Eruption Time, Min",main="Waiting vs. Eruption Time")
abline(lm.eru)
dev.off()

### b)
var<-rep(seq(1,5,1),each=1000)
testcoef<-function(e) lm(eruptions~I(waiting+sample((-e):e,size=n,replace=TRUE)))$coef[[2]]
slopes<-unlist(mclapply(var,testcoef,mc.cores=4))
betas<-c(lm.eru$coef[[2]],tapply(slopes,var,mean))
var<-((.5+seq(0,5,1))^2)/3
lm.beta<-lm(betas~var)
pdf("hw5_2_b.pdf")
plot(var,betas,xlab="Assumed Variance",ylab="Beta",main="Model Beta Versus Variance In Preidctor",xlim=c(0,max(var)),ylim=c(min(betas),lm.beta$coef[1]))
abline(lm.beta)
dev.off()

### d)
pdf("hw5_2_d.pdf")
plot(predict(lm.eru),residuals(lm.eru),xlab="Predicted Values",ylab="Residuals",main="Residual Plot")
abline(h=0)
dev.off()
higher<-waiting>68
waiting.low<-ifelse(higher,0,waiting)
waiting.high<-ifelse(higher,waiting,0)
lm.eru2<-lm(eruptions~higher + waiting.low + waiting.high)
ana<-anova(lm.eru,lm.eru2)
print(xtable(ana),type="latex",file="hw5_2_d.tex")

#############################
# Problem 3 
#############################
data('savings')
data('star')

b<-1000

lm.sav<-lm(sr~pop15+pop75+dpi+ddpi,savings)
lm.star<-lm(light ~ temp, data=star)

pred.sav<-predict(lm.sav)
pred.star<-predict(lm.star)

sig.sav<-summary(lm.sav)$sigma
sig.star<-summary(lm.star)$sigma

n.sav<-nrow(savings)
n.star<-nrow(star)

### a)
m.rs.sav<-max(abs(rstudent(lm.sav)))
m.rs.star<-max(abs(rstudent(lm.star)))

sav.rs.boot<-function(a) max(abs(rstudent(lm(I(pred.sav+rnorm(n.sav,sig.sav))~pop15+pop75+dpi+ddpi,savings))))
star.rs.boot<-function(a) max(abs(rstudent(lm(I(pred.star+rnorm(n.star,sig.star))~temp,star))))

sav.rs.boots<-mclapply(seq(1,b),sav.rs.boot,mc.cores=4)
star.rs.boots<-mclapply(seq(1,b),star.rs.boot,mc.cores=4)

p.rs.sav<-sum(sav.rs.boots>m.rs.sav)/b
p.rs.star<-sum(star.rs.boots>m.rs.star)/b

### b)
m.ck.sav<-max(abs(cooks.distance(lm.sav)))
m.ck.star<-max(abs(cooks.distance(lm.star)))

sav.ck.boot<-function(a) max(abs(cooks.distance(lm(I(pred.sav+rnorm(n.sav,sig.sav))~pop15+pop75+dpi+ddpi,savings))))
star.ck.boot<-function(a) max(abs(cooks.distance(lm(I(pred.star+rnorm(n.star,sig.star))~temp,star))))

sav.ck.boots<-mclapply(seq(1,b),sav.ck.boot,mc.cores=4)
star.ck.boots<-mclapply(seq(1,b),star.ck.boot,mc.cores=4)

p.ck.sav<-sum(sav.ck.boots>m.ck.sav)/b
p.ck.star<-sum(star.ck.boots>m.ck.star)/b

##########################
# Problem 5
##########################
data('longley')
attach(longley)
lm.long<-lm(Employed~GNP.deflator+GNP+Unemployed+Armed.Forces+Population+Year)
x<-model.matrix(lm.long)[,-1]

### a)
e<-eigen(t(x)%*%x)$val
cond<-sqrt(e/min(e[e>0]))
print(xtable(t(as.matrix(cond))),type="latex",file="hw5_6_a.tex")

### b)
cor<-round(cor(longley[,-7]),3)
print(xtable(cor),type="latex",file="hw5_6_b.tex")

### c)
vif<-vif(x)
print(xtable(t(as.matrix(vif))),type="latex",file="hw5_6_c.tex")



