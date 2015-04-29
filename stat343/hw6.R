library('parallel')
library('alr3')
library('faraway')
library('car')
library('MASS')
library('locfit')
library('xtable')
library('nlme')
library('quantreg')

#digits of precision
d<- 3

#############################
# Problem 1
#############################

data(ufcgf)
attach(ufcgf)

### 5.3.1.1
fir.lm<-lm(Height~Dbh)
fir.loess<-loess(Height~Dbh)

G<-sum((predict(fir.lm)-predict(fir.loess))^2)/(summary(fir.lm)$sigma^2)

### 5.3.1.2-4
B<-1000
Yhat<-predict(fir.lm)
ehat<-residuals(fir.lm)

res.boot<- function(c,x,y,e) {
    ystar<-y+sample(e,length(y),replace=TRUE)
    lm.m    <-lm(ystar~x)
    loess.m <-loess(ystar~x)
    Gstar<-sum((predict(lm.m)-predict(loess.m))^2)/(summary(lm.m)$sigma^2)
    return(Gstar)
}
res.boots<-mcmapply(res.boot,c=1:B,MoreArgs=list(x=Dbh,y=Yhat,e=ehat),mc.cores=4)
fir.p<-sum(res.boots>G)/B

#############################
# Problem 3 
#############################

data(divusa)
attach(divusa)
n<-nrow(divusa)
div.lm<-lm(divorce~unemployed+femlab+marriage+birth+military)

### a)
pdf("hw6_3_a_res.pdf")
plot(year,residuals(div.lm),ylab="Residual",xlab="Year",main="Residual by Year")
dev.off()

pdf("hw6_3_a_auto.pdf")
plot(tail(residuals(div.lm),n-1)~head(residuals(div.lm),n-1),xlab=expression(hat(epsilon)[i]),ylab=expression(hat(epsilon)[i+1]),main="Residual vs Residual With Lag of 1")
dev.off()

### b)
div.alm<-gls(divorce~unemployed+femlab+marriage+birth+military,correlation=corAR1(form=~year),method="ML")

print(xtable(summary(div.alm)$tTable,digits=d),type="latex",file="hw6_3_alm_coef.tex")
print(xtable(intervals(div.alm)$corStruct,digits=d),type="latex",file="hw6_3_alm_acor.tex")
print(xtable(summary(div.lm),digits=d),type="latex",file="hw6_3_lm_coef.tex")

#############################
# Problem 4 
#############################
data(cars)
attach(cars)
car.lm<-lm(dist~speed)

pdf("hw6_4_res.pdf")
plot(fitted(car.lm),residuals(car.lm),xlab="Fitted",ylab="Residuals",main="Residual Plot")
abline(h=0)
dev.off()

pdf("hw6_4_sqrtres.pdf")
sqrtres<-lm(sqrt(abs(residuals(car.lm)))~fitted(car.lm))
plot(fitted(car.lm),sqrt(abs(residuals(car.lm))),xlab="Fitted",ylab="Square Root Abs Residuals",main="Square Root Transformed Residual Plot")
abline(sqrtres)
dev.off()

var.test(residuals(car.lm)[speed<14],residuals(car.lm)[speed>=14])

pdf("hw6_4_qq.pdf")
qqnorm(residuals(car.lm),ylab="Residuals",main="Residual vs Normal Qunatiles")
qqline(residuals(car.lm))
dev.off()

pdf("hw6_4_ker.pdf")
scatter.smooth(speed,dist,main="Linear vs. Kernel Regression")
abline(car.lm,lty=2)
legend('topleft', legend=c('Linear','Kernel'), lty=c(2,1))
dev.off()

s<-20
higher<-speed>s
speed.low<-ifelse(higher,0,speed-s)
speed.high<-ifelse(higher,speed-s,0)
car.lm2<-lm(dist~speed.low + speed.high)
ana<-anova(car.lm,car.lm2)
print(xtable(ana,digits=3),type="latex",file="hw6_4.tex")

#############################
# Problem 5 
#############################
data(stackloss)

sl.lm<-lm(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss)
print(xtable(summary(sl.lm)$coefficients,digits=d),type="latex",file="hw6_5_lm.tex")

sl.hub<-rlm(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss)
print(xtable(summary(sl.hub)$coefficients,digits=d),type="latex",file="hw6_5_hub.tex")

sl.lts<-ltsreg(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss,nsamp="exact")
print(xtable(as.data.frame(coef(sl.lts)),digits=d),type="latex",file="hw6_5_lts.tex")

sl.lad<-rq(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss)
print(xtable(summary(sl.lad)$coefficients,digits=d),type="latex",file="hw6_5_lad.tex")

names<-row.names(stackloss)
n<-nrow(stackloss)
lp<-4
pdf("hw6_5_halfnorm.pdf")
halfnorm(cooks.distance(sl.lm),lp,labs=names,ylab="Cook's Distance",main="Cook's Distance vs. Half-Normal Quantiles")
dev.off()
vars<-c("Intercept",attr(terms(sl.lm),"term.labels"))
p.var.names <- c("Intercept","Acid Concentration","Water Temp","Air Flow")
ab.var.names <- c("int","ac","wt","af")
for (i in 1:4) {
    minmax<-c(min(dfbeta(sl.lm)[,i]),max(dfbeta(sl.lm)[,i]))
    cutoff<-sort(abs(dfbeta(sl.lm)[,i]))[n-lp]
    extreme<-abs(dfbeta(sl.lm)[,i])>=cutoff
    pdf(paste("hw6_5_",ab.var.names[i],".pdf",sep=""))
    plot(dfbeta(sl.lm)[!extreme,i],ylab=paste("Change in",p.var.names[i],"coef"),xlab="Index",main=paste("Influence of Observations on",p.var.names[i],"coef"),xlim=c(0,n),ylim=minmax)
    text((1:n)[extreme],dfbeta(sl.lm)[extreme,i],names[extreme])
    dev.off()
}

stud<-rstudent(sl.lm)
tval<-qt(.05/42,42-4)

sl.lm.no<-lm(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss[-21,])
print(xtable(summary(sl.lm.no)$coefficients,digits=d),type="latex",file="hw6_5_lm_nout.tex")

pdf("hw6_5_halfnorm_nout.pdf")
halfnorm(cooks.distance(sl.lm.no),lp,labs=names,ylab="Cook's Distance",main="Cook's Distance vs. Half-Normal Quantiles")
dev.off()

sl.hub.no<-rlm(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss[-21,])
print(xtable(summary(sl.hub.no)$coefficients,digits=d),type="latex",file="hw6_5_hub_nout.tex")

sl.lts.no<-ltsreg(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss[-21,],nsamp="exact")
print(xtable(as.data.frame(coef(sl.lts.no)),digits=d),type="latex",file="hw6_5_lts_nout.tex")

sl.lad.no<-rq(stack.loss~Acid.Conc.+Water.Temp+Air.Flow,data=stackloss[-21,])
print(xtable(summary(sl.lad.no)$coefficients,digits=d),type="latex",file="hw6_5_lad_nout.tex")






    



