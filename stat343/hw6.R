library('parallel')
library('alr3')
library('faraway')
library('car')
library('MASS')
library('locfit')
library('xtable')
library('nlme')

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

d<- 3
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
line(x)
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






    



