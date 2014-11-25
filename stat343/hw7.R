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





