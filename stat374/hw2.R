library('parallel')
library('locfit')
library('gam')

############################################
# Problem 1
############################################
### part c)
# Input data and transform
Y<-c(.4,.378,.356,.333,.311,.311,.289,.267,.244,.244,.222,.222,.222,.222,.222,.200,.178,.156)
p<-c(.346,.298,.276,.222,.273,.270,.263,.201,.269,.230,.264,.256,.303,.264,.226,.285,.316,.200)
Z<-sqrt(45)*asin(2*Y-1)
th<-sqrt(45)*asin(2*p-1)

# Compute estimators and calculate risk
n<-length(Y)
mean <- mean(Z)
mle <- Z 
JS0 <- Z*(1-(n-2)*1/(sum(Z^2)))
JSm <- mean + (Z-mean)*(1-(n-2)*1/(sum((Z-mean)^2)))
risk <- c(sum((mle-th)^2),sum((JS0-th)^2),sum((JSm-th)^2))

############################################
# Problem 2 
############################################

test<-read.table('test.txt',header=TRUE)
train<-read.table('train.txt',header=TRUE)
attach(train)

cas.log<-log(casual+1)
reg.log<-log(registered+1)
ct.log<-log(count+1)
season.f<-factor(season)
holiday.f<-factor(holiday)
workingday.f<-factor(workingday)
month.f<-factor(month)
weather.f<-factor(weather)
training<-day<=15
validate<-day<=19 & day>15
hour.mod<-abs(hour-17)*abs(hour-8)
hour.mod.2<-abs(hour-17)

my.data<-data.frame(train,cas.log,reg.log,ct.log,season.f,holiday.f,workingday.f,month.f,weather.f,training,validate,hour.mod,hour.mod.2)

reg.lm <- lm(reg.log~season.f*workingday.f+weather.f+hour.mod.2+daylabel, data=my.data[training,])
cas.lm <- lm(cas.log~season.f*workingday.f+weather.f+hour.mod.2+daylabel, data=my.data[training,])

reg.er<-sum((my.data$reg.log[my.data$validate] -predict(reg.lm,my.data[validate,]))^2)
cas.er<-sum((my.data$cas.log[my.data$validate] -predict(cas.lm,my.data[validate,]))^2)

### part b)

reg.dm<-tapply(my.data$reg.log[training],my.data$daylabel[training],mean)
cas.dm<-tapply(my.data$cas.log[training],my.data$daylabel[training],mean)


reg.loc<-locfit(reg.log ~ daylabel, data=my.data[training,])
pdf('hw2_2b_reg_plot.pdf')
plot(as.integer(names(reg.dm)),reg.dm, main="Mean Log Registered by Day",xlab="Day Label",ylab="Log Registered")
plot(reg.loc,add=TRUE)
dev.off()

reg.res <- residuals(reg.loc, type="raw", data=my.data[training,])
reg.lm.res <- lm(reg.res~season.f*workingday.f+weather.f+hour.mod.2+daylabel, data=my.data[training,])
reg.loc.pred <-predict(reg.loc, my.data[validate,]) + predict(reg.lm.res,my.data[validate,])
reg.loc.er<-sum((my.data$reg.log[my.data$validate] -reg.loc.pred)^2)


cas.loc<-locfit(cas.log ~ daylabel, data=my.data[training,])
pdf('hw2_2b_cas_plot.pdf')
plot(as.integer(names(cas.dm)),cas.dm, main="Mean Log Casual by Day",xlab="Day Label",ylab="Log Registered")
plot(cas.loc,add=TRUE)
dev.off()

cas.res <- residuals(cas.loc, type="raw", data=my.data[training,])
cas.lm.res <- lm(cas.res~season.f*workingday.f+weather.f+hour.mod.2+daylabel, data=my.data[training,])
cas.loc.pred <-predict(cas.loc, my.data[validate,]) + predict(cas.lm.res,my.data[validate,])
cas.loc.er<-sum((my.data$cas.log[my.data$validate] -cas.loc.pred)^2)

### part c)
reg.gam<-gam(reg.log~s(daylabel)+s(humidity)+s(temp)+s(windspeed)+s(hour),data=my.data[training,])
reg.gam.pred <-predict(reg.gam, my.data[my.data$validate,])
reg.gam.er<-sum((my.data$reg.log[my.data$validate] -reg.gam.pred)^2)

### part d)
reg.fin<-gam(reg.log~s(daylabel)+s(humidity)+s(temp)+s(windspeed)+s(hour)+ workingday.f + weather.f*season.f,data=my.data[training,])
reg.fin.pred <-predict(reg.fin, my.data[my.data$validate,])
reg.fin.er<-sum((my.data$reg.log[my.data$validate] -reg.fin.pred)^2)

cas.fin<-gam(cas.log~s(daylabel)+s(humidity)+s(temp)+s(windspeed)+s(hour)+ workingday.f + weather.f*season.f,data=my.data[training,])
cas.fin.pred <-predict(cas.fin, my.data[my.data$validate,])
cas.fin.er<-sum((my.data$cas.log[my.data$validate] -cas.fin.pred)^2)

ct.fin<-gam(ct.log~s(daylabel)+s(humidity)+s(temp)+s(windspeed)+s(hour)+ workingday.f + weather.f*season.f,data=my.data[training,])
ct.fin.pred <-predict(ct.fin, my.data[my.data$validate,])
ct.fin.er<-sum((my.data$ct.log[my.data$validate] -ct.fin.pred)^2)

sum.fin.er<-sum((my.data$ct.log[my.data$validate]-reg.fin.pred -cas.fin.pred)^2)

### part d
final.model<-gam(ct.log~s(daylabel)+s(humidity)+s(temp)+s(windspeed)+s(hour)+ workingday.f + weather.f*season.f,data=my.data)

attach(test)

cas.log<-log(casual+1)
reg.log<-log(registered+1)
ct.log<-log(count+1)
season.f<-factor(season)
holiday.f<-factor(holiday)
workingday.f<-factor(workingday)
month.f<-factor(month)
weather.f<-factor(weather)
training<-day<=15
validate<-day<=19 & day>15
hour.mod<-abs(hour-17)*abs(hour-8)
hour.mod.2<-abs(hour-17)

final<-exp(predict(final.model,newdata=test))-1
write.table(final,file="assn2-probx-ajmaurer.txt", row.names=FALSE, col.names=FALSE)


############################################
# Problem 4
############################################
glass<-read.table('glass.dat',header=TRUE)

### a
hist.den<-function(X,b) {
    min<-min(X)-1e-5
    w<-(max(X)+1e-5-min)/b
    n<-nrow(X)
    den<-sapply(1:b, function(bn) sum((w*bn+min>=X) & (X>w*bn-w+min))/(n))
    return(den)
}

cv.hist<- function(X,B) {
    n <- nrow(X)
    nb <- length(B)
    span<-max(X)-min(X)-1e-5
    score<-function(buc) {
        p<-hist.den(X,buc)
        h<-1/buc
        return( 2/(n*h-h)-(n+1)*sum(p^2)/(n*h-h))
    }
    cv.score<-sapply(B,score)
    return(cv.score)
}

if (FALSE) {
    for (i in 1:(ncol(glass)-1)) {
        cv.score<-cv.hist(glass[i],1:200)
        b<-which.min(cv.score)
        pdf(paste("hw2_4a_",names(glass)[i],"_hist.pdf",sep=""))
        hist(as.matrix(glass[i]),breaks=b-1, freq=FALSE, plot=TRUE, main=paste("Histogram of ",names(glass)[i],", ",toString(b)," bins", sep=""), xlab=names(glass)[i])
        dev.off()
    }
    
### b    

for (i in 1:(ncol(glass)-1)) {
    gau<-density(as.matrix(glass[i]),bw="ucv",kernel="gaussian")
    X<-as.matrix(glass[i])
    J=function(h){
        fhat=Vectorize(function(x) density(X,from=x,to=x,n=1,bw=h,kernel=ker)$y)
        fhati=Vectorize(function(i) density(X[-i],from=X[i],to=X[i],n=1,bw=h)$y)
        F=fhati(1:length(X))
        return(integrate(function(x) fhat(x)^2,-Inf,Inf,subdivisions=300,rel.tol = 1e-2)$value-2*mean(F))
    }
    steps<-40
    delta<-min(.015,gau$bw/3)
    h<-max(gau$bw-steps*delta/2,delta)+ delta*(0:steps)

    ker<-"rectangular"
    boxbw<-h[which.min(sapply(h,J))]
    ker<-"epanechnikov"
    epabw<-h[which.min(sapply(h,J))]

    box<-density(as.matrix(glass[i]),bw=boxbw,kernel="rectangular")
    epa<-density(as.matrix(glass[i]),bw=epabw,kernel="epanechnikov")
    pdf(paste("hw2_4b_",names(glass)[i],"_den.pdf",sep=""))
    plot(box,main=paste(names(glass)[i],"Density Estimates"),xlab=names(glass)[i], ylab="Density",lty=1)
    lines(gau, lty=2)
    lines(epa, lty=3)
    legend('topright', legend=c(paste("Boxcar, BW=",format(box$bw,digits=3)),paste("Gaussian, BW=",format(gau$bw,digits=3)),paste("Epanechnikov, BW=",format(epa$bw,digits=3))), lty=c(1,2,3))
    dev.off()
}
}


