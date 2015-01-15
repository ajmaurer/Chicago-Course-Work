library('xtable')
library('stats4')
library('locfit')


### Problem 1
tb1 <-data.frame(row.names=c(1,2,3,4,5,6))
tb1$sz  <-c(  35,  31,  33,  34,  28,  34)
tb1$ym  <-c( 610, 685, 694, 750, 732, 710)
tb1$yse <-c(  25,  39,  67,  41,  39,  54)
tb1$tm  <-c( 478, 442, 414, 375, 395, 411)
tb1$tse <-c(  26,  37,  64,  44,  28,  53)
tb1$n   <-c(  16,  14,  17,  14,   9,  14) 
tb1$p   <-c(.039,.011,.031,.002,.008,.026)

# b
tb1$dm<-tb1$ym-tb1$tm
tb1$dse<-sqrt(tb1$yse^2+tb1$tse^2)

b.out<-t(as.matrix(tb1$dse))
colnames(b.out)<-1:6
rownames(b.out)<-"SE"
print(xtable(b.out,digits=1),type="latex",file="hw1/p_1_a.tex")

# c
c.SE<-sqrt(sum((tb1$yse^2+tb1$tse^2)*tb1$n^2))/sum(tb1$n)
c.m<-mean(tb1$dm)
c.CI<-c.m+1.96*c.SE*c(-1,1)

### 2
s<-list()
s$s1<-c(76 ,82,83,54 ,35 ,46 ,87 ,68)
s$s2<-c(87 ,95,98,100,109,109,100,81 ,75 ,68 ,67)
s$s3<-c(105,83,76,75 ,51 ,76 ,93 ,75 ,62 )   
s$s4<-c(95 ,90,76,76 ,87 ,79 ,77 ,71 )       
s$s5<-c(76 ,76,78,79 ,72 ,68 ,75 ,78 )       
s$s6<-c(78 ,78,78,86 ,87 ,81 ,73 ,67 ,75 ,82 ,83 )
s$s7<-c(82 ,79,81,79 ,77 ,79 ,79 ,78 ,79 ,82 ,76 ,73 ,64)
s$s8<-c(84 ,86,85,82 ,77 ,76 ,77 ,80 ,83 ,81 ,78 ,78 ,78)

# b
s.ul <- unlist(s)
n<-length(s.ul)
k<-8
df1 <- k-1
df2 <- n-k 
ss1 <- 0
ss2 <- 0
s$mean <- mean(s.ul)
s$means<- rep(NA,k)
for (i in 1:k) {
    s$means[i] <- mean(s[[i]])
    ss1<-ss1+ length(s[[i]])*(s$mean-s$means[i])^2
    ss2<-ss2+sum((s[[i]]-s$means[i])^2)
}
f<-(ss1/df1)/(ss2/df2)
p<-1-pf(f,df1,df2)

# c
s$sd<-sd(s.ul)
s$sds<- rep(NA,k)

LL0<-function(mu,sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8) {
    LL <- 0
    lsigma<-c(sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8)
    for (i in 1:8) {
        LL <- LL -sum(log(dnorm(s[[i]],mu,exp(lsigma[i]))))
    }
    return(LL)
}
LL1<-function(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8,sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8) {
    LL <- 0
    lmu<-c(mu1,mu2,mu3,mu4,mu5,mu6,mu7,mu8)
    lsigma<-c(sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8)
    for (i in 1:8) {
        LL <- LL -sum(log(dnorm(s[[i]],lmu[i],exp(lsigma[i]))))
    }
    return(LL)
}

w<-log(s$sd)

mle0<-mle(LL0,start=list(mu=s$mean,sig1=w,sig2=w,sig3=w,sig4=w,sig5=w,sig6=w,sig7=w,sig8=w))
mle0.coef<-c(coef(mle0)[1],exp(coef(mle0)[-1]))
print(xtable(t(as.matrix(mle0.coef)),digits=1),type="latex",file="hw1/p_2_c0.tex")
mle0.ll <- -summary(mle0)@m2logL

mle1<-mle(LL1,start=list(mu1=s$mean,mu2=s$mean,mu3=s$mean,mu4=s$mean,mu5=s$mean,mu6=s$mean,mu7=s$mean,mu8=s$mean,sig1=w,sig2=w,sig3=w,sig4=w,sig5=w,sig6=w,sig7=w,sig8=w))
mle1.coef<-rbind(coef(mle1)[1:k],exp(coef(mle1)[(k+1):(2*k)]))
colnames(mle1.coef)<-1:8
rownames(mle1.coef)<-c('mu','sigma')
print(xtable(mle1.coef,digits=1),type="latex",file="hw1/p_2_c1.tex")
mle1.ll <- -summary(mle1)@m2logL

LLd<- mle1.ll-mle0.ll
p.2c<-1-pchisq(LLd,9-2) 

#d
sn<-data.frame(row.names=1:13)
sn$s1<-c(76 ,82,83,54 ,35 ,46 ,87 ,68 ,NA ,NA ,NA ,NA ,NA)
sn$s2<-c(87 ,95,98,100,109,109,100,81 ,75 ,68 ,67 ,NA ,NA)
sn$s3<-c(105,83,76,75 ,51 ,76 ,93 ,75 ,62 ,NA ,NA ,NA ,NA)
sn$s4<-c(95 ,90,76,76 ,87 ,79 ,77 ,71 ,NA ,NA ,NA ,NA ,NA)
sn$s5<-c(76 ,76,78,79 ,72 ,68 ,75 ,78 ,NA ,NA ,NA ,NA ,NA)
sn$s6<-c(78 ,78,78,86 ,87 ,81 ,73 ,67 ,75 ,82 ,83 ,NA ,NA)
sn$s7<-c(82 ,79,81,79 ,77 ,79 ,79 ,78 ,79 ,82 ,76 ,73 ,64)
sn$s8<-c(84 ,86,85,82 ,77 ,76 ,77 ,80 ,83 ,81 ,78 ,78 ,78)

pdf('hw1/p_2_d.pdf')
matplot(matrix(cumsum(!is.na(sn)),nrow=13),sn,type="l",lty=rep(1,8),col=c(1:6,8:9),xlab="Index",ylab=expression(paste("Deviations of ",10^-5," From 9.8 m/s"^2)),main="Measurements of g")
legend('bottomright',paste('Series',1:8),lty=rep(1,8),col=c(1:6,8:9))
dev.off()

#e
sink('hw1/p2e.txt')
summary(gls(s.ul~1,weight=varConstPower(1,form=~A),correlation=corAR1(form=~A)))
sink()

