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


### Problem 1
# Data
dose        <- c(   2, 2.64,3.48,4.59,6.06,   8)
ddt_num     <- c(   3,    5,  19,  19,  24,  35)
ddt_denom   <- c(  50,   49,  47,  38,  49,  50)
bhc_num     <- c(   2,   14,  20,  27,  41,  40)
bhc_denom   <- c(  50,   49,  50,  50,  50,  50)
both_num    <- c(  28,   37,  46,  48,  48,  50)
both_denom  <- c(  50,   50,  50,  50,  50,  50)

insect      <- data.frame(
                            dose        = rep(dose,3),
                            insecticide = rep(c('DDT','BHC','DDT+BHC'),each=6),
                            died        = c(ddt_num,bhc_num,both_num),
                            pop         = c(ddt_denom,bhc_denom,both_denom)
                         )
insect$lived <- insect$pop - insect$died

## (a) & (b)
ins.logit.b <- glm(cbind(died,lived)~insecticide+dose,family=binomial,data=insect)
n           <- 500
xseq        <- seq(0,9,length.out=n)

pdf('hw3/A1_ab.pdf')
matplot(dose,cbind(ddt_num/ddt_denom,bhc_num/bhc_denom,both_num/both_denom),lty=c(1,1,1),col=2:4,pch=c(1,1,1),main="Insect Death Rate With Logistic Fit",ylab="Death Rate",xlab=expression(paste('Dose Insecticide (mg/10c',m^2,')')))
points(xseq,predict(ins.logit.b,data.frame(dose=xseq,insecticide=rep('DDT',n)    ),type='response'),type='l',col=2)
points(xseq,predict(ins.logit.b,data.frame(dose=xseq,insecticide=rep('BHC',n)    ),type='response'),type='l',col=3)
points(xseq,predict(ins.logit.b,data.frame(dose=xseq,insecticide=rep('DDT+BHC',n)),type='response'),type='l',col=4)
legend('bottomright',c('DDT',expression(paste(gamma,'-BHC')),expression(paste('DDT + ',gamma,'-BHC'))),col=2:4,pch=c(1,1,1))
dev.off()

## (c)
ins.logit.c     <- glm(cbind(died,lived)~insecticide+log(dose),family=binomial,data=insect)
ins.logit.int.c <- glm(cbind(died,lived)~insecticide*log(dose),family=binomial,data=insect)
lrtest(ins.logit.c,ins.logit.int.c)

## (e) - Fieller's method
ins.logit.e     <- glm(cbind(died,lived)~insecticide+log(dose)-1,family=binomial,data=insect)
var1 = 'insecticideDDT+BHC'
vecs <- list(c(-1,0,1,0),c(0,-1,1,0))
CIS  <- list(NULL,NULL)
for (i in 1:2){
    a <- vecs[[i]] %*% coef(ins.logit.e)
    b <- coef(ins.logit.e)['log(dose)']
    sa <- t(vecs[[i]]) %*% vcov(ins.logit.e) %*% vecs[[i]]
    sab <- t(vecs[[i]]) %*%vcov(ins.logit.e) %*% c(0,0,0,1)
    sb <- vcov(ins.logit.e)['log(dose)','log(dose)']
    est <-a/b
    lb  <-((a*b - 1.644^2*sab) - sqrt((a*b - 1.644^2*sab)^2-(b^2-1.644^2*sb)*(a^2-1.644^2*sa)))/(b^2-1.644^2*sb)
    ub  <-((a*b - 1.644^2*sab) + sqrt((a*b - 1.644^2*sab)^2-(b^2-1.644^2*sb)*(a^2-1.644^2*sa)))/(b^2-1.644^2*sb)
    CIS[[i]]<-exp(c(lb,est,ub))
}
    
## (f) - Try other links
ins.ll.f  <- glm(cbind(lived,died)~insecticide+log(dose)-1,family=binomial(link='cloglog'),data=insect)
ins.cll.f <- glm(cbind(died,lived)~insecticide+log(dose)-1,family=binomial(link='cloglog'),data=insect)
ins.pro.f <- glm(cbind(died,lived)~insecticide+log(dose)-1,family=binomial(link='probit' ),data=insect)

var1 = 'insecticideDDT+BHC'
vecs <- list(c(-1,0,1,0),c(0,-1,1,0))
CIS2 <- list(NULL,NULL)
for (i in 1:2){
    a <- vecs[[i]] %*% coef(ins.cll.f)
    b <- coef(ins.cll.f)['log(dose)']
    sa <- t(vecs[[i]]) %*% vcov(ins.cll.f) %*% vecs[[i]]
    sab <- t(vecs[[i]]) %*%vcov(ins.cll.f) %*% c(0,0,0,1)
    sb <- vcov(ins.cll.f)['log(dose)','log(dose)']
    est <-a/b
    lb  <-((a*b - 1.644^2*sab) - sqrt((a*b - 1.644^2*sab)^2-(b^2-1.644^2*sb)*(a^2-1.644^2*sa)))/(b^2-1.644^2*sb)
    ub  <-((a*b - 1.644^2*sab) + sqrt((a*b - 1.644^2*sab)^2-(b^2-1.644^2*sb)*(a^2-1.644^2*sa)))/(b^2-1.644^2*sb)
    CIS2[[i]]<-exp(c(lb,est,ub))
}

## (g) estimate dose for 99% kill rate
ins.logit.g <- glm(cbind(died,lived)~insecticide+log(dose)-1,offset=rep(-logit(.99),18),family=binomial(link='probit'),data=insect)
a <- coef(ins.logit.g)[3]
b <- coef(ins.logit.g)['log(dose)']
sa <- vcov(ins.logit.g)[3,3]
sab <- vcov(ins.logit.g)[3,4]
sb <- vcov(ins.logit.g)['log(dose)','log(dose)']
est <-a/b
lb  <-((a*b - 1.644^2*sab) - sqrt((a*b - 1.644^2*sab)^2-(b^2-1.644^2*sb)*(a^2-1.644^2*sa)))/(b^2-1.644^2*sb)
ub  <-((a*b - 1.644^2*sab) + sqrt((a*b - 1.644^2*sab)^2-(b^2-1.644^2*sb)*(a^2-1.644^2*sa)))/(b^2-1.644^2*sb)
CIS3<-exp(c(lb,est,ub))

### Problem 2
children <- read.csv('hw3/hw3_2.csv')
child.poi <- glm(round(Count*Mean)~Age+as.factor(Education)+Rural,offset=log(Count),data=children,family='poisson')
child.poi2 <- glm(round(Count*Mean)~Age+Rural,offset=log(Count),data=children,family='poisson')
child.poi3 <- glm(round(Count*Mean)~Age*Rural+as.factor(Education),offset=log(Count),data=children,family='poisson')

# 95% for urban woman, upper elementary, 10 years
testa <- data.frame(Age='10–14',Rural=0,Education=3,Count=1)
predictiona <- predict(child.poi,testa,se.fit=T)
CI2a <- exp(c(-1.96,0,1.96)*predictiona$se.fit + predictiona$fit)

# rural woman, secondary education ,lifetime
testb <- data.frame(Age='25+',Rural=1,Education=4,Count=1)
predictionb <- predict(child.poi,testb,se.fit=T)
CI2b <- exp(c(-1.64,0,1.64)*predictionb$se.fit + predictionb$fit)


