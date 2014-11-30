###### Header
library('parallel')
library('locfit')
library('gam')
library('MASS')
library('xtable')
library('reshape2')
library('chron')
library('graphics')
library('Iso')
library('wavelets')

### Problems
problem1<-F
problem2<-F
problem3<-T

### Options
procs<-2

############# Problem 1
if (problem1) {
    ### a
    n<-500
    th<-2
    b<-100
    risk<-function(n,th) {
        z<-rnorm(n)+th
        hat<-pava(z)
        return(c(n*(th-mean(z))^2,sum((th-hat)^2)))
    }
    est<-c(0,0)
    for (i in 1:b){
        est<-est+risk(n,th)
    }
    est<-est/b
    
    ### b
    n<-seq(20,200,4)
    sig<-.25
    b<-500 
    risk.step<-function(n,sig) {
        step<-c(rep(0,n/2),rep(1,n/2))
        z<-rnorm(n)*sig+step
        hat<-pava(z)
        return(sum((hat-step)^2))
    }
    risk.line<-function(n,sig) {
        line<-1/(n:1)
        z<-rnorm(n)*sig+line
        hat<-pava(z)
        return(sum((hat-line)^2))
    }
    risk<-matrix(0,nrow=length(n),ncol=2)
    for (i in 1:b) {
        risk[,1]<- risk[,1] + mcmapply(risk.step,n=n,MoreArgs=list(sig=sig),mc.cores=procs)
        risk[,2]<- risk[,2] + mcmapply(risk.line,n=n,MoreArgs=list(sig=sig),mc.cores=procs)
    }
    risk<-risk/b
    pdf('hw4/1_b.pdf')
    matplot(x=n,y=risk,col=c(1,1),xlab="N",type="l",ylab=paste("Average Risk Over",b,"Simulations"),main="Monotone Regression Risk")
    legend('bottomright',c("Step Function","Line"),lty=c(1,2))
    dev.off()
    
}
############# Problem 2 
if (problem2) {
    doppler<-function(x) sqrt(x*(1-x))*sin(2.1*pi/(x+.05))
    soft.thresh<-function(x,t) sign(x)*ifelse(abs(x)-t>0,abs(x)-t,0)
    n<-1024
    x<-(1:n)/n
    ### a)
    error<-list(rnorm(n)*.01,rnorm(n)*.1,rnorm(n)*ifelse(runif(n)>.05,.1,4))
    sig<-c(.01,.1,sqrt(.95*.01^2+.05*4^2))
    names<-c("pt01","pt1","bimodal")
    titles<-c("Variance of .01","Variance of .1","Variance of .01 or 4")
    models<-c("Actual","James Stein","Universal Threshold","SureShrink","Local Linear")
    msqe<-matrix(,nrow=3,ncol=5)
    rownames(msqe)<-titles
    colnames(msqe)<-models
    fitmats<-list()


    for (i in 1:3) {
        y<-doppler(x)+error[[i]]
        fit<-matrix(,nrow=n,ncol=5)
        fit[,1]<-doppler(x)

        # local linear
        cv.score<-function(h) {
            cv.fit<-fitted(locfit(y~lp(x,deg=1,h=h),kern="gauss",ev=dat()),what="coef",cv=TRUE,maxk=2000)
            return(sum((y-cv.fit)^2)/n)
        }
        hseq<-seq(.002,.02,.002)
        cv.mat.w<-mcmapply(cv.score,h=hseq,mc.cores=procs)
        h<-hseq[cv.mat.w==min(cv.mat.w)]
        ll<-locfit(y~lp(x,deg=1,h=h),kern="gauss",ev=dat())
        fit[,5]<-predict(ll)

        # wavelets
        wav.tr<-dwt(y)

        wav.us<-wav.tr
        lambda.us<-sig[i]*sqrt(2*log(n)/n)
        wav.js<-wav.tr
        wav.ss<-wav.tr
        sure.ss<-function(lam,n,sig,beta) sum(sig^2/n-2*sig^2/n*(abs(beta)<=lam)+min(beta^2,lam^2))
        lambda.ss<-NULL
        for (j in 1:wav.tr@level) {
            nj<-length(wav.tr@W[[j]])
            # universal shrinkage
            wav.us@W[[j]]<-soft.thresh(wav.tr@W[[j]],lambda.us)
            # james stein
            wav.js@W[[j]]<-(1-((n-2)*sig[i]/n)/sum(wav.tr@W[[j]]^2))*wav.tr@W[[j]]
            # sure shrink
            lambda.ss[[j]]<-optim(par=lambda.us,fn=sure.ss,lower=0,upper=(sig[i]/sqrt(nj))*sqrt(2*log(nj))^3,method="L-BFGS-B",n=n,sig=sig[i],beta=wav.tr@W[[j]])$par
            wav.ss@W[[j]]<-soft.thresh(wav.tr@W[[j]],lambda.ss[[j]])
        }
        fit[,2]<-idwt(wav.js)
        fit[,3]<-idwt(wav.us)
        fit[,4]<-idwt(wav.ss)

        # plot 
        pdf(paste('hw4/2_ab_',names[i],".pdf",sep=""))
        matplot(x=x,y=fit[,c(5,4)],col=rep(1,2),lty=c(1,3),type="l",main=titles[i],xlab="x",ylab="y")
        legend('bottomright',models[c(5,4)],lty=c(1,3))
        dev.off()

        msqe[i,]<-colSums((fit-fit[,1])^2)
        fitmats[[i]]<-fit
    
    }
    print(xtable(msqe[,-1],digits=3),type="latex",file="hw4_2.tex")
    
}
    
############# Problem 3
if (problem3) {
}


