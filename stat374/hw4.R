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
library('wavethresh')

### Problems
problem1<-F
problem2<-F
problem3<-F
problem4<-F
problem5<-T

### Options
procs<-4

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
    eps<-.001
    flies<-read.table("hw4/flies.dat",header=T)
    flies$logm<-log(flies$mort.rate+eps)

    ### a
    iso.mat<-matrix(,nrow=90,ncol=5)
    iso.mat[,1]<-flies$logm[1:90]
    days<-c(25,50,75,85)
    for (d in 1:4) {
        func<-pava(flies$logm[1:days[d]],long.out=T,stepfun=T)$h
        iso.mat[,(d+1)]<-func(0:90)[-1]
    }

    pdf('hw4/3_a.pdf')
    matplot(y=iso.mat,x=flies$day[1:90],type="l",col=rep(1,5),lty=1:5,xlab="Day",ylab="Log Mortality Rate",main="Isotonic Fit")
    legend('bottom',c("Actual",paste(days,"Day Fit")),lty=1:5)
    dev.off()
 
    ### b
    cv.score<-function(h) {
        cv.fit<-fitted(locfit(logm~lp(day,deg=1,h=h),kern="gauss",ev=dat(),data=flies),what="coef",cv=TRUE,maxk=2000)
        return(sum((y-cv.fit)^2)/n)
    }
    hseq<-exp(seq(-2,3,.25))
    cv.mat.w<-mcmapply(cv.score,h=hseq,mc.cores=procs)
    h<-hseq[cv.mat.w==min(cv.mat.w)]

    days<-c(25,50,75,85,171)
    for (d in days) {
        predmat.ll<-matrix(,nrow=d+1,ncol=4)
        m.ll<-predict(locfit(logm~lp(day,deg=1,h=h),kern="gauss",ev=dat(),data=flies[1:(d+1),]),se.fit=T)
        predmat.ll<-cbind(flies$logm[1:(d+1)],m.ll$fit,m.ll$fit-crit(rob.ll.w)$crit.val*m.ll$se.fit,m.ll$fit+crit(rob.ll.w)$crit.val*m.ll$se.fit)
        pdf(paste('hw4/3_b',d,'day.pdf',sep=""))
        matplot(x=flies$day[1:(d+1)],y=predmat.ll,col=rep(1,4),type="l",lty=c(1,2,3,3),xlab="Day",,ylab="Log Mortality Rate",main=paste(d,"Day Local Linear Fit"))
        legend('bottom',c("Actual","Predicted","95% CI"),lty=1:3)
        dev.off()
    }

    ### c
    soft.thresh<-function(x,t) sign(x)*ifelse(abs(x)-t>0,abs(x)-t,0)
    J<-7
    n<-2^J

    fs<-flies[1:n,]
    alpha<-mean(fs$logm)

    Z<-list()
    haar<-list()
    siv<-list()
    for (j in 0:(J-1)) {
        # The full sequence of haar waveletes at a particular level
        haar[[j+1]]<-2^(j/2)*rep(c(rep(-1,2^(J-j-1)),rep(1,2^(J-j-1))),2^j)
        # This identifies disitinct haar waveletes within the full sequence
        siv[[j+1]]<-ceiling(1:n/(2^(J-j)))
        # Generates the initial estimate of the coefficients 
        Z[[j+1]]<-(1/n)* aggregate(fs$logm*haar[[j+1]],by=list(siv[[j+1]]),FUN=sum)$x
    }
    # MAD
    sigma<-sqrt(n)*median(abs(Z[[J]]-median(Z[[J]])))/.6745
    # use universal threshold
    lambda<-sigma*sqrt(2*log(n)/n)
    beta<-list()
    Yhat<-rep(alpha,n)
    for (j in 0:(J-1)) {
        #apply threshold to get estimates
        beta[[j+1]]<-soft.thresh(Z[[j+1]],lambda)
        # Add the jth level of haar wavelets on
        Yhat<-Yhat+haar[[j+1]]*rep(beta[[j+1]],rep(2^(J-j),2^j))
    }

    wav.mat<-cbind(fs$logm,Yhat)
    pdf('hw4/3_c.pdf')
    matplot(y=wav.mat,x=fs$day,type="l",col=c(1,1),lty=c(1,2),xlab="Day",ylab="Log Mortality Rate",main="Haar Wavelet Fit")
    legend('bottom',c("True Rate","Haar Wavelet Fit"),lty=c(1,2))
    dev.off()
    
    ### d
    da.wav<-wd(fs$logm,family="DaubLeAsymm",filter.number=8)
    da.wav.thr<-threshold(da.wav,type="soft",policy="universal")
    da.wav.mat<-cbind(fs$logm,wr(da.wav.thr))
    pdf('hw4/3_d.pdf')
    matplot(y=da.wav.mat,x=fs$day,col=c(1,1),type="l",lty=c(1,2),xlab="Day",ylab="Log Mortality Rate",main="Daubechies Wavelet Fit")
    legend('topright',c("True Rate","Daubechies Wavelet Fit"),lty=c(1,2))
    dev.off()
    
}
############# Problem 4 
if (problem4) {
    n<-10
    d<-1000
    b<-1000
    est<-matrix(,nrow=b,ncol=2)
    for (i in 1:b) {
        x<-replicate(n,rnorm(d))
        xbar<-apply(x,1,mean)
        est[i,]<-sum(xbar^2)+c(d/n,-d/n)  
    }
    pdf('hw4/4_bayes.pdf')
    lim<-c(min(est),max(est))
    hist(est[,1],main="Bayes Estimate",xlab="Prediction",xlim=lim)
    dev.off()
    pdf('hw4/4_freq.pdf')
    hist(est[,2],main="Frequentist Estimate",xlab="Prediction",xlim=lim)
    dev.off()
}

############# Problem 5 
if (problem5) {
    dirc<-function(rfunc,alpha,n,...) {
        step.mat<-matrix(,nrow=n,ncol=2)
        step.mat[,1]<-rfunc(n,...)
        V<-rbeta(n,1,alpha) 
        step.mat[,2]<-cumprod(c(1,1-V))[-(n+1)]*V
        order.mat <- step.mat[order(step.mat[,1]),] 
        return(stepfun(order.mat[,1],cumsum(c(0,order.mat[,2]))))
    }
    
    ### d
    alpha<-c(1,1,1,5,5,5,20,20,50,50)
    n<-200
    pp<-seq(-4,4,length.out=300)
    mat<-matrix(,nrow=300,10)
    for (i in 1:10) {
        func<-dirc(rnorm,alpha[i],n)
        mat[,i]<-func(pp)
    }
    pdf('hw4/5_d.pdf')
    matplot(x=pp,y=mat,type="l",col=rep(1,10),lty=c(1,1,1,2,2,2,3,3,4,4),main="Dirchlet Processes",xlab="x",ylab="Cumulative Probability")
    legend('bottomright',paste("Alpha = ",c(1,5,20,50)),lty=1:4)
    dev.off()

    ### e
    pp<-seq(-7,17,length.out=300)
    type2<-list()
    for (n in c(10,25,100)) {
        #i
        mat<-matrix(,nrow=300,ncol=3)
        X<-sort(rnorm(n,5,3))
        mat[,1]<-stepfun(X,(0:n)/n)(pp)
        e<-sqrt(log(2/.05)/(2*n))
        mat[,2]<-pmax(mat[,1]-e,0)
        mat[,3]<-pmin(mat[,1]+e,1)
        pdf(paste('hw4/5_ei',n,'.pdf',sep=""))
        matplot(mat,pp,type="l",col=rep(1,3),lty=c(1,3,3),main=paste("Emperical Distribution & DKW, N=",n),xlab="x",ylab="Cumulative Probability")
        dev.off()
        
        #ii 
        alpha<-100
        posterior<-function(x,X,func,alpha,N,...) {
            return((N/(N+alpha))*stepfun(X,(0:N)/N)(x)+(alpha/(N+alpha))*func(x,...))
        }
        rposterior<-function(num,X,func,alpha,N,...) {
            return(ifelse(runif(n)>(N/(N+alpha)),sample(X,n,replace=TRUE),func(n,...)))
        }
        post<-posterior(pp,X,pnorm,alpha,n)
        test<-replicate(100,dirc(rposterior,alpha,n=100,num=100,X=X,func=rnorm,alpha=alpha,mean=0,sd=1)(pp))
        ub<-apply(test,1,quantile,.95)
        lb<-apply(test,1,quantile,.05)
        mat<-cbind(post,ub,lb,replicate(4,dirc(rposterior,alpha,100,num=100,X=X,func=rnorm,alpha=alpha,mean=0,sd=1)(pp)))
        for (i in 2:ncol(mat)) mat[,i]<-mat[,i]/max(mat[,i])
        pdf(paste('hw4/5_eii',n,'.pdf',sep=""))
        matplot(y=mat,x=pp,type="l",col=rep(1,3),lty=c(1,3,3,rep(2,4)),main=paste("Posterior, N=",n),xlab="x",ylab="Cumulative Probability")
        legend('bottomright',c("Posterior","Additional Draws","CI"),lty=1:3)
        dev.off()
    }
        
        
        
        
}

