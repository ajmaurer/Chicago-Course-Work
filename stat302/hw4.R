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
library("qvalue")

### Caller
problem1 <- FALSE 
problem2 <- FALSE
problem3 <- TRUE

### Problem 1
## ii - simulate simple hierarchical model
p.values <- function(x){
    return(t.test(x)$p.value)
}

simulateD <- function(pi0,sigb,n=10,m=1000,sigj=1){
    output <- list()
    output$beta     <- ifelse(rbinom(m,1,pi0),0,rnorm(m,0,sigb))
    output$D        <- replicate(n,rnorm(m,output$beta,sigj))
    output$p.values <- apply(output$D,1,p.values)
    return(output)
}
   
if (problem1){
    ## iii - generate simple p-values
    scenarios <- list(c(1,1),c(.5,3),c(0,3))    # each is a pi0, sigb pair
    simulations <- list()
    for (i in 1:3){
        pi0  <- scenarios[[i]][1]
        sigb <- scenarios[[i]][2]
        simulations[[i]] <- simulateD(pi0,sigb)

        pdf(paste('hw4/pr1_iii_hist_',letters[i],'.pdf',sep=''))
        hist(simulations[[i]]$p.values,breaks=20,freq=F,main=bquote("Histogram, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab='p-value',ylab='Density')
        dev.off()
    }
}

## iv - implement Benjamini-Hochberg
BH <- function(pvals,alpha){
    m    <- length(pvals)
    ord  <- order(pvals)
    test <- pvals[ord]< (1:m)*alpha/m
    test[ord] <- (1:m)<Position(identity,test,right=T,nomatch=0)
    return(test)
}

## v - implement FDR
FDR <- function(beta,gamma,pFDR=FALSE){
    R <- sum(gamma)
    V <- sum(gamma & beta==0)
    FDR <- ifelse(R>0,V/R,ifelse(pFDR,NA,0))
    return(FDR)
}

if (problem1){
    ## vi - simulate actual FDR
    k <- 1000
    scenarios <- list(c(1,1),c(.5,3),c(0,3))    # each is a pi0, sigb pair
    alphas    <- seq(.05,.5,.05)
    FDRs      <- list()
    pFDRs     <- list()
    simulations.k <- list()
    for (i in 1:3){
        simulations.k[[i]] <-list()
        pi0  <- scenarios[[i]][1]
        sigb <- scenarios[[i]][2]
        sims <- mclapply(1:k,function(x){simulateD(pi0,sigb)},mc.cores=4)
        simulations.k[[i]]$p.values <- mclapply(sims,function(x){x$p.values},mc.cores=4)
        simulations.k[[i]]$betas    <- mclapply(sims,function(x){x$beta}    ,mc.cores=4)
        simulations.k[[i]]$gammas<-list()
        simulations.k[[i]]$FDRs  <-NULL
        simulations.k[[i]]$pFDRs <-NULL
        for (j in 1:length(alphas)) {
            simulations.k[[i]]$gammas[[j]] <- mclapply(simulations.k[[i]]$p.values,BH,alpha=alphas[j],mc.cores=4)
            # FDR
            FDRs <- mcmapply(FDR,simulations.k[[i]]$betas,simulations.k[[i]]$gammas[[j]],mc.cores=4)
            simulations.k[[i]]$FDRavgs[j]  <- mean(unlist(FDRs))       
            # pFDR
            pFDRs <- mcmapply(FDR,simulations.k[[i]]$betas,simulations.k[[i]]$gammas[[j]],TRUE,mc.cores=4)
            simulations.k[[i]]$pFDRavgs[j] <- mean(unlist(pFDRs),na.rm=T)       
        }
        # FDR
        pdf(paste('hw4/pr1_vi_',letters[i],'.pdf',sep=''))
        plot(alphas,simulations.k[[i]]$FDRavgs,type='b',main=bquote("Simulated FDR, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='FDR',ylim=c(0,1))
        dev.off()

        # pFDR
        pdf(paste('hw4/pr1_vii_',letters[i],'.pdf',sep=''))
        plot(alphas,simulations.k[[i]]$pFDRavgs,type='b',main=bquote("Simulated pFDR, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='pFDR',ylim=c(0,1))
        dev.off()


    }
}

### Problem2 
qgammas <- function(pvals,alpha){
    qvals <- qvalue(pvals)$qvalues
    return(qvals < alpha)
}

if (problem2){
    part1 <- FALSE
    part2 <- TRUE
    if (part1) {
        ## i - simulate actual FDR using qvals
        k <- 1000
        scenarios <- list(c(1,1),c(.5,3),c(0,3))    # each is a pi0, sigb pair
        alphas    <- seq(.05,.5,.05)
        FDRs      <- list()
        pFDRs     <- list()
        simulations.k <- list()
        for (i in 1:3){
            simulations.k[[i]] <-list()
            pi0  <- scenarios[[i]][1]
            sigb <- scenarios[[i]][2]
            sims <- mclapply(1:k,function(x){simulateD(pi0,sigb)},mc.cores=4)
            simulations.k[[i]]$p.values <- mclapply(sims,function(x){x$p.values},mc.cores=4)
            simulations.k[[i]]$betas    <- mclapply(sims,function(x){x$beta}    ,mc.cores=4)
            simulations.k[[i]]$gammas<-list()
            simulations.k[[i]]$FDRs  <-NULL
            simulations.k[[i]]$pFDRs <-NULL
            for (j in 1:length(alphas)) {
                simulations.k[[i]]$gammas[[j]] <- mclapply(simulations.k[[i]]$p.values,qgammas,alpha=alphas[j],mc.cores=4)
                # FDR
                FDRs <- mcmapply(FDR,simulations.k[[i]]$betas,simulations.k[[i]]$gammas[[j]],mc.cores=4)
                simulations.k[[i]]$FDRavgs[j]  <- mean(unlist(FDRs))       
                # pFDR
                pFDRs <- mcmapply(FDR,simulations.k[[i]]$betas,simulations.k[[i]]$gammas[[j]],TRUE,mc.cores=4)
                simulations.k[[i]]$pFDRavgs[j] <- mean(unlist(pFDRs),na.rm=T)       
            }
            # FDR
            pdf(paste('hw4/pr2_i_FDR_',letters[i],'.pdf',sep=''))
            plot(alphas,simulations.k[[i]]$FDRavgs,type='b',main=bquote("Simulated FDR using q values, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='FDR',ylim=c(0,1))
            dev.off()
 
            # pFDR
            pdf(paste('hw4/pr2_i_pFDR_',letters[i],'.pdf',sep=''))
            plot(alphas,simulations.k[[i]]$pFDRavgs,type='b',main=bquote("Simulated pFDR using q values, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='pFDR',ylim=c(0,1))
            dev.off()
 
        }
    }

    if (part2) {
        k <- 100
        pi <- seq(0,1,length.out=9)
        sig <- c(.5,1,2,4)
        nicesig <-c('pt5','1','2','4')
        outmat <- matrix(,nrow=length(sig),ncol=length(pi))
        for (i in 1:length(sig)) {
            for (j in 1:length(pi)) {
                pi_gen <- function(x) {
                    return(qvalue(simulateD(pi[j],sig[i])$p.values)$pi0)
                }
                pi0s <- mclapply(1:k,pi_gen,mc.cores=4)
                outmat[i,j] <- mean(unlist(pi0s))
            }
            pdf(paste('hw4/pr2_ii_',nicesig[i],'.pdf',sep=''))
            plot(pi,outmat[i,],type='b',main=bquote("True vs. Estimated "~ pi[0] ~', '~ sigma[b]== .(sig[i])),xlab=expression(pi[0]),ylab=expression(hat(pi)[0]),ylim=c(0,1),xlim=c(0,1))
            abline(a=0,b=1,lty=2)
            dev.off()

        }

        #print(xtable(outmat),type="latex",file="hw4/pr2_ii.tex")

    }

}


### Problem 3
if (problem3){
    ##part iii - r function to compute log-likelihood for hyperparameters

    # compute log-likelihood for single observation in the hierarchical model
    hloglik.sin <- function(Dbar,pi,sig.b,sig.n=1/sqrt(10)) {

        # fraction part second term in log
        frac  <- (1-pi)*sig.n/sqrt(sig.b^2+sig.n^2)
        # exponential part second term in log
        expon <- exp(sig.b^2*Dbar^2/(2*sig.n^4 + 2*sig.b^2)) 

        # actual log-likelihood
        loglik <- log(pi + frac*expon)  
        return(loglik)
    }

    # Vectorize the log-likelihood for a single term
    v.hloglik.sin <- Vectorize(hloglik.sin,c('pi','sig.b'))

    # function to return the entire log-likelihood
    hloglik <- function(Dbar,pi,sig.b) {
       return(sum(v.hloglik.sin(Dbar,pi=pi,sig.b=sig.b)))
    }

    # wrapper of above for optim
    hloglik.optim <- function(param,Dbar){
        # convert to form in equation
        pi    <- inv.logit(param[1])
        sig.b <- exp(param[1])
        return(hloglik(Dbar,pi,sig.b))
    }

    # converts output optim back into normal form
    conv.param <- function(param) {
        return(c(inv.logit(param[1]),exp(param[2])))
    }

    Dbar <- apply(simulateD(.2,2)$D,1,mean)
    opt.param <- conv.param(optim(c(.5,2),hloglik.optim,Dbar=Dbar)$par)

    ## part v - estimate pvalues with posterior
    post.pval <- function(Dbar,pi,sig.b,sig.n=1/sqrt(10)){
        denom <- exp(hloglik.sin(Dbar,pi,sig.b,sig.n=1/sqrt(10)))
        return(pi/denom)
    }
    v.post.pval <- Vectorize(post.pval,vectorize.args=c('pi','sig.b'))

    simulateDbay <- function(pi0,sigb,n=10,m=1000,sigj=1){
        output <- list()
        output$beta     <- ifelse(rbinom(m,1,pi0),0,rnorm(m,0,sigb))
        output$Dbar     <- apply(replicate(n,rnorm(m,output$beta,sigj)),1,mean)
        output$p.values <- v.post.pval(output$Dbar,pi=pi0,sig.b=sigb)
        return(output)
    }

    k <- 1000
    scenarios <- list(c(.1,.5),c(.3,1),c(.4,1))    # each is a pi0, sigb pair
    alphas    <- seq(.05,.5,.05)
    FDRs      <- list()
    pFDRs     <- list()
    simulations.k <- list()
    for (i in 1:3){
        simulations.k[[i]] <-list()
        pi0  <- scenarios[[i]][1]
        sigb <- scenarios[[i]][2]
        sims <- mclapply(1:k,function(x){simulateDbay(pi0,sigb)},mc.cores=4)
        simulations.k[[i]]$p.values <- mclapply(sims,function(x){x$p.values},mc.cores=4)
        simulations.k[[i]]$betas    <- mclapply(sims,function(x){x$beta}    ,mc.cores=4)
        simulations.k[[i]]$gammas<-list()
        simulations.k[[i]]$FDRs  <-NULL
        simulations.k[[i]]$pFDRs <-NULL
        for (j in 1:length(alphas)) {
            simulations.k[[i]]$gammas[[j]] <- mclapply(simulations.k[[i]]$p.values,function(x){x<alphas[j]},mc.cores=4)
            # FDR
            FDRs <- mcmapply(FDR,simulations.k[[i]]$betas,simulations.k[[i]]$gammas[[j]],mc.cores=4)
            simulations.k[[i]]$FDRavgs[j]  <- mean(unlist(FDRs))       
            # pFDR
            pFDRs <- mcmapply(FDR,simulations.k[[i]]$betas,simulations.k[[i]]$gammas[[j]],TRUE,mc.cores=4)
            simulations.k[[i]]$pFDRavgs[j] <- mean(unlist(pFDRs),na.rm=T)       
        }
        # FDR
        pdf(paste('hw4/pr3_v_FDR_',letters[i],'.pdf',sep=''))
        plot(alphas,simulations.k[[i]]$FDRavgs,type='b',main=bquote("Simulated FDR using Empirical Bayes, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='FDR',ylim=c(0,1))
        dev.off()
 
        # pFDR
        pdf(paste('hw4/pr3_v_pFDR_',letters[i],'.pdf',sep=''))
        plot(alphas,simulations.k[[i]]$pFDRavgs,type='b',main=bquote("Simulated pFDR using Empirical Bayes, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='pFDR',ylim=c(0,1))
        dev.off()
 
    }
    

}

