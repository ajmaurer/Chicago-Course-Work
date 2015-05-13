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

### Caller
problem1 <- TRUE

### Problem 1
if (problem1){
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
        pdf(paste('hw4/pr1_vi_hist_',letters[i],'.pdf',sep=''))
        plot(alphas,simulations.k[[i]]$FDRavgs,type='b',main=bquote("Simulated FDR, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='FDR',ylim=c(0,1))
        dev.off()

        # pFDR
        pdf(paste('hw4/pr1_vii_hist_',letters[i],'.pdf',sep=''))
        plot(alphas,simulations.k[[i]]$pFDRavgs,type='b',main=bquote("Simulated pFDR, "~pi[0] == .(pi0)~','~ sigma[b]== .(sigb)),xlab=expression(alpha),ylab='pFDR',ylim=c(0,1))
        dev.off()

    }



            





}


