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
if (problem1) {
    ## part (b) - Write a function to calculate credible intervals
    drawTheta <- function(b,X,n,a) {
        p <- rbeta(3,a+X,a+n-X)
        return(p[1]*p[2]/(p[1]*p[2]+(1-p[1])*p[3]))
    }
    bootThetaCI <- function(b,X,n,a) {
        vals<-sort(unlist(mclapply(1:b,drawTheta,X<-X,n<-n,a<-a,mc.cores=4)))
        return(c(vals[ceiling(b*.025)],vals[ceiling(b*.975)]))
    }

    ## part (c)/(d) - generate credible intervals for a series of inputs

    formList <- function(li) {
        return(paste('(',paste(li,collapse=','),')',sep=''))
    }
    formCI <- function(ci,d=3) {
        return(formList(format(ci,digits=d)))
    }

    scenarios <- list( n=list(20,20,80,80),X=list(c(2,18,2),c(10,18,0),c(20,60,20),c(40,72,8)),CIpt5=list(),CI1=list(),CI2=list())    

    for (i in 1:4) {
        scenarios$CIpt5[[i]] <- bootThetaCI(10000,scenarios$X[[i]],scenarios$n[[i]],.5)
        scenarios$CI1[[i]]   <- bootThetaCI(10000,scenarios$X[[i]],scenarios$n[[i]],1)
        scenarios$CI2[[i]]   <- bootThetaCI(10000,scenarios$X[[i]],scenarios$n[[i]],2)
    }

    # build the tables
    table.1.cd = data.frame(
        n=unlist(scenarios$n),
        x=unlist(lapply(scenarios$X,formList)),
        CIpt5=unlist(lapply(scenarios$CIpt5,formCI)),
        CI1  =unlist(lapply(scenarios$CI1,formCI)),
        CI2  =unlist(lapply(scenarios$CIpt5,formCI))
        )
    colnames(table.1.cd) <- c("n0=n1=n2","(x0,x1,x2)","Alpha=.5 95% CI","Alpha=1 95% CI","Alpha=2 95% CI")
    
    print(xtable(table.1.cd),type="latex",file="hw3/pr1_cd.tex")

    ## part (e) - simulate the frequentist coverage properties
    simFreqCov <- function(i,b,n,p,a) {
        X <- rbinom(3,n,p)
        th <- p[1]*p[2]/(p[1]*p[2]+(1-p[1])*p[3]) 
        CI <- bootThetaCI(b,X,n,a)
        return(ifelse(th<CI[1],-1,ifelse(th>CI[2],1,0)))
    }

    # parameters
    n <- 20
    b <- 1000
    ps <- list(c(.5,.5,.5),c(.2,.6,.7),c(.5,.1,.9),c(.95,.95,.05),c(.2,.1,.9))
    as <- list(.5,1,2)
    leftTail <- matrix(,nrow=length(ps),ncol=length(as))
    rightTail <- matrix(,nrow=length(ps),ncol=length(as))
    for (i in 1:length(ps)) {
        for (j in 1:length(as)) {
            outcome <- unlist(lapply(1:b,simFreqCov,b<-b,n<-n,p<-ps[[i]],a<-as[[j]]))
            leftTail[i,j] <- sum(outcome==-1)/b
            rightTail[i,j] <- sum(outcome==1)/b
        }
    }
    coverageData <- data.frame(
        p=unlist(lapply(ps,formList)),
        CIpt5Low = leftTail[,1],
        CIpt5Hi = rightTail[,1],
        CI1Low = leftTail[,2],
        CI1Hi = rightTail[,2],
        CI2Low = leftTail[,3],
        CI2Hi = rightTail[,3]
        )
    print(xtable(coverageData),type="latex",file="hw3/pr1_e.tex")
}


