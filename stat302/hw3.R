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
problem1 <- FALSE
problem5 <- TRUE

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
        CI2  =unlist(lapply(scenarios$CI2,formCI))
        )
    colnames(table.1.cd) <- c("n0=n1=n2","(x0,x1,x2)","a=.5 95% CI","a=1 95% CI","a=2 95% CI")
    
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
    colnames(coverageData) <- c("(p0,p1,p2)",paste("a=",rep(c(.5,1,2),times=1,each=2),c(" th<LB"," UB<th"),sep=''))
    
    print(xtable(coverageData,digits=3),type="latex",file="hw3/pr1_e.tex")
}

### Problem 5
if (problem5) {

    # Call given code
    source('official/exercises/seeb/train_test.R')
    subpops <- unique(test$Population)
    loci <- names(test)[3:26]
    loci <- loci[!grepl("\\.1",loci)] # get rid of .1 version
    alpha=2

    ## part (a)
    # Function to summarize frequency at locus by population (similar to trainc)
    # We assume a dirchlet prior with each alpha identical
    compute_freq <- function(data,locus,alpha){
        counts <- table(data[,locus],data$Population) + 
                    table(data[,paste(locus,".1",sep="")],data$Population) + alpha
        return(counts/colSums(counts))
    }

    calculate_post <-function(alpha,train,test,subpops,loci){
        # Get frequency at each locus
        train_freq <- list()
        for (loc in loci) {
            train_freq[[loc]] <- as.data.frame.matrix(compute_freq(train,loc,alpha))
            train_freq[[loc]]$allele <- as.factor(rownames(train_freq[[loc]]))
        }
 
        # Set uniform prior
        priors<-paste(subpops,"prior",sep="_") 
        test[,priors]<-.25
 
        # Calculate Log-Likelihood
        log_lks <-paste(subpops,"loglk",sep="_")
        test[,log_lks]<-0
        for (loc in loci) {
            for (e in c("",".1")) {
                locus <- paste(loc,e,sep="")
                test <- merge(test,train_freq[[loc]],by.x=locus,
                              by.y="allele",all.x=T,sort=F)
                # Missing data won't effect likelihood
                test[is.na(test[,locus]),subpops] <- 1  
                test[,log_lks] <- test[,log_lks] + log(test[,subpops])
                test<-test[,!(names(test) %in% subpops)]
 
            }
        }
 
        # Calculate Posterior
        posteriors<-paste(subpops,"post",sep="_")
        test[,posteriors] <- test[,priors]*exp(test[,log_lks]) 
        test[,posteriors] <- test[,posteriors]/apply(test[,posteriors],1,sum) # normalize

        # Predict population with highest posterior
        test$post_max <- apply(test[,posteriors],1,max)
        test$predict <- ""
        for (pop in subpops) {
            test[test$post_max==test[,paste(pop,"post",sep="_")],"predict"] <- pop
        }

        # Calculate total log likelihood
        tot_log_lk <- sum(log(apply(exp(test[,log_lks]),1,sum)))

        # Accuracy
        error_rate <- mean(test$predict!=test$Population)

        return(list(tot_log_lk=tot_log_lk,error_rate=error_rate,dataset=test))
    }

    ## part(b) 
    alist <- c(.001,.01,.2,.5,1,2,5,10,100)
    alist_error_rate<-mclapply(alist,function(a){calculate_post(a,train,test,subpops,loci)$error_rate},mc.cores=4)

    table.5.b <- rbind(alist,format(alist_error_rate,digits=3))
    rownames(table.5.b) <- c("Alpha Values","Error Rate")
    colnames(table.5.b) <- paste("Value",1:length(alist))

    print(xtable(table.5.b),type="latex",file="hw3/pr5_b.tex")

    ## part(d) 
    opt_param <- optim(1,function(a){-calculate_post(a,train,test,subpops,loci)$tot_log_lk},method="L-BFGS-B",lower=1E-10,upper=Inf)

    ## part(e) - we assume a discrete uniform prior on the hyperparameter a
    a_prior <- c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4)
    len <- length(a_prior)

    # few groups of  variable names
    log_lks <-paste(subpops,"loglk",sep="_")
    tot_lks <-paste(subpops,"tot_lk",sep="_")
    posteriors<-paste(subpops,"post",sep="_")

    # Generate data for each prior and copy test
    a_datasets <- mclapply(a_prior,function(a){calculate_post(a,train,test,subpops,loci)$dataset},mc.cores=4)
    test.e <-test
    test.e[,tot_lks] <- 0

    # Set uniform prior
    priors<-paste(subpops,"prior",sep="_") 
    test.e[,priors]<-.25

    # Calculate total likelihood
    for (i in 1:len) {
        test.e <- merge(test.e,a_datasets[[i]][,c('Individual',log_lks)],by.x="Individual",by.y="Individual",all.x=T,sort=F)
        test.e[,tot_lks] <- test.e[,tot_lks] + exp(test.e[,log_lks])/len
        test.e<-test.e[,!(names(test.e) %in% log_lks)]
    }

    # Then the posteriors
    test.e[,posteriors] <- test.e[,priors]*test.e[,tot_lks] 
    test.e[,posteriors] <- test.e[,posteriors]/apply(test.e[,posteriors],1,sum) # normalize

    # Predict population with highest posterior
    test.e$post_max <- apply(test.e[,posteriors],1,max)
    test.e$predict <- ""
    for (pop in subpops) {
        test.e[test.e$post_max==test.e[,paste(pop,"post",sep="_")],"predict"] <- pop
    }
 
    # Accuracy
    error_rate.e <- mean(test.e$predict!=test.e$Population)

}



