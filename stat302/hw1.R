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
problem2 <- FALSE 
problem5 <- TRUE

### Problem 2
# We are simulating a number of p-values in the context of experiments where treatment effect is estimated as compared to a control
# with probability p, treatment has effect beta~exp(l), otherwise beta=0
# n units receive treatment, n receive control
# If treatment group, y~N(beta,1). Thus, beta_hat ~ N(beta,n^-.5) and sigma_hat ~ chisq(2n-2)/(2n-2) (these two statistics are independent random variables, even when generated from the same data, so can be drawn independently).
# Thus, beta_hat/sqrt(sigma_hat/n) ~ t_(2n-2), giving p-values
if (problem2) {
    ps  <- c(.1,.3,.5,.7,.9)    # portion positive treatment effect (pi)
    ls  <- c(.125,.5,2,8,16,32) # defines exponential which beta is drawn from
    n  <- 10                    # sample units receive control/treatment each
    boots <- 100000

    # Generates pvalue for each experiment. 1st variable indicates positive effect,
    # 2nd parameter for exponential, 3rd number experimental units
    pval_boot <- function(true,lambda,size) {
        beta      <- ifelse(true,rexp(1,lambda),0)
        beta_hat  <- rnorm(1,beta,n^-.5)
        sigma_hat <- rchisq(1,2*n-2)/(2*n-2)
        tscore    <- beta_hat/sqrt(sigma_hat/n)
        return(pt(tscore,2*n-2,lower.tail=F))
    }

    # For each compination of pi and lambda, generate portion of p-values between 
    # .04 and .05 which are true positives over given number boots
    nl<-length(ls)
    np<-length(ps)
    evidence<-matrix(,nrow=nl,ncol=np)
    for (i in 1:nl) {
        for (j in 1:np) {
            l <- ls[i]
            p <- ps[j]
            trues <- runif(boots)<p
            pvals <- unlist(mclapply(trues,pval_boot,lambda=l,size=n,mc.cores=4))
            evidence[i,j] <- mean(trues[.04<=pvals & pvals<=.05])
        }
    }

    rownames(evidence)<-paste("lambda =",ls)
    colnames(evidence)<-paste("pi =",ps)
    print(xtable(evidence,digits=2),type="latex",file="hw1/pr2.tex")
}

### Problem 5
if (problem5) {
    subpops <- c("EelR","FeatherHfa","FeatherHsp","KlamathRfa")
    # Call given code
    source('official/exercises/seeb/train_test.R')

    # Function to summarize frequency at locus by population (similar to trainc)
    # Add 1 so that no Allele is impossible in each subpopulation
    compute_freq <- function(data,locus){
        counts <- table(data[,1+2*locus],data$Population) + 
                    table(data[,2+2*locus],data$Population) + 1
        return(counts/sum(counts))
    }

    # Get frequency at each locus
    train_freq <- list()
    for (i in 1:12) {
        train_freq[[i]] <- as.data.frame.matrix(compute_freq(train,i))
        train_freq[[i]]$allele <- as.factor(rownames(test_freq[[i]]))
    }

    # Set uniform prior
    priors<-paste(subpops,"prior",sep="_") 
    test[,priors]<-.25

    # Calculate Log-Likelihood
    log_lks <-paste(subpops,"loglk",sep="_")
    test[,log_lks]<-0
    for (i in 1:12) {
        for (j in 1:2) {
            test <- merge(test,train_freq[[i]],by.x=names(test)[j+i*2],
                          by.y="allele")
            # Missing data won't effect likelihood
            test[is.na(test[,j+i*2]),subpops] <- 1  
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

    # Accuracy
    accuracy <- mean(test$predict==test$Population)
        
}
        
        
        
 

