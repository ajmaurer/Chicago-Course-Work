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

### callers
problem0 <- FALSE 
problem3 <- TRUE

### A0
if (problem0) {
    #Mating counts used in the PNAS 2010 paper by
    #Sharon, Segal, Ringo, Hefetz, Zilber-Rosenberg and Rosenberg
    #"Commensal bacteria play a role in mating preference of Drosophila melanogaster"
    #Data provided by the authors on request
    # First column is generation no 
    # Last four columns are mating counts by type CxC CxS SxC SxS
    # Last three rows taken from a separate expt of a similar design)
    matings <- matrix(c(
    2, 12, 8, 9, 16,
    6, 10, 5, 9, 10,
    7, 17, 9, 9, 15,
    9, 8, 7, 7, 9,
    10, 18, 13, 5, 12,
    11, 12, 5, 7, 14,
    13, 14, 9, 8, 12,
    15, 18, 9, 7, 15,
    16, 14, 5, 5, 10,
    17, 31, 22, 12, 27,
    20, 23, 13, 10, 20,
    21, 13, 7, 5, 14,
    26, 30, 19, 12, 21,
    31, 9, 7, 3, 10,
    37, 20, 14, 11, 17,
    11, 18, 11, 7, 16,
    12, 16, 11, 8, 15,
    13, 22, 13, 8, 13), ncol=5, byrow=T)
    colnames(matings) <- c("gen", "CxC", "CxS", "SxC", "SxS")
    
    wells <- c(NA,NA,NA,24, 39, 20, 24, 36, 23, 70, 46, 24, 45, 23, 48, 48, 48, 48)
    totmat <- apply(matings[,2:5],1,sum)
    SII <- (matings[,2]+matings[,5]-matings[,3]-matings[,4])/totmat

    matings <- data.frame(cbind(wells,matings,totmat,SII))

    genpoly = poly(matings$gen,2)
    genlm <- lm(matings$SII~genpoly)
    xseq <- seq(0,40,length.out=500)
    polyseq <- predict(genpoly,xseq) %*% coef(genlm)[2:3] + coef(genlm)[1]

    pdf('hw2/A0_SII.pdf')
    plot(SII~gen,data=matings,main='SII by Generation With Quadratic Fit',xlab='Generation',ylab='SII')
    points(xseq,polyseq,type='l')
    dev.off()

    pdf('hw2/A0_wells.pdf')
    plot(totmat~wells,data=matings[1:15,],main='Matings vs. Wells',xlab='Number of Wells',ylab='Total Number of Matings')
    welllm = lm(totmat~wells,data=matings[1:15,]) 
    abline(welllm)
    points(matings$wells[16:18],matings$totmat[16:18],col=2)
    legend('topleft',lty=c(0,0,1),pch=c(1,1,NA),col=c(1,2,1),legend=c('First 15','Last 3','LM First 15'))
    dev.off()

    X = model.matrix(welllm)
    x0 = c(1,48)
    se = summary(welllm)$sigma * sqrt(1+t(x0) %*% solve(t(X) %*% X) %*% x0)
    tval = sqrt(3)*(x0 %*% coef(welllm) - mean(matings$totmat[16:18]))/se
    pval = 2*pt(-tval,10)
}

### A3
if (problem3) {
    dlog2  <- seq(0,5)
    deaths <- c(0,2,3,5,7,10)
    n      <- c(7,9,8,7,9,11)
    mor.logit <- glm(cbind(deaths,n-deaths)~dlog2,family=binomial)

    xseq = seq(0,10,length.out=500)
    pdf('hw2/A3_plot1.pdf')
    plot(deaths/n~dlog2,main="Mortality by Log Dosage",xlab=expression('Log'[2]*'Dose'),ylab='Mortality Rate')
    points(xseq,predict(mor.logit,data.frame(dlog2=xseq),type='response'),type='l')
    dev.off()

    print(xtable(summary(mor.logit)),type="latex",file="hw2/A3_logit.tex")
    a <- coef(mor.logit)[1]
    b <- coef(mor.logit)[2]
    sa <- vcov(mor.logit)[1,1]
    sab <- vcov(mor.logit)[1,2]
    sb <- vcov(mor.logit)[2,2]
    l1 <- ((a*b - 1.96^2*sab) - sqrt((a*b - 1.96^2*sab)^2-(b^2-1.96^2*sb)*(a^2-1.96^2*sa)))/(b^2-1.96^2*sb)
    l2 <- ((a*b - 1.96^2*sab) + sqrt((a*b - 1.96^2*sab)^2-(b^2-1.96^2*sb)*(a^2-1.96^2*sa)))/(b^2-1.96^2*sb)

    mor.logit4 <- glm(cbind(deaths,n-deaths)~I(dlog2-4)-1,family=binomial)
    print(xtable(summary(mor.logit4)),type="latex",file="hw2/A3_logit4.tex")

    lr4 <- 2*logLik(mor.logit) - 2*logLik(mor.logit4)
    lr4p <- 1-pchisq(lr4,1)

    xseqlk <- seq(0,6,length.out=100)
    loglks <- unlist(lapply(xseqlk,function(x) logLik(glm(cbind(deaths,n-deaths)~I(dlog2-x)-1,family=binomial))))
    cutoff <- logLik(mor.logit) - .5*qchisq(.95,1)
    lblk   <- min(xseqlk[loglks>cutoff])
    ublk   <- max(xseqlk[loglks>cutoff])

    pdf('hw2/A3_lkl.pdf')
    plot(loglks~xseqlk,type='l',xlab=expression('Hypothesized log'[2]*'LD'[50]),ylab="Restricted Log Likelihood")
    abline(h=cutoff,lty=2)
    dev.off()



}


