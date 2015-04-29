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
problem2 <- TRUE 

### Problem 2
if (problem2) {
    ## part iii)

    # compute numerical integral
    integrand <- function(th){ pnorm(th-2.32)*dnorm(th) }
    numerator <- integrate(integrand,lower<- -Inf, upper<- 0)$value
    denominator <- integrate(integrand,lower<- -Inf, upper<- Inf)$value
    portion_int <- numerator/denominator

    # Simulation
    n <- 100000
    th <- rnorm(n)
    X <- rnorm(n,mean=th)
    portion_sim <- sum(X>2.32 & th<0)/sum(X>2.32)


}


