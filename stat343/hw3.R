library('faraway')

#############################
# Problem 1
#############################

# Load the data and create the matrices 
data('savings')
attach(savings)
n<-nrow(savings)
x<-as.matrix(data.frame(rep(1,n),savings[-1]))
y<-as.matrix(savings[1])

### Part c)
# Loop over epsilons and try methods until failure
eps <- .1
norm.diff <- NULL
continue <- TRUE
while (continue) {
    eps <- eps *.1
    # Add pop75 + noise
    x.i<-cbind(x,as.matrix(x[,3]+ eps*rnorm(n))) 

    #See if the QR method still works
    QR.success <-tryCatch(
        {
            q<-qr.Q(qr(x.i))
            r<-qr.R(qr(x.i))
            beta.qr <- solve(r, t(q) %*% y)
            TRUE
        },
        error=function(e) return(FALSE),
        warning=function(e) return(FALSE)
        )

    # See if the normal method still works
    Normal.success <-tryCatch(
        {
            beta.x.i  <- solve(t(x.i)%*%x.i) %*% t(x.i) %*% y
            TRUE
        },
        error=function(e) return(FALSE),
        warning=function(e) return(FALSE)
        )

    # Keep going and calculate norm of difference if both succeed 
    continue <- QR.success & Normal.success
    if (continue) norm.diff <- c(norm.diff,norm(beta.qr-beta.x.i))
}
if (!QR.success) print(paste('The QR method failed at an epsilon of ',toString(eps)))
if (!Normal.success) print(paste('The normal equation failed at an epsilon of ',toString(eps)))

### Part d)
model<-lm(y~x.i)
summary(model)
beta.qr

#############################
# Problem 5
#############################

data(prostate)

### 3.1.a)
m<-lm(lpsa~lcavol + lweight + age +lbph + svi + lcp + gleason + pgg45, data=prostate)
coe<-summary(m)$coef
ci.9 <- c(coe[4,1] + qt(.05,88)*coe[4,2],coe[4,1] - qt(.05,88)*coe[4,2])
ci.95 <- c(coe[4,1] + qt(.025,88)*coe[4,2],coe[4,1] - qt(.025,88)*coe[4,2])

### 3.1.c)
nb<-1000
n<-nrow(prostate)
x<-as.matrix(prostate[1:8])
y<-as.matrix(prostate[9])

# Thus function pulls the t score from a regression predicting the y permuted
est.t.age <- function(n) summary(lm(y[sample(n,rep=FALSE)]~x))$coef[4,3]

# Which we now run on a vector of the draw size
t.age <- sapply(rep(n,nb),est.t.age)
pval<- mean(abs(t.age)>abs(coe[4,3]))

### 3.1.d)
#only lcavol, lweight, and svi had a pvalue of less than .05
rm<-lm(lpsa~lcavol + lweight + svi, data=prostate)
anova(m,rm)

### 3.2.a)
data(cheddar)
summary(lm(taste ~ Acetic + H2S + Lactic, data=cheddar))

### 3.2.b)
summary(lm(cheddar$taste ~ exp(cheddar$Acetic) + exp(cheddar$H2S) + cheddar$Lactic))

