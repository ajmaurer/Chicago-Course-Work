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

### Problem 1
data.1   <- melt(read.table('hw2/pr1.dat',header=T),id.vars="T",variable='t',value.name='proth') 
data.1$t <- as.numeric(gsub("X","",data.1$t))
data.1   <- data.1[!is.na(data.1$proth),]
data.1$R <- data.1$T - data.1$t

data.1$tfac <- as.factor(data.1$t)
data.1$Tfac <- as.factor(data.1$T)
data.1$Rfac <- as.factor(data.1$R)

an.T.R.t <- anova(lm(proth~Tfac+Rfac+tfac+tfac:Rfac,data=data.1))[,c(2,1,3)]
an.t.T.R <- anova(lm(proth~tfac+Tfac+Rfac+tfac:Rfac,data=data.1))[,c(2,1,3)]
an.R.t.T <- anova(lm(proth~Rfac+tfac+Tfac+tfac:Rfac,data=data.1))[,c(2,1,3)]

output.1     <- an.T.R.t[1,]
output.1[2,] <- an.t.T.R[1,]
output.1[3,] <- an.R.t.T[1,]
output.1[4,] <- an.t.R.T[3,]
output.1[5,] <- an.T.R.t[3,]
output.1[6,] <- an.t.T.R[3,]
output.1[7,] <- an.t.T.R[4,]

vars <- c('T','t','R')
numer <- '(Fac(R)+Fac(t)+Fac(T))'
denoms <- c('(Fac(t)+Fac(R))','(Fac(R)+Fac(T))','(Fac(t)+Fac(T))')

rownames(output.1) <- c(paste('Fac(',vars,')/1',sep=''),paste(numer,'/',denoms,sep=''),paste('(Fac(t):Fac(R)/',numer,sep=''))

print(xtable(output.1),type="latex",file="hw2/pr1.tex")

      
      
      
      
      
      
      
      
      


