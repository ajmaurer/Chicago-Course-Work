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

## part d
A <- as.factor(c(3,2,2,3,1,2,3,1,1))
B <- as.factor(c(1,1,1,2,2,2,3,3,3))
y <- c(6.74,6.34,5.57,10.86,.23,4.73,8.57,3.05,2.18)
data.1.d<-data.frame(A<-A,B<-B,y<-y)

# part iii
models <- list(lm(y~1,data=data.1.d)) 
models[[2]] <- lm(y~A,data=data.1.d)
models[[3]] <- lm(y~B,data=data.1.d)
models[[4]] <- lm(y~A+B,data=data.1.d)
models[[5]] <- lm(y~A*B,data=data.1.d)

output.1.d <- matrix(,nrow=5,ncol=3)
rownames(output.1.d) <- c("1","A","B","A+B","A*B")
colnames(output.1.d) <- c("Dimension","Sqrd Len Proj","Sqrd Len Comp")
for (m in 1:5) {
    output.1.d[m,1] <- models[[m]]$rank
    output.1.d[m,2] <- sum(predict(models[[m]])^2)
    output.1.d[m,3] <- sum(summary(models[[m]])$residual^2)
}

print(xtable(output.1.d,digits=c(0,0,2,2)),type="latex",file="hw1/pr1_d_iii.tex")

### Problem 2
fungus <- read.table("http://www.stat.uchicago.edu/~pmcc/datasets/fungus.txt",header=T)
fungus <- melt(fungus,id.vars="Tmp",value.name="logwt",variable.name="Isolate")
fungus$Tmp_f <- as.factor(fungus$Tmp)

## part i
anova.2.i <-anova(lm(logwt~Isolate+Tmp_f,data=fungus))
rownames(anova.2.i)<-c("Isolate","Temp","Resid")
print(xtable(anova.2.i,digits=c(0,0,4,4,4,4)),type="latex",file="hw1/pr2_i.tex")

## part ii
# The problem requires stitching together some godawful version of an anova table
anova.2.ii.main <-anova(lm(logwt~poly(Tmp,1)+poly(Tmp,2)+poly(Tmp,3)+poly(Tmp,4)+Tmp_f,data=fungus))[1:5,1:3]
anova.2.ii.sub  <-anova(lm(logwt~poly(Tmp,3)+Tmp_f,data=fungus))[2,1:3]

anova.2.ii <- rbind(anova.2.ii.main,anova.2.ii.sub,anova.2.i[2:3,1:3])
rownames(anova.2.ii) <- c("P1/1","P2/P1","P3/P2","P4/P3","Temp/P4","Temp/P3","Total Temp","Resid")
F <- c(anova.2.ii[1:7,3]/anova.2.ii[8,3],NA)
pvalue<-rep(NA,8)
for (i in 1:8) {
    pvalue[i] <- signif(1 - pf(F[i],anova.2.ii[i,1],anova.2.ii[8,1]),4)
}

anova.2.ii <- cbind(anova.2.ii,F,pvalue)
rownames(anova.2.ii)[4:5]<-c("F value","Pr(>F)")
print(xtable(anova.2.ii,digits=c(0,0,4,4,4,4)),type="latex",file="hw1/pr2_ii.tex")

## part iv
poly.2.iv <- poly(fungus$Tmp,3)
beta.2.iv <- coef(lm(logwt~poly.2.iv+Isolate,data=fungus))[2:4]
temp.2.iv <- seq(50,90,length.out=5000)
curve.2.iv <- predict(poly.2.iv,temp.2.iv) %*% beta.2.iv

opt.temp <- temp.2.iv[sort(curve.2.iv,decreasing=T,index.return=T)$ix[1]]

pdf("hw1/pr2_iv.pdf")
plot(temp.2.iv,curve.2.iv,type="l",main="3rd Degree Polynomial on Temperature",xlab="Temperature",ylab="Predicted Effect on Growth")
dev.off()

### Problem 3
elem <- c("Mg","Zn","Fe","Pb","Cu")
elec <- c("O","A","K")
elect<-read.table("http://www.stat.uchicago.edu/~pmcc/glm/electro_chem.dat")
colnames(elect)<-elem
elect$row<-elem

elect$electrolyte <- rep(elec,each=5)
elect_long <- melt(elect,id.vars=c("electrolyte","row"),value.name="voltage",variable.name="col")
elect_long$row <- factor(elect_long$row, levels=elem)
elect_long$col <- factor(elect_long$col, levels=elem)

## part ii
alpha <- matrix(,nrow=3,ncol=4)
include <- elem[2:5]
for (i in 1:3) {
    mm<-model.matrix(voltage~row+col+0,data=elect_long[elect_long$electrolyte==elec[i],])
    X <- mm[,paste("row",include,sep="")] - mm[,paste("col",include,sep="")] 
    alpha[i,] <- coef(lm(elect_list[[i]]$voltage ~ X+0))
}

colnames(alpha) <- include
rownames(alpha) <- names(elect_list)

print(xtable(alpha,digits=4),type="latex",file="hw1/pr3_ii.tex")

## part iii
mm<-model.matrix(voltage~electrolyte+row+col+0,data=elect_long)
X <- mm[,paste("row",include,sep="")] - mm[,paste("col",include,sep="")]

lm.3.iii.1 <- lm(voltage~electrolyte:X,data=elect_long)
lm.3.iii.2 <- lm(voltage~X,data=elect_long)

print(xtable(anova(lm.3.iii.2,lm.3.iii.1)),type="latex",file="hw1/pr3_iii.tex")

### Problem 5
temp <- as.factor(rep(c(0,10,20),each=4))
weeks <- rep(c(2,4,6,8),3)
ascor <- c(45,47,46,46, 45,43,41,37, 34,28,21,16)
beans <- data.frame(temp,weeks,ascor)

lm.5 <- lm(log(ascor) ~ temp:weeks,data=beans)

se<-summary(lm.5)$coefficients[,2]
interval<-(cbind(coef(lm.5)+qt(.025,8)*se,coef(lm.5),coef(lm.5)-qt(.025,8)*se)/log(.5))^-1
