library("alr3")
library("FNN")
data("htwt")
sink("hw1.Rout", append=FALSE, split=TRUE)

readline("Problem 2.1.1")
jpeg('hw1_2_1_1.jpg')
plot(htwt$Ht,htwt$Wt)
dev.off()
readline("I would argue that this data, though seeming to indicate some sort of correlation between height and weight befiting a linear model, is too small considering the amount of noise for a linear regression to make sense. The estimates would be far too noisy to read anything into it. That said, if there was more data it looks like a linear model would make sense.")

readline("Problem 2.1.2")
n<-length(htwt$Ht)
xbar<-sum(htwt$Ht)/n
xbar
ybar<-sum(htwt$Wt)/n
ybar
sxx<-sum((htwt$Ht-xbar)^2)
sxx
syy<-sum((htwt$Wt-ybar)^2)
syy
sxy<-sum((htwt$Ht-xbar)*(htwt$Wt-ybar))
sxy
beta1 <- sxy/sxx
beta1
beta0 <- ybar - beta1*xbar
beta0
jpeg('hw1_2_1_2.jpg')
plot(htwt$Ht,htwt$Wt)
curve(beta0 + x*beta1, add=TRUE)
dev.off()

readline("Problem 2.1.3")
rss<-syy-sxy^2/sxx
sigmasq<-rss/(n-2)
sigmasq
se_beta1 <- sqrt(sigmasq/sxx)
se_beta1
se_beta0 <- sqrt(sigmasq*(1/n + xbar^2/sxx))
se_beta0
cov_betas <- -sigmasq*xbar/sxx
cov_betas
t_beta1 <- beta1/se_beta1
t_beta1
p_val_beta1 <- 2*pt(-abs(t_beta1),n-2)
p_val_beta1
t_beta0 <- beta0/se_beta1
t_beta0
p_val_beta0 <- 2*pt(-abs(t_beta0),n-2)
p_val_beta0

readline("Problem 2.1.4")
df<-c(1,n-2,n-1)
ss<-c(syy-rss,rss,syy)
ms<-c((syy-rss)/1,sigmasq,NA)
f<-c((syy-rss)/1/sigmasq,NA,NA)
pval<-pf(-abs(f),1,n-2)
anova<-data.frame(df,ss,ms,f,pval, row.names = c("Regression","Residual","Total"))
anova
t_beta1^2

readline("Problem 4.a")
wine<-read.csv(file="winequalityRed.csv",head=TRUE,sep=";")
summary(wine)
plot(wine)

jpeg('hw1_4b_alcohol.jpg')
hist(wine$alcohol)
dev.off()
jpeg('hw1_4b_ph.jpg')
hist(wine$pH)
dev.off()
jpeg('hw1_4den_to_acid.jpg')
plot(wine$density,wine$fixed.acidity)
dev.off()
jpeg('hw1_4alc_to_cit.jpg')
plot(wine$alcohol,wine$citric.acid)
dev.off()

readline("Problem 4c")
subset<-subset(wine, select=c(fixed.acidity,pH,citric.acid,density))
jpeg('hw1.4c.jpg')
plot(subset)
dev.off()
test<-rnorm(length(subset$pH))
jpeg('hw1_4c_fixed_acidity.jpg')
qqplot(test,wine$fixed.acidity)
dev.off()
jpeg('hw1.4c.pH.jpg')
qqplot(test,wine$pH)
dev.off()
jpeg('hw1_4c_citric_acid.jpg')
qqplot(test,wine$citric.acid)
dev.off()
jpeg('hw1_4c_density.jpg')
qqplot(test,wine$density)
dev.off()

readline("Problem 4.d")
phlm<-lm(wine$pH ~ wine$citric.acid)
summary(phlm)
jpeg('hw1_4d_pHresiduals.jpg')
plot(fitted(phlm),residuals(phlm))
dev.off()
denlm<-lm(wine$density ~ wine$fixed.acidity)
summary(denlm)
jpeg('hw1_4d_denresiduals.jpg')
plot(fitted(denlm),residuals(denlm))
dev.off()

readline("Problem 4.e")
phnnm<- knn.reg(wine$citric.acid, test=NULL, wine$pH, k=50, algorithm=c("kd_tree","cover_tree","brute"))
jpeg('hw1_4e_pHfit.jpg')
plot(wine$citric.acid,wine$pH, pch=3)
points(wine$citric.acid,phlm$fitted.values, pch=18)
points(wine$citric.acid,phnnm$pred, pch=20)
dev.off()
dennnm<- knn.reg(wine$fixed.acidity, test=NULL, wine$density, k=50, algorithm=c("kd_tree","cover_tree","brute"))
jpeg('hw1_4e_denfit.jpg')
plot(wine$fixed.acidity,wine$density, pch=3)
points(wine$fixed.acidity,denlm$fitted.values, pch=18)
points(wine$fixed.acidity,dennnm$pred, pch=20)
dev.off()

sink()
