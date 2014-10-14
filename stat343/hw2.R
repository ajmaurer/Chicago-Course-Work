wine<-read.csv(file="hw1/winequalityRed.csv",head=TRUE,sep=";")
lm<-lm(pH ~ fixed.acidity, data=wine)

# Generate confidence intervals
confint(lm,level=.66)
confint(lm,level=.9)
confint(lm,level=.95)

# Generate variance estimate
c<-c(1,10)
n<-nrow(wine)
cov<-summary(lm)$cov.unscaled
sig2 <- summary(lm)$sigma^2
varpred<- sig2 * t(c) %*% cov %*% c 
varpred

# Make prediction
predict(lm,newdata=data.frame(fixed.acidity=13.1), interval="confidence")
predict(lm,newdata=data.frame(fixed.acidity=13.1), interval="prediction")

# Make and plot 100 predictions
points<-data.frame(fixed.acidity =seq(min(wine$fixed.acidity),max(wine$fixed.acidity),length.out=100))
predictions<-predict(lm,newdata=points, interval="confidence")
pdf('hw2_plot.pdf')
plot(points$fixed.acidity,predictions[,1],cex=.5,main="CI of Mean pH Estimate", xlab="Fixed Acidity", ylab="pH")
points(points$fixed.acidity,predictions[,2],cex=.5, pch=18)
points(points$fixed.acidity,predictions[,3],cex=.5, pch=18)
dev.off()
