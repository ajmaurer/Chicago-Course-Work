library("alr3")

###########################
# Problem 4
###########################
data(wm1)
attach(wm1)

### 2.13
# 1
pdf("4.2.13.1.pdf")
plot(RSpd,CSpd)
dev.off
# 2
reg<-lm(CSpd~RSpd)
summary(reg)
# 3
predict(reg,data.frame(CSpd=NA,RSpd=7,4285), interval="confidence")
