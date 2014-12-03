library('parallel')
library('alr3')
library('faraway')
library('car')
library('MASS')
library('locfit')
library('xtable')
library('nlme')
library('quantreg')
library('splines')
library('rpart')
library('randomForest')
library('leaps')

prob3<-F
prob4<-F
prob5<-T

d<-3

if (prob3) {
    data(teengamb)

    ### a)
    crit<-.2
    exclude<-NULL
    excl.pval<-NULL
    go<-TRUE
    p.lm<-lm(gamble~.,data=teengamb)

    while (go) {
        p.lm<-lm(gamble~.,data=teengamb[!(names(teengamb) %in% exclude)])
        pvals<-summary(p.lm)$coefficients[-1,4]
        max.ind <- which.max(pvals)
        max.pval <- pvals[max.ind]
        if (max.pval>crit) {
            excl.pval<-c(excl.pval,max.pval)
            exclude<-c(exclude,names(max.ind))
        }
        else { 
            go<-FALSE
        }
    }
    ecl.ma<-cbind(1:length(excl.pval),excl.pval)
    colnames(ecl.ma)<-c("Order","p-value")
    print(xtable(summary(p.lm),digits=d),type="latex",file="hw8/hw8_3_a_lm.tex")
    print(xtable(ecl.ma,digits=d),type="latex",file="hw8/hw8_3_a_excl.tex")
        
    ### b,c,d)
    subs<-regsubsets(gamble~.,data=teengamb)
    sum<-summary(subs)
    varname<-colnames(sum$which)
    print(xtable(sum$which,digits=d),type="latex",file="hw8/hw8_3_which.tex")

    aic<-nrow(teengamb)*log(sum$rss/nrow(teengamb)) + (2:ncol(teengamb))*2
    aic.min<-min(aic)
    aic.excl<-!sum$which[which.min(aic),]
    aic.lm<-lm(gamble~.,data=teengamb[!names(teengamb)%in% varname[aic.excl]])
    aic.tbl<-rbind(1:(length(varname)-1),aic)
    rownames(aic.tbl)<-c("Predictors","AIC")
    print(xtable(aic.tbl,digits=d),type="latex",file="hw8/hw8_3_b.tex")
    
    adjr2.excl<-!sum$which[which.max(sum$adjr2),]
    adjr2.max<-max(sum$adjr2)
    adjr2.lm<-lm(gamble~.,data=teengamb[!names(teengamb)%in% varname[adjr2.excl]])
    adjr2.tbl<-rbind(1:(length(varname)-1),sum$adjr2)
    rownames(adjr2.tbl)<-c("Predictors","Adjusted R Squared")
    print(xtable(adjr2.tbl,digits=d),type="latex",file="hw8/hw8_3_c.tex")

    any.cp <- sum(sum$cp<2:length(varname))
    cp.mdl <- ifelse(any.cp,which(sum$cp<2:length(varname))[1],min(sum$cp+2:length(varname)))
    cp.opt <- sum$cp[cp.mdl]
    cp.excl<-!sum$which[cp.mdl,]
    cp.lm<-lm(gamble~.,data=teengamb[!names(teengamb)%in% varname[cp.excl]])
    cp.tbl<-rbind(2:(length(varname)),sum$cp)
    rownames(cp.tbl)<-c("Parameters","Mallow's Cp")
    print(xtable(cp.tbl,digits=d),type="latex",file="hw8/hw8_3_d.tex")


}

if (prob4) {
    data(stackloss)
    n<-nrow(stackloss)
    sl.lm<-lm(stack.loss~.,data=stackloss)
    print(xtable(summary(sl.lm),digits=d),type="latex",file="hw8/hw8_4_lm.tex")

    subs<-regsubsets(stack.loss~.,data=stackloss)
    print(xtable(summary(subs)$which,digits=d),type="latex",file="hw8/hw8_4_small_which.tex")
    print(xtable(as.matrix(summary(subs)$adjr2),digits=d),type="latex",file="hw8/hw8_4_small_adjr2.tex")
    # remove Acid.Conc
    sl.s.lm <- lm(stack.loss ~ Air.Flow + Water.Temp, data=stackloss)
    print(xtable(summary(sl.s.lm),digits=d),type="latex",file="hw8/hw8_4_small_lm.tex")

    # check for outliers
    pdf('hw8/hw8_4_small_rstudent.pdf')
    qqplot(rt(1000,df=n-3),rstudent(sl.s.lm),xlab="T Distribution With 18 Degrees of Freedom", ylab="Studentized Residuals", main="QQ Plot Studentized Resisduals vs. T Distribution")
    abline(0,1)
    dev.off()
    ind.l.rstud.s<-which.max(abs(rstudent(sl.s.lm)))
    l.rstud.s<-rstudent(sl.s.lm)[ind.l.rstud.s]
    crit.s <-qt(.025/n,n-3)

    pdf('hw8/hw8_4_small_cook.pdf')
    halfnorm(cooks.distance(sl.s.lm),3,ylab="Cook's Distance",main="Halfnorm Cook's Distances")
    dev.off()

    # Do the same for big model
    pdf('hw8/hw8_4_rstudent.pdf')
    qqplot(rt(1000,df=n-3),rstudent(sl.lm),xlab="T Distribution With 18 Degrees of Freedom", ylab="Studentized Residuals", main="QQ Plot Studentized Resisduals vs. T Distribution")
    abline(0,1)
    dev.off()
    ind.l.rstud<-which.max(abs(rstudent(sl.lm)))
    l.rstud<-rstudent(sl.lm)[ind.l.rstud]
    crit <-qt(.025/(n-1),(n-1)-3)

    pdf('hw8/hw8_4_cook.pdf')
    halfnorm(cooks.distance(sl.lm),3,ylab="Cook's Distance",main="Halfnorm Cook's Distances")
    dev.off()

    subs.no<-regsubsets(stack.loss~.,data=stackloss[-21,])
    print(xtable(summary(subs.no)$which,digits=d),type="latex",file="hw8/hw8_4_no_which.tex")
    print(xtable(as.matrix(summary(subs.no)$adjr2),digits=d),type="latex",file="hw8/hw8_4_no_adjr2.tex")
}

if (prob5) {
    data(seatpos)
    pr.seat<-prcomp(seatpos[-9])
    X<-scale(seatpos[-9],scale=F)

    ### a
    sing.vals <- sort(sqrt(abs(eigen(t(X)%*%X)$values)),decreasing=TRUE)
    pdf('hw8/5_a.pdf')
    plot(sing.vals,ylab="Singular Value",xlab="Rank of Singular Value",main="Singular Values of Centered Design Matrix")
    dev.off()

    ### b
    pdf('hw8/5_b.pdf')
    matplot(pr.seat$x[,c(1,2)],pch=c("1","2"),col=c(1,1),xlab="Vector Component",ylab="Component Value",main="Values of Two Largest Singular Vectors")
    dev.off()

    ### c
    PC<-pr.seat$x[,c(1,2)]
    PC.lm <-lm(hipcenter~PC,data=seatpos)
    all.lm<-lm(hipcenter~.,data=seatpos)
    sm.lm<-lm(hipcenter~Age+Ht,data=seatpos)

    print(xtable(summary(PC.lm),digits=d),type="latex",file="hw8/hw8_5_pc_lm.tex")
    print(xtable(summary(all.lm),digits=d),type="latex",file="hw8/hw8_5_all_lm.tex")
    print(xtable(summary(sm.lm),digits=d),type="latex",file="hw8/hw8_5_sm_lm.tex")
    
}
    
