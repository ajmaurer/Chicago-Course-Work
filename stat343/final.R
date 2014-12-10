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
library('mgcv')
options(width=200)

prob1<-T
prob2<-F
prob3<-F

res.qq<-function(model,file,title) {
    if (length(residuals(model))<5000) {
        p.val<-round(shapiro.test(residuals(model))$p.value,5)
        sub.text<-paste("Shapiro Test P-Value:",round(p.val,4))
    }
    else sub.text<-NULL   
    pdf(file)
    qqnorm(residuals(model),ylab="Residual Quantiles",main=title,sub=sub.text)
    qqline(residuals(model))    
    dev.off()
}

stud.res.qq<-function(model,file,title) {
    rstud<-rstudent(model)
    max<-abs(rstud[which.max(abs(rstud))])
    n<-length(rstud)
    p<-length(coef(model))
    bf.crit<-abs(qt(.05/(n*2),n-p))
    pdf(file)
    qqnorm(rstudent(model),ylab="Residual Quantiles",main=title,sub=paste("Bonferroni Critical Value:",round(bf.crit,4),"Most Extreme Studentized Residual:",round(max,4)))
    qqline(rstudent(model))    
    dev.off()
}


res.plot<-function(model,file,title) {
    pdf(file)
    n<-length(predict(model))
    size<-max(.05,min(.5,5/sqrt(n)))
    plot(residuals(model)~predict(model),ylab="Residuals",xlab="Predicted",main=title,cex=size)
    plot(locfit(residuals(model)~lp(predict(model),nn=.2,deg=1)),add=T)
    legend('topleft',"Local Linear Fit",lty=1)
    dev.off()
}

var.res.plot<-function(model,data,stub,title,terms) {
    n<-length(predict(model))
    size<-max(.05,min(.5,5/sqrt(n)))
    for (i in terms) {
        name<-names(data)[i]
        pdf(paste(stub,name,'.pdf',sep=""))
        plot(residuals(model)~data[,i],ylab="Residuals",xlab=name,main=paste(title,", ",name,sep=""),cex=size)
        plot(locfit(residuals(model)~lp(data[,i],nn=.3,deg=1)),add=T)
        legend('topleft',"Local Linear Fit",lty=1)
        dev.off()
    }
}
        
sc.res.plot<-function(model,file,title) {
    pdf(file)
    n<-length(predict(model))
    size<-max(.05,min(.5,5/sqrt(n)))
    plot(sqrt(abs(residuals(model)))~predict(model),ylab=expression(sqrt("|Residuals|")),xlab="Predicted",main=title,cex=size)
    plot(locfit(sqrt(abs(residuals(model)))~lp(predict(model),nn=.25,deg=1)),add=T)
    legend('topleft',"Local Linear Fit",lty=1)
    dev.off()
}

if (prob1) {
    sink("final/1_out.txt",split=T,type="output")
    af<-read.delim('final/airfoil_self_noise.dat',col.names=c("freq","ang","cl","vel","thi","pre"))
    mins<-apply(af,2,min)
    maxs<-apply(af,2,max)
    af.s<-data.frame(scale(af,center=mins,scale=maxs-mins))

    ### part a 
    print("PART A")
    af.lm<-lm(pre~.,data=af)
    print(summary(af.lm))

    # Analyze various problems
    # Residual Normality
    res.qq(af.lm,'final/1a_qq.pdf','Basic Linear QQ')

    # Constant variance
    sc.res.plot(af.lm,'final/1a_sqrt_res_plot.pdf','Basic Linear Scaled Residuals')
    print(var.test(residuals(af.lm)[af$ang==0],residuals(af.lm)[af$ang>0],data=af))

    # Solution
    # first, normalize all predictors to 0-1
    mins<-apply(af,2,min)
    maxs<-apply(af,2,max)
    af.s<-data.frame(scale(af,center=mins,scale=maxs-mins))

    # Then, square root transformations all around
    af.t.lm<-lm(pre~freq+sqrt(1-ang)+sqrt(1-cl)+sqrt(vel)+sqrt(1-thi),data=af.s)
    print("After transforming once")
    print(summary(af.t.lm))
    res.qq(af.t.lm,'final/1a_qq_lm2.pdf','Transformed Linear Model Res QQ')
    sc.res.plot(af.t.lm,'final/1a_sqrt_res_plot_lm2.pdf','Transformed Linear Model Scaled Residuals')
    print(var.test(residuals(af.t.lm)[af$ang==0],residuals(af.t.lm)[af$ang>0],data=af))
                
    # This transformation is ridiculous, not using it
    # This gets the error approximately normally distributed
    #trans1<-function(y,s,c,p) (y-c)/exp(abs(s*(y-c))^p)
    #af.t.lm<-lm(trans1(pre,.065,125,50)~.,data=af)
    #res.qq(af.t.lm,'final/1a_qq_tran.pdf','Extreme Transformed Linear Model QQ')
    #sc.res.plot(af.t.lm,'final/1a_sqrt_res_plot_tran.pdf','Extreme Transformed Linear Model Residuals')
    #print(summary(af.t.lm))
    #print(var.test(residuals(af.t.lm)[af$ang==0],residuals(af.t.lm)[af$ang>0],data=af))

    ### b
    # look for nonlinearity in last model
    res.plot(af.t.lm,'final/1b_res_plot_lm2.pdf','Transformed Linear Model Residuals')
    var.res.plot(af.t.lm,data=af.s,stub='final/1b_res_plot_tran_',title='Transformed Linear Model Residuals',terms=c(1,2,5))

    # Try something more extreme 
    af.x.lm<-lm(sqrt(pre)~polym(freq,ang,cl,degree=4)+sqrt(vel)+sqrt(1-thi),data=af.s)
    print(summary(af.x.lm))
    res.plot(af.x.lm,'final/1b_res_plot_int.pdf','Interacted Linear Model Residuals')
    var.res.plot(af.x.lm,data=af.s,stub='final/1b_res_plot_int_',title='Interacted Linear Model Residuals',terms=c(1,2,5))
    sc.res.plot(af.x.lm,'final/1b_sqrt_res_plot_int.pdf','Interacted Linear Model Scaled Residuals')

    # everything as factor 
    #af.alt.lm<-lm(pre~as.factor(freq)+as.factor(ang)+as.factor(cl)+as.factor(vel)+as.factor(thi),data=af)
    #print(summary(af.alt.lm))
    #res.plot(af.alt.lm,'final/1b_res_plot_int.pdf','All Factor Model Residuals')
    #res.qq(af.alt.lm,'final/1b_qq_int.pdf','All Factor Model QQ')

    # Best fit... though completely ridiculous
    #af.wtf.lm<-lm(pre^-5~poly(freq,cl,degree=2)*as.factor(ang)*vel*thi,data=af)
    #print(summary(af.wtf.lm))
    #res.plot(af.wtf.lm,'final/1b_res_plot_wtf.pdf','Heavily Interacted Model Residuals')
    #res.qq(af.wtf.lm,'final/1b_qq_wtf.pdf','Heavily Interacted Model QQ')

    ### part d
    af.tree<-rpart(pre~freq+ang+cl+vel+thi,data=af,method='anova',control=rpart.control(cp=0,minsplit=2,minbucket=1,maxdepth=30,xval=10))
    af.pr.tree<-prune(af.tree,cp=af.tree$cptable[which.min(af.tree$cptable[,"xerror"]),"CP"])
    res.plot(af.pr.tree,'final/1d_res_plot_tree.pdf','Regression Tree Residuals')
    res.qq(af.pr.tree,'final/1d_qq_tree.pdf','Regression Tree QQ')
    TSS<-sum((af$pre-mean(af$pre))^2)
    RSS<-sum(residuals(af.pr.tree)^2)
    n<-nrow(af)
    p<-1+max(af.pr.tree$cptable[,'nsplit'])
    print(paste("Complexity parameter is",af.tree$cptable[which.min(af.tree$cptable[,"xerror"])],"chosen at an estimated cross validation error of ",min(af.tree$cptable[,"xerror"])))
    print(paste("n is",n))
    print(paste("p is",p))
    print(paste("R2 for Cart:",round(1-RSS/TSS,4)))
    print(paste("Adjusted R2 for Cart:",round(1-(RSS/(n-p))/(TSS/(n-1)),4)))

    sink()
}

if (prob2) {
    sink("final/2_out.txt",split=T,type="output")
    cch<-read.csv(file='final/CCH.csv',header=T)

    ### a)
    cch.lm<-lm(PE~.,data=cch)
    print(summary(cch.lm))

    # i. normality
    res.qq(cch.lm,'final/2ai_qq.pdf','Residuals QQ')

    # ii. outliers
    pdf('final/2aii_cook.pdf')
    halfnorm(cooks.distance(cch.lm),5,ylab="Cook's Distance",main="Cook's Distances Half Normal")
    dev.off()
    stud.res.qq(cch.lm,'final/2aii_stud_qq.pdf','Studetnized Residuals QQ')

    # iii. structure model
    res.plot(cch.lm,'final/2aiii_res_plot.pdf','Residuals')
    names<-names(cch.lm$coefficients[-1])
    for (i in 1:4) {
        pdf(paste('final/2aiii_pr_',names[i],'.pdf',sep=""))
        termplot(cch.lm,partial.resid=T,terms=i,col.res=1,cex=.2,col.smth=1)
        dev.off()
    }

    ### b) bootstrap to determine whether q most extreme points are outliers
    q <- 10
    b <- 2000
    sigma <- summary(cch.lm)$sigma
    Yhat <- predict(cch.lm)
    X <- as.matrix(cch[1:4])

    # i. using cook distance
    cook.boot  <- function(sigma,Yhat,X,q) {
        cooks <- sort(cooks.distance(lm(rnorm(length(Yhat),Yhat,sigma)~X)),decreasing=T)
        return(cooks[1:q])
    }
    cook.mat <- t(replicate(b,cook.boot(sigma,Yhat,X,q)))
    sort.cook.mat <- apply(cook.mat,2,sort)
    x <- 1-(1:b)/b
    ymax <- sort.cook.mat[ceiling(.995*b),1]
    pdf('final/2bi_plot.pdf')
    matplot(y=sort.cook.mat,x=x,lty=rep(2:6,2),col=rep(1,10),type="l",xlim=c(0,.1),ylim=c(0,ymax),xlab="P-Value",ylab="Critical Value",main="Critical Value, q-th Largest Cook's Distance")
    xseq <- seq(.01,.1,.01)
    yseq <- rep(NULL,10)
    for (i in 1:length(xseq)) {
        yseq[i] <- sort.cook.mat[ceiling((1-xseq[i])*b),i]
    }
    text(x=xseq,y=yseq,pos=rep(4,10),labels=c("1st","2nd","3rd",paste(4:10,"th",sep="")),offset=-.6)
    dev.off()

    cook.crits <- sort.cook.mat[ceiling(b*(1-c(.01,.05,.1))),]
    rownames(cook.crits) <- c(".01",".05",".1")
    colnames(cook.crits) <- 1:10
    cook.lar <- sort(cooks.distance(cch.lm),decreasing=T)[1:10]
    names(cook.lar) <- 1:10
    
    print("10 largest cook values")
    print(round(cook.lar,4))
    print("Critical values")
    print(round(cook.crits,4))
    print("Cook Value - Do Extreme Points Fall Above Critical Values?")
    print(t(cook.crits)<cook.lar)


    # ii. using studetnized residual
    rstud.boot <- function(sigma,Yhat,X,q) {
        rstuds <- sort(abs(rstudent(lm(rnorm(length(Yhat),Yhat,sigma)~X))),decreasing=T)
        return(rstuds[1:q])
    }
    rstud.mat <- t(replicate(b,rstud.boot(sigma,Yhat,X,q)))
    sort.rstud.mat <- apply(rstud.mat,2,sort)
    x <- 1-(1:b)/b
    ymax <- sort.rstud.mat[ceiling(.995*b),1]
    pdf('final/2bii_plot.pdf')
    matplot(y=sort.rstud.mat,x=x,lty=rep(2:6,2),col=rep(1,10),type="l",xlim=c(0,.1),ylim=c(0,ymax),xlab="P-Value",ylab="Critical Value",main="Critical Value, q-th Largest Absolute Studentized Residual")
    xseq <- seq(.01,.1,.01)
    yseq <- rep(NULL,10)
    for (i in 1:length(xseq)) {
        yseq[i] <- sort.rstud.mat[ceiling((1-xseq[i])*b),i]
    }
    text(x=xseq,y=yseq,pos=rep(4,10),labels=c("1st","2nd","3rd",paste(4:10,"th",sep="")),offset=-.6)
    dev.off()

    rstud.crits <- sort.rstud.mat[ceiling(b*(1-c(.01,.05,.1))),]
    rownames(rstud.crits) <- c(".01",".05",".1")
    colnames(rstud.crits) <- 1:10
    lar.rstud <- sort(abs(rstudent(cch.lm)),decreasing=T)[1:10]
    names(lar.rstud) <- 1:10

    print("10 largest absolute studetnized residuals")
    print(round(lar.rstud,4))
    print("Critical values")
    print(round(rstud.crits,4))
    print("Studentized Residual - Do Extreme Points Fall Above Critical Values?")
    print(t(rstud.crits)<lar.rstud)

    sink()
}

if (prob3) {
    sink('final/3_out.txt')
    library('chemometrics')
    library('pls')
    library('lars')
    data('NIR')

    set.seed(1968)
    nir<-cbind(NIR$yGlcEtOH$Glucose,NIR$xNIR)   
    names(nir)[1]<-"gluc"
    nir<-nir[sample(1:nrow(nir),nrow(nir)),]
    nir.tr<-nir[1:126,]
    nir.te<-nir[127:166,]
    np<-c(nrow(nir.tr[-1]),ncol(nir.tr[-1]))

    k<-10

    rmse<-function(x,y) sqrt(mean((x-y)^2))

    ### b)
    print("PART B")
    # check for outliers ... or not so much 
    #nir.rob<-cov.rob(nir.tr[-1])
    #md<-mahalanobis(nir,center=nir.rob$center,cov=nir.rob$cov)
    #plot(qchisq(1:n/(n+1),p),sort(md))

    nir.pcr.k<-pcr(gluc~.,data=nir.tr,ncomp=100,validation="CV",segments=k)
    nir.pcr.k.cv<-RMSEP(nir.pcr.k,estimate="CV")
    nir.pcr.k.cv.comp<-which.min(nir.pcr.k.cv$val)
    pdf('final/3b_pcr_cv.pdf')
    plot(nir.pcr.k.cv,main=paste(k,"Fold Cross Validation, PCR"))
    dev.off()
    print(paste(k,"fold CV for PCR minimum at",nir.pcr.k.cv.comp,"Components, estimated RMSE of",nir.pcr.k.cv$val[nir.pcr.k.cv.comp]))
    print(paste("Test Set RMSE of",rmse(nir.te[1],predict(nir.pcr.k,newdata=nir.te,ncomp=nir.pcr.k.cv.comp))))

    # Scree plot
    nir.pca<-prcomp(nir.tr[-1])
    pdf('final/3b_pcr_scree.pdf')
    plot(nir.pca$sdev[1:100],type="l",ylab="Standard Deviation of Primary Component",xlab="Primary Component Number",main="PCR Scree Plot")
    dev.off()
    nir.scree.comp<-10
    print(paste("Select ",nir.scree.comp,"Components from Scree for PCR"))
    print(paste("Test Set RMSE of",rmse(nir.te[1],predict(nir.pcr.k,newdata=nir.te,ncomp=nir.scree.comp))))

    ### c)
    print("PART C")
    nir.pls.k<-plsr(gluc~.,data=nir.tr,ncomp=100,validation="CV",segments=k)
    nir.pls.k.cv<-RMSEP(nir.pls.k,estimate="CV")
    nir.pls.k.cv.comp<-which.min(nir.pls.k.cv$val)
    pdf('final/3c_pls_cv.pdf')
    plot(nir.pls.k.cv,main=paste(k,"Fold Cross Validation, PLS"))
    dev.off()
    
    print(paste(k,"fold CV for PLS minimum at",nir.pls.k.cv.comp,"Components, estimated RMSE of",nir.pls.k.cv$val[nir.pls.k.cv.comp]))
    print(paste("Test Set RMSE of",rmse(nir.te[1],predict(nir.pls.k,newdata=nir.te,ncomp=nir.pls.k.cv.comp))))

    ### d)
    print("PART D")

    nir.pcr.loo<-pcr(gluc~.,data=nir.tr,ncomp=100,validation="LOO")
    nir.pcr.loo.cv<-RMSEP(nir.pcr.loo,estimate="CV")
    nir.pcr.loo.comp<-which.min(nir.pcr.loo.cv$val)
    pdf('final/3d_pcr_loo.pdf')
    plot(nir.pcr.loo.cv,main=paste("Leave One Out Cross Validation, PCR"))
    dev.off()
    print(paste("PCR LOO CV minimum at",nir.pcr.loo.comp,"Components, estimated RMSE of",nir.pcr.loo.cv$val[nir.pcr.loo.comp]))
    print(paste("Test Set RMSE of",rmse(nir.te[1],predict(nir.pcr.loo,newdata=nir.te,ncomp=nir.pcr.loo.comp))))

    nir.pls.loo<-plsr(gluc~.,data=nir.tr,ncomp=100,validation="LOO")
    nir.pls.loo.cv<-RMSEP(nir.pls.loo,estimate="CV")
    nir.pls.loo.comp<-which.min(nir.pls.loo.cv$val)
    pdf('final/3d_pls_loo.pdf')
    plot(nir.pls.loo.cv,main=paste("Leave One Out Cross Validation, pls"))
    dev.off()
    print(paste("pls LOO CV minimum at",nir.pls.loo.comp,"Components, estimated RMSE of",nir.pls.loo.cv$val[nir.pls.loo.comp]))
    print(paste("Test Set RMSE of",rmse(nir.te[1],predict(nir.pls.loo,newdata=nir.te,ncomp=nir.pls.loo.comp))))

    ### e)
    print("PART E")

    nir.lars<-lars(as.matrix(nir.tr[,-1]),nir.tr$gluc)
    nir.lars.cv<-cv.lars(as.matrix(nir.tr[,-1]),nir.tr$gluc,K=k,index=seq(0,.2,length.out=200))
    nir.lars.t<-nir.lars.cv$index[which.min(nir.lars.cv$cv)]
    nir.lars.mincv<-sqrt(min(nir.lars.cv$cv))
    pdf('final/3e_cv_plot.pdf')
    plot(nir.lars.cv$index,sqrt(nir.lars.cv$cv),type="l",ylab="RSMEP",xlab="|beta|/max|beta|",main=paste(k,"Fold Cross Validation, LASSO"))
    dev.off()
    pdf('final/3e_lasso_plot.pdf')
    plot(nir.lars,xlim=c(0,.001),ylim=c(-75,75))
    dev.off()
    print(paste(k,"fold CV for LASSO minimum at",nir.lars.t,"Components, estimated RMSE of",nir.lars.mincv))
    print(paste(sum(predict(nir.lars,s=nir.lars.t,type="coef",mode="fraction")$coefficients>0),"positive coefficients"))
    print(paste("Test Set RMSE of",rmse(nir.te[1],predict(nir.lars,newx=as.matrix(nir.te[,-1]),s=nir.lars.t,mode="fraction")$fit)))















    sink()
}
