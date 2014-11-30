###### Header
library('parallel')
library('locfit')
library('gam')
library('MASS')
library('xtable')
library('reshape2')
library('chron')
library('graphics')

# Turns on different parts of the code
clean.data <- F 
data.sum   <- F 
par.reg    <- F 
nonpar.reg <- F 
eval.fit   <- T 

# par.reg options
cross.validate <- T
procs <- 2

# global function
# I'd be shocked if there isn't a time/date function to convert military time to 12 hour time with am/pm, but I'm too lazy to find it
oclock<-function(h) paste(((as.numeric(h)-1)%%12)+1,ifelse(as.numeric(h) %/%12,"pm","am"),sep="")
# Returns a sequence of colors along gray scale, since color printing is a scam
gray.scale<-function(n) gray(seq(1,0,length.out=n))

###### Data cleaning
if (clean.data) {
    ### Generate day/hour counts of robberies
    raw.robberies<-read.csv("final_project/robberies_2011_to_2013.csv",header=TRUE,as.is=TRUE)
    rob.datetime<-as.POSIXct(raw.robberies$DateTime, "CST", format="%m/%d/%Y %I:%M:%S %p")
    hourly.robberies <- aggregate(1:length(rob.datetime),by=list(hours(rob.datetime),as.Date(rob.datetime)),FUN=length)
    colnames(hourly.robberies)<-c("hour","date","robberies")
    rm(raw.robberies,rob.datetime)

    ### Clean up weather info.
    raw.weather<-read.csv("final_project/ohare_hrly_weather.csv",header=TRUE,as.is=TRUE)
    wea.datetime<-as.POSIXct(paste(raw.weather$Date, raw.weather$LocalStdTime), "CST", format="%Y-%m-%d %H:%M:%S")
    hourly.weather<- data.frame(datetime=wea.datetime, date=as.Date(wea.datetime), hour=hours(wea.datetime), temp=as.numeric(ifelse(raw.weather$AirTemp==" M ",NA,raw.weather$AirTemp)), precip=as.numeric(ifelse(raw.weather$PrecipTotal=="-",0,ifelse(raw.weather$PrecipTotal==" T",.005,raw.weather$PrecipTotal))))
    # We want to kick out observations with missing temp 
    hourly.weather<-hourly.weather[!is.na(hourly.weather$temp),]
    # We also want only one observation per hour, perferably the standard one taken at 51min after the hour. Otherwise, take first in the hour or first value of dupplicated times
    hourly.weather <- hourly.weather[order(hourly.weather$date,hourly.weather$hour,minutes(hourly.weather$datetime)==51,hourly.weather$datetime),]
    hourly.weather <- hourly.weather[!duplicated(hourly.weather[,c("date","hour")]),-1]
    rm(raw.weather,wea.datetime)
    
    ### Merge the two data sets. Weather is only missing for a few hours, so we just drop them
    rob.data <- merge(hourly.robberies,hourly.weather,by=c("date","hour"),all.y=TRUE,all.x=FALSE)
    rob.data$robberies[is.na(rob.data$robberies)] <- 0
    # Mark weekend/holidays
    rob.data$workday <- !(is.weekend(rob.data$date) | is.holiday(rob.data$date))
    # Pick a 20% subset as our test data set
    rob.data$test.set <- FALSE
    rob.data$test.set[sample.int(nrow(rob.data),size=floor(nrow(rob.data)*.2))] <- TRUE 
    rm(hourly.weather,hourly.robberies)
    save(rob.data,file="final_project/hourly_robberies.RData")

}
###### Data Summary
if (data.sum) {
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }
    rob.data$temp.buc<-cut(rob.data$temp,breaks=c(-100,30,50,70,200))
    rob.data$posprecip<-rob.data$precip>.0

    rd.w <-rob.data[ rob.data$workday,]
    rd.nw<-rob.data[!rob.data$workday,]

    # Rate of robberies at different temperature ranges
    tb.w<-tapply(rd.w$robberies,list(rd.w$hour,rd.w$temp.buc),FUN=mean)
    pdf("final_project/tempbuc_w.pdf")
    matplot(0:23,tb.w,type="l",col=rep(1,4),xlab="Hour (Military Time)",ylab="Mean Robberies Per Hour",main="Rate of Robberies at Different Temperatures, Workdays",xaxp=c(0,24,8),ylim=c(0,3),yaxp=c(0,3,6))
    legend('bottomright',c("Below 30 Degrees","30-50 Degrees","50-70 Degrees","Above 70 Degrees"),lty=1:4)
    dev.off()

    tb.nw<-tapply(rd.nw$robberies,list(rd.nw$hour,rd.nw$temp.buc),FUN=mean)
    pdf("final_project/tempbuc_nw.pdf")
    matplot(0:23,tb.nw,type="l",col=rep(1,4),xlab="Hour (Military Time)",ylab="Mean Robberies Per Hour",main="Rate of Robberies at Different Temperatures, Weekend/Holiday",xaxp=c(0,24,8),ylim=c(0,3),yaxp=c(0,3,6))
    legend('bottomright',c("Below 30 Degrees","30-50 Degrees","50-70 Degrees","Above 70 Degrees"),lty=1:4)
    dev.off()

    # Bar plot of counts at different times
    hour.seq<-seq(2,23,8)
    rob.ct.w<-table(rd.w$hour[rd.w$hour %in% hour.seq],rd.w$robberies[rd.w$hour %in% hour.seq])
    rob.ct.nw<-table(rd.nw$hour[rd.nw$hour %in% hour.seq],rd.nw$robberies[rd.nw$hour %in% hour.seq])
    pdf("final_project/ct_dist_w.pdf")
    barplot(rob.ct.w,main="Distribution of Robbery Counts, Workdays",xlab="Number of Robberies During Hour",ylab="Occurrences in Data",col=gray.scale(length(hour.seq)),legend=oclock(rownames(rob.ct.w)),beside=TRUE)
    dev.off()
    pdf("final_project/ct_dist_nw.pdf")
    barplot(rob.ct.nw,main="Distribution of Robbery Counts, Workdays",xlab="Number of Robberies During Hour",ylab="Occurrences in Data",col=gray.scale(length(hour.seq)),legend=oclock(rownames(rob.ct.nw)),beside=TRUE)
    dev.off()

    # Effect of precipitation
    precip.w<-ifelse(rob.data$posprecip,ifelse(rob.data$workday,1,2),ifelse(rob.data$workday,3,4)) 
    precip.tab<-tapply(rob.data$robberies,list(rob.data$hour,precip.w),FUN=mean)
    pdf("final_project/precip.pdf")
    matplot(0:23,precip.tab,type="l",col=rep(1,4),xlab="Hour (Military Time)",ylab="Mean Robberies Per Hour",main="Rate of Robberies With/Without Precipitation",xaxp=c(0,24,8),ylim=c(0,3),yaxp=c(0,3,6),lty=c(1,3,2,4))
    legend('bottomright',c("Precipitation, Workday","Precipitation, Weekend/Holiday","No Precipitation, Workday","No Precipitation, Weekend/Holiday"),lty=c(1,3,2,4))
    dev.off()

}

###### Parametric Model
if (par.reg) {
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }

    ### Poisson model, selected based on Akaike Criteria
    rob.poi<-glm(robberies~as.factor(hour)+as.factor(floor(hour/3))*workday+as.factor(floor(hour/3))*log(temp+50),family=poisson(link=log),data=rob.data,subset=!test.set)

    save(rob.poi,file="final_project/poisson_model.RData")
}

###### Non-Parametric Model
if (nonpar.reg) {
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }
    
    # Data to run regression on
    rd<-rob.data[!rob.data$test.set,]
    
    # Transform the count data
    rd$logrob<-log(rd$robberies+1)

    ### Don't want a discontinuity at 12am, so add back in data from prior and next day to get from -24 to 48 hours
    rd.m <- rd
    rd.m$date<-rd.m$date+1
    rd.m$hour<-rd.m$hour-24
    rd.p <- rd
    rd.p$date<-rd.p$date-1
    rd.p$hour<-rd.p$hour+24
    rd<-rbind(rd.m,rd,rd.p)
    rd$workday <- !(is.weekend(rd$date) | is.holiday(rd$date))
    rm(rd.m,rd.p)

    # Want to select bandwith and scale for second variable (effectively two bandwiths). Test a grid of values and pick the best one. Don't care about hours outside 0-23. 
    param.ll<-list() # This will be the list of parameters for fitting the model
    if (cross.validate) {
        sub<-rd$hour>=0 & rd$hour<24
        cv.score<-function(nn,s2,data,subset) {
            cv.fit<-fitted(locfit(logrob~lp(temp,hour,deg=1,nn=nn,scale=c(1,s2)),kern="gauss",data=data),what="coef",cv=TRUE)[subset]
            return(sum((data$logrob[subset]-cv.fit)^2)/sum(subset))
        }
        nnvals<-seq(.01,.2,.01)
        svals<-seq(.01,.2,.01)
        cv.mat.w<-NULL
        cv.mat.nw<-NULL
        for (s in svals) {
            cv.mat.w<-cbind(cv.mat.w,mcmapply(cv.score,nn=nnvals,MoreArgs=list(s2=s,data=rd,subset=sub & rd$workday),mc.cores=procs))
            cv.mat.nw<-cbind(cv.mat.nw,mcmapply(cv.score,nn=nnvals,MoreArgs=list(s2=s,data=rd,subset=sub & !rd$workday),mc.cores=procs))
        }
        # fitting parameters of the model
        ind.w<-which(cv.mat.w==min(cv.mat.w),arr.ind=TRUE)
        param.ll$nnval.w<-nnvals[ind.w[1]]
        param.ll$s2.w<-svals[ind.w[2]]
        ind.nw<-which(cv.mat.nw==min(cv.mat.nw),arr.ind=TRUE)
        param.ll$nnval.nw<-nnvals[ind.nw[1]]
        param.ll$s2.nw<-svals[ind.nw[2]]
        save(param.ll,file="final_project/ll_param.RData")
    }
    else {
        # These values were found from previous cross validation
        load(file="final_project/ll_param.RData")
    }


    # Fit local linear model 
    rob.ll.w<-locfit(logrob~lp(temp,hour,deg=1,nn=param.ll$nnval.w,scale=c(1,param.ll$s2.w)),kern="gauss",data=rd,subset=workday)
    rob.ll.nw<-locfit(logrob~lp(temp,hour,deg=1,nn=param.ll$nnval.nw,scale=c(1,param.ll$s2.nw)),kern="gauss",data=rd,subset=!workday)

    save(rob.ll.w,rob.ll.nw,file="final_project/local_linear_model.RData")

}
if (eval.fit) {
    # Load data 
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }
    # Load models
    if (!exists("rob.ll.w")) {
        load(file="final_project/local_linear_model.RData")
    }
    if (!exists("rob.poi")) {
        load(file="final_project/poisson_model.RData")
    }

    ### Generate contours of fit. 

    # First generate data set to predict on
    ntemp<-250
    nhour<-250
    rtemp<-c(-10,100)
    rhour<-c(0,23.999)
    temp.seq<-seq(rtemp[1],rtemp[2],length.out=ntemp)
    hour.seq<-seq(rhour[1],rhour[2],length.out=nhour)
    contour.data <-data.frame(hour=rep(hour.seq,ntemp),temp=(as.numeric(gl(ntemp,nhour))-1)*(rtemp[2]-rtemp[1])/(ntemp-1)+rtemp[1])

    # Make matrix of predicted values for weekend/weekdays for each model
    # Local linear
    predmat.ll.w<-matrix(exp(predict(rob.ll.w,newdata=contour.data))-1,nrow=ntemp)
    predmat.ll.nw<-matrix(exp(predict(rob.ll.nw,newdata=contour.data))-1,nrow=ntemp)
    # poisson
    contour.data$hour<-floor(contour.data$hour)
    contour.data$workday<-TRUE
    predmat.poi.w<-matrix(predict(rob.poi,newdata=contour.data,type="response"),nrow=ntemp)
    contour.data$workday<-FALSE
    predmat.poi.nw<-matrix(predict(rob.poi,newdata=contour.data,type="response"),nrow=ntemp)

    # make the contour plots
    pred.mats<-list(predmat.poi.w,predmat.ll.w,predmat.poi.nw,predmat.ll.nw)
    names(pred.mats)<-c("poi_w","ll_w","poi_nw","ll_nw")
    pred.mats$titles<-c("Poisson Model, Workdays","Local Linear Model, Workdays","Poisson Model, Weekend/Holiday","Local Linear Model, Weekend/Holiday")
    # parameters for plots
    max.pred<-max(pred.mats[[1]],pred.mats[[2]],pred.mats[[3]],pred.mats[[4]])
    lvls<-pretty(c(0,max.pred),20)
    for (i in 1:4) {
        for (fill in c(T,F)) {
            pdf(paste("final_project/",names(pred.mats)[i],"_contour",ifelse(fill,"_fill",""),".pdf",sep=""))
            if (fill) filled.contour(x=hour.seq,y=temp.seq,z=pred.mats[[i]],levels=lvls,color.palette=gray.scale)
            else    contour(x=hour.seq,y=temp.seq,z=pred.mats[[i]],levels=lvls)
            title(main=paste("Predicted Robberies,", pred.mats$titles[i]),xlab="Hour (Military Time)", ylab="Temperature")
            par(xaxp=c(round(rhour),8),yaxp=c(round(rtemp),11))
            dev.off()
        }
    }

    ### Generate confidence bands arounds lines showing robberies at different temperatures for a few set hours
    # Data to predict on - just use a few sets of hours
    hour.seq<-seq(2,23,8)
    nhour<-length(hour.seq)
    cb.data <-data.frame(hour=rep(hour.seq,ntemp),temp=(as.numeric(gl(ntemp,nhour))-1)*(max(temp.seq)-min(temp.seq))/(ntemp-1)+rtemp[1])

    # Matrices of predicted values and bands
    td<-function(d) matrix(d,ncol=nhour,byrow=TRUE)
    # Local linear
    ll.w.pred<-predict(rob.ll.w,newdata=cb.data,se.fit=TRUE,band="local")
    predmat.ll.w<-exp(cbind(td(ll.w.pred$fit),td(ll.w.pred$fit-crit(rob.ll.w)$crit.val*ll.w.pred$se.fit),td(ll.w.pred$fit+crit(rob.ll.w)$crit.val*ll.w.pred$se.fit))) -1
    ll.nw.pred<-predict(rob.ll.nw,newdata=cb.data,se.fit=TRUE,band="local")
    predmat.ll.nw<-exp(cbind(td(ll.nw.pred$fit),td(ll.nw.pred$fit-crit(rob.ll.nw)$crit.val*ll.nw.pred$se.fit),td(ll.nw.pred$fit+crit(rob.ll.nw)$crit.val*ll.nw.pred$se.fit))) -1
    # poisson
    cb.data$hour<-floor(cb.data$hour)
    cb.data$workday<-TRUE
    poi.w.pred<-predict(rob.poi,newdata=cb.data,type="link",se.fit=TRUE)
    predmat.poi.w<-exp(cbind(td(poi.w.pred$fit),td(poi.w.pred$fit-1.96*poi.w.pred$se.fit),td(poi.w.pred$fit+1.96*poi.w.pred$se.fit)))    
    cb.data$workday<-FALSE
    poi.nw.pred<-predict(rob.poi,newdata=cb.data,type="link",se.fit=TRUE)
    predmat.poi.nw<-exp(cbind(td(poi.nw.pred$fit),td(poi.nw.pred$fit-1.96*poi.nw.pred$se.fit),td(poi.nw.pred$fit+1.96*poi.nw.pred$se.fit)))    

    # make the plots
    pred.mats<-list(predmat.poi.w,predmat.ll.w,predmat.poi.nw,predmat.ll.nw)
    names(pred.mats)<-c("poi_w","ll_w","poi_nw","ll_nw")
    pred.mats$titles<-c("Poisson Model, Workdays","Local Linear Model, Workdays","Poisson Model, Weekend/Holiday","Local Linear Model, Weekend/Holiday")
    # parameters for plots
    max.pred<-max(pred.mats[[1]][,1:3],pred.mats[[2]][,1:3],pred.mats[[3]][,1:3],pred.mats[[4]][,1:3])
    min.pred<-min(pred.mats[[1]][,1:3],pred.mats[[2]][,1:3],pred.mats[[3]][,1:3],pred.mats[[4]][,1:3])
    ytic<-pretty(c(max.pred,min.pred))
    lwd.ci<-.4
    for (i in 1:4) {
        pdf(paste("final_project/",names(pred.mats)[i],"_ci.pdf",sep=""))
        matplot(x=temp.seq,y=pred.mats[[i]],type="l",lty=c(rep(1:nhour,3)),lwd=c(rep(1,nhour),rep(lwd.ci,2*nhour)),col=rep(1,nhour*3),xlab="Temperature",ylab="Mean Robberies Per Hour",main=paste("Predicted Rate of Robberies,",pred.mats$titles[i]),ylim=c(min(ytic),max(ytic)))
        legend('bottomright',c(oclock(hour.seq),rep("95% CI",nhour)),lty=rep(1:nhour,2),lwd=c(rep(1,nhour),rep(lwd.ci,nhour)),ncol=2)
        dev.off()
    }

    ### Table of 0 vs 60 deg fits and CIs
    temp.seq<-c(0,60)
    hour.seq<-0:23
    td <-data.frame(hour=rep(hour.seq,length(temp.seq)),temp=c(rep(temp.seq[1],length(hour.seq)),rep(temp.seq[2],length(hour.seq))))
    td=rbind(td,td)
    orig.len <-length(temp.seq)*length(hour.seq)
    td$workday<-c(rep(T,orig.len),rep(F,orig.len))
    poi.pred<-predict(rob.poi,newdata=td,type="link",se.fit=TRUE)
    ll.w.pred<-predict(rob.ll.w,newdata=td,se.fit=TRUE,band="local")
    ll.nw.pred<-predict(rob.ll.nw,newdata=td,se.fit=TRUE,band="local")
    td$poi.lk.fit <- poi.pred$fit
    td$poi.se.fit <- poi.pred$se.fit
    td$poi.crit   <- 1.96
    td$ll.lk.fit  <- ifelse(td$workday,ll.w.pred$fit,ll.nw.pred$fit)
    td$ll.se.fit  <- ifelse(td$workday,ll.w.pred$se.fit,ll.nw.pred$se.fit)
    td$ll.crit    <- ifelse(td$workday,crit(rob.ll.w)$crit.val,crit(rob.ll.nw)$crit.val)
    td[,paste("poi.",c("fit","lb","ub"),sep="")]<-exp(cbind(td$poi.lk.fit,td$poi.lk.fit-td$poi.crit*td$ll.se.fit,td$poi.lk.fit+td$poi.crit*td$ll.se.fit))
    td[,paste("ll.",c("fit","lb","ub"),sep="")]<-exp(cbind(td$ll.lk.fit,td$ll.lk.fit-td$ll.crit*td$ll.se.fit,td$ll.lk.fit+td$ll.crit*td$ll.se.fit)) - 1

    mats<-list(NULL,NULL,NULL,NULL)
    mats$names<-c("poi_w","poi_nw","ll_w","ll_nw")
    mats$wk<-c(TRUE,FALSE,TRUE,FALSE)
    mats$pf<-c("poi","poi","ll","ll")
    prefixes<-c("Lower Bound,","Estimate,","Upper Bound")
    endings<-paste(temp.seq,"Degrees")
    for (i in 1:4) {
        mats[[i]]<-t(as.matrix(td[td$workday==mats$wk[i],paste(mats$pf[i],c(".lb",".fit",".ub"),sep="")]))
        mats[[i]]<-rbind(mats[[i]][,1:length(hour.seq)],mats[[i]][,(length(hour.seq)+1):(2*length(hour.seq))])
        mats[[i]]<-rbind(mats[[i]],mats[[i]][5,]-mats[[i]][2,])
        rownames(mats[[i]])<-c(outer(prefixes,endings,FUN=paste),"Estimated Difference")
        colnames(mats[[i]])<-oclock(0:23)
        print(xtable(mats[[i]],digits=3),type="latex",file=paste("final_project/",mats$names[i],"_0_60_fit.tex",sep=""))
    }

    # Root mean squared error of estimates
    cell.mean<-lm(robberies~workday*hour*temp,data=rob.data,subset=!rob.data$test)
    test<-rob.data[rob.data$test,]
    test$logrob<-log(test$robberies+1)
    test$ll.lpred<-ifelse(test$workday,predict(rob.ll.w,newdata=test),predict(rob.ll.nw,newdata=test))
    test$poi.lpred<-log(predict(rob.poi,newdata=test,type="response")+1)
    test$cell.mean.lpred<-log(predict(cell.mean,newdata=test)+1)
    poi.rmsle<-sqrt(sum((test$poi.lpred-test$logrob)^2)/nrow(test))
    ll.rmsle<-sqrt(sum((test$ll.lpred-test$logrob)^2)/nrow(test))
    cell.mean.rmsle<-sqrt(sum((test$cell.mean.lpred-test$logrob)^2)/nrow(test))
    rmsle<-c(poi.rmsle,ll.rmsle,cell.mean.rmsle)
    names(rmsle)<-c("Poisson RMSLE","Local Linear RMSLE","Cell Mean")
    write.table(rmsle,file='final_project/rmsle.txt')
        
}
