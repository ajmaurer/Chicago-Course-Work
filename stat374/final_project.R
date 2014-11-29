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
data.sum   <- T
par.reg    <- F
nonpar.reg <- F 
eval.fit   <- T

# par.reg options
cross.validate <- F

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
    rob.poi<-glm(glm(robberies~as.factor(hour)+as.factor(floor(hour/3))*workday+as.factor(floor(hour/3))*log(temp+50),family=poisson(link=log),data=rob.data,subset=test.set))

    save(rob.poi,file="final_project/poisson_model.RData")
}

###### Non-Parametric Model
if (nonpar.reg) {
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }
    
    # Data to run regression on
    rd<-rob.data[rob.data$test.set,]
    
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
        nnvals<-seq(.02,.1,.01)
        svals<-seq(.06,.2,.01)
        cv.mat.w<-NULL
        cv.mat.nw<-NULL
        for (s in svals) {
            cv.mat.w<-cbind(cv.mat.w,mcmapply(cv.score,nn=nnvals,MoreArgs=list(s2=s,data=rd,subset=sub & rd$workday),mc.cores=4))
            cv.mat.nw<-cbind(cv.mat.nw,mcmapply(cv.score,nn=nnvals,MoreArgs=list(s2=s,data=rd,subset=sub & !rd$workday),mc.cores=4))
        }
        # fitting parameters of the model
        ind.w<-which(cv.mat.w==min(cv.mat.w),arr.ind=TRUE)
        param.ll$nnval.w<-nnvals[ind.w[1]]
        param.ll$s2.w<-svals[ind.w[2]]
        ind.nw<-which(cv.mat.nw==min(cv.mat.nw),arr.ind=TRUE)
        param.ll$nnval.nw<-nnvals[ind.nw[1]]
        param.ll$s2.nw<-svals[ind.nw[2]]
    }
    else {
        # These values were found from previous cross validation
        param.ll$nnval.w<-.03
        param.ll$s2.w<-.09
        param.ll$nnval.nw<-.03
        param.ll$s2.nw<-.14
    }

    # Fit local linear model 
    rob.ll.w<-locfit(logrob~lp(temp,hour,deg=1,nn=param.ll$nnval.w,scale=c(1,param.ll$s2.w)),kern="gauss",data=rd,subset=workday)
    rob.ll.nw<-locfit(logrob~lp(temp,hour,deg=1,nn=param.ll$nnval.nw,scale=c(1,param.ll$s2.nw)),kern="gauss",data=rd,subset=!workday)

    save(param.ll,rob.ll.w,rob.ll.nw,file="final_project/local_linear_model.RData")

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
    predmat.poi.w<-matrix(predict(rob.poi,newdata=contour.data),nrow=ntemp)
    contour.data$workday<-FALSE
    predmat.poi.nw<-matrix(predict(rob.poi,newdata=contour.data),nrow=ntemp)

    # make the contour plots
    pred.mats<-list(predmat.poi.w,predmat.ll.w,predmat.poi.nw,predmat.ll.nw)
    names(pred.mats)<-c("poi_w","ll_w","poi_nw","ll_nw")
    pred.mats$titles<-c("Poisson Model, Workdays","Local Linear Model, Workdays","Poisson Model, Weekend/Holiday","Local Linear Model, Weekend/Holiday")
    # parameters for plots
    max.pred<-max(pred.mats[[1]],pred.mats[[2]],pred.mats[[3]],pred.mats[[4]])
    lvls<-pretty(c(0,max.pred),20)
    gray.scale<-function(n) gray(seq(1,0,length.out=n))
    for (i in 1:4) {
        pdf(paste("final_project/",names(pred.mats)[i],"_contour.pdf",sep=""))
        #filled.contour(x=hour.seq,y=temp.seq,z=pred.mats[[i]],levels=lvls,color.palette=gray.scale)
        contour(x=hour.seq,y=temp.seq,z=pred.mats[[i]],levels=lvls)
        title(main=paste("Predicted Robberies,", pred.mats$titles[i]),xlab="Hour (Military Time)", ylab="Temperature")
        par(xaxp=c(round(rhour),8),yaxp=c(round(rtemp),11))
        dev.off()
    }

    ### Generate confidence bands arounds lines showing robberies at different temperatures for a few set hours
    # Data to predict on - just use a few sets of hours
    hour.seq<-seq(2,23,8)
    nhour<-length(hour.seq)
    cb.data <-data.frame(hour=rep(hour.seq,ntemp),temp=(as.numeric(gl(ntemp,nhour))-1)*(max(temp.seq)-min(temp.seq))/(ntemp-1)+rtemp[1])

    # Matrices of predicted values and bands
    # Local linear
    ll.w.pred<-predict(rob.ll.w,newdata=cb.data,se.fit=TRUE,band="local")
    predmat.ll.w<-exp(matrix(cbind(ll.w.pred$fit,ll.w.pred$fit-crit(rob.ll.w)$crit.val*ll.w.pred$se.fit,ll.w.pred$fit+crit(rob.ll.w)$crit.val*ll.w.pred$se.fit),ncol=3*nhour,byrow=TRUE)) -1
    ll.nw.pred<-predict(rob.ll.nw,newdata=cb.data,se.fit=TRUE,band="local")
    predmat.ll.nw<-exp(matrix(cbind(ll.nw.pred$fit,ll.nw.pred$fit-crit(rob.ll.nw)$crit.val*ll.nw.pred$se.fit,ll.nw.pred$fit+crit(rob.ll.nw)$crit.val*ll.nw.pred$se.fit),ncol=3*nhour,byrow=TRUE))-1
    # poisson
    cb.data$hour<-floor(cb.data$hour)
    cb.data$workday<-TRUE
    poi.w.pred<-predict(rob.poi,newdata=cb.data,type="link",se.fit=TRUE)
    predmat.poi.w<-exp(matrix(cbind(poi.w.pred$fit,poi.w.pred$fit-1.96*poi.w.pred$se.fit,poi.w.pred$fit+1.96*poi.w.pred$se.fit),ncol=3*nhour,byrow=TRUE))
    cb.data$workday<-FALSE
    poi.nw.pred<-predict(rob.poi,newdata=cb.data,type="link",se.fit=TRUE)
    predmat.poi.nw<-exp(matrix(cbind(poi.nw.pred$fit,poi.nw.pred$fit-1.96*poi.nw.pred$se.fit,poi.nw.pred$fit+1.96*poi.nw.pred$se.fit),ncol=3*nhour,byrow=TRUE))

    # make the plots
    pred.mats<-list(predmat.poi.w,predmat.ll.w,predmat.poi.nw,predmat.ll.nw)
    names(pred.mats)<-c("poi_w","ll_w","poi_nw","ll_nw")
    pred.mats$titles<-c("Poisson Model, Workdays","Local Linear Model, Workdays","Poisson Model, Weekend/Holiday","Local Linear Model, Weekend/Holiday")
    # parameters for plots
    max.pred<-max(pred.mats[[1]],pred.mats[[2]],pred.mats[[3]],pred.mats[[4]])
    min.pred<-min(pred.mats[[1]],pred.mats[[2]],pred.mats[[3]],pred.mats[[4]])
    ytic<-pretty(c(max.pred,min.pred))
    omit.3<-function(n) c(1:2,4:(n+1))
    for (i in 1:4) {
        pdf(paste("final_project/",names(pred.mats)[i],"_ci.pdf",sep=""))
        matplot(x=temp.seq,y=pred.mats[[i]],type="l",lty=c(omit.3(nhour),rep(3,nhour*2)),lwd=c(rep(1,nhour),rep(.5,2*nhour)),col=rep(1,nhour*3),xlab="Temperature",ylab="Mean Robberies Per Hour",main=paste("Predicted Rate of Robberies at Different Temperatures,",pred.mats$titles[i]),ylim=c(min(ytic),max(ytic)))
        legend('bottomright',c(paste(((hour.seq-1)%%12)+1,ifelse(hour.seq %/%12,"pm","am")),"95% CI"),lty=c(omit.3(nhour),3),lwd=c(rep(1,nhour),.5))
        dev.off()
    }



     
}
    
