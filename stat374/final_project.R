###### Header
library('parallel')
library('locfit')
library('gam')
library('MASS')
library('xtable')
library('reshape2')
library('chron')

# Turns on different parts of the code
clean.data <- FALSE 
data.sum   <- FALSE
par.reg    <- FALSE 
nonpar.reg <- FALSE
eval.fit   <- TRUE 

# par.reg options
cross.validate <- FALSE

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
    matplot(tb.w,type="l",col=rep(1,4),xlab="Hour (Military Time)",ylab="Mean Robberies Per Hour",main="Rate of Robberies at Different Temperatures, Workdays",xaxp=c(0,24,8),ylim=c(0,3),yaxp=c(0,3,6))
    legend('bottomright',c("Below 30 Degrees","30-50 Degrees","50-70 Degrees","Above 70 Degrees"),lty=1:4)
    dev.off()

    tb.nw<-tapply(rd.nw$robberies,list(rd.nw$hour,rd.nw$temp.buc),FUN=mean)
    pdf("final_project/tempbuc_nw.pdf")
    matplot(tb.nw,type="l",col=rep(1,4),xlab="Hour (Military Time)",ylab="Mean Robberies Per Hour",main="Rate of Robberies at Different Temperatures, Weekend/Holiday",xaxp=c(0,24,8),ylim=c(0,3),yaxp=c(0,3,6))
    legend('bottomright',c("Below 30 Degrees","30-50 Degrees","50-70 Degrees","Above 70 Degrees"),lty=1:4)
    dev.off()

    # Effect of precipitation
    precip.w<-ifelse(rob.data$posprecip,ifelse(rob.data$workday,1,2),ifelse(rob.data$workday,3,4)) 
    precip.tab<-tapply(rob.data$robberies,list(rob.data$hour,precip.w),FUN=mean)
    pdf("final_project/precip.pdf")
    matplot(precip.tab,type="l",col=rep(1,4),xlab="Hour (Military Time)",ylab="Mean Robberies Per Hour",main="Rate of Robberies With/Without Precipitation",xaxp=c(0,24,8),ylim=c(0,3),yaxp=c(0,3,6),lty=c(1,3,2,4))
    legend('bottomright',c("Precipitation, Workday","Precipitation, Weekend/Holiday","No Precipitation, Workday","No Precipitation, Weekend/Holiday"),lty=c(1,3,2,4))
    dev.off()

}

###### Parametric Model
if (par.reg) {
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }

    ### Poisson model, selected based on Akaike Criteria
    rob.poi<-glm(glm(robberies~as.factor(hour)+as.factor(floor(hour/3))*workday+as.factor(floor(hour/3))*log(temp+10),family=poisson(link=log),data=rob.data,subset=test.set))

    save(rob.poi,file="final_project/poisson_model.RData")
}

###### Non-Parametric Model
if (nonpar.reg) {
    if (!exists("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }
    
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

    # Want to select bandwith and scale for second variable (effectively two bandwiths). Test a grid of values and pick the best one. Don't care about hours outside 0-23. 
    if (cross.validate) {
        sub<-rd$hour>=0 & rd$hour<24
        cv.score<-function(nn,s2,data,subset) {
            cv.fit<-fitted(rob.ll<-locfit(logrob~lp(temp,hour,deg=1,nn=nn,scale=c(1,s2)),kern="gauss",data=data),what="coef",cv=TRUE)[subset]
            return(sum((data$logrob[subset]-cv.fit)^2)/sum(subset))
        }
        nnvals<-seq(.02,.06,.005)
        svals<-seq(.02,.1,.01)
        cv.mat<-NULL
        for (s in svals) {
            cv.mat<-cbind(cv.mat,mcmapply(cv.score,nn=nnvals,MoreArgs=list(s2=s,data=rd,subset=sub),mc.cores=4))
        }
        ind<-which(cv.mat==min(cv.mat),arr.ind=TRUE)
        nnval<-nnvals[ind[1]]
        s2<-svals[ind[2]]
    }
    else {
        # These values were found from previous cross validation
        nnval<-.045
        s2<-.06
    }

    # Fit local linear model 
    rob.ll<-locfit(logrob~lp(temp,hour,deg=1,nn=nnval,scale=c(1,s2)),kern="gauss",data=rd)

    # Confidence intervals
    #crit(rob.ll,...)

    save(rob.ll,file="final_project/local_linear_model.RData")

}
if (eval.fit) {

}
    
