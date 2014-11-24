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
data.sum   <- TRUE 
par.reg    <- FALSE
nonpar.reg <- FALSE

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
    rm(hourly.weather,hourly.robberies)
    save(rob.data,file="final_project/hourly_robberies.RData")
    
}
###### Data Summary
if (data.sum) {
    if (!exits("rob.data")) {
        load(file="final_project/hourly_robberies.RData")
    }
    #tapply(robberies,list(hour,month),FUN=mean)
}
    
