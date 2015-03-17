library('xtable')
library('stats4')
library('locfit')
library('coefplot')
library('SDMTools')
library('plotrix')
d=3

### Problem 2
o<-"O"
a<-"A"
b<-"B"
c<-"C"
blk<-as.factor(c(rep(4.1,4),rep(4.2,4),rep(3.1,3),rep(3.2,3),2.1,2.1,2.2,2.2,2.3,2.3))

# a
trt<-as.factor(c(o,a,b,c, o,a,b,c, a,b,c, a,b,c, o,a, o,b, o,c))
trt.2a<-trt
#trt.2a<-as.factor(c(o,a,b,c, o,a,b,c, o,a,b, c,o,a, b,c, o,b, a,c))
mm.2a<-model.matrix(~trt+blk)
cov.2a<-solve(t(mm.2a)%*%mm.2a)
tc<-cov.2a[2:4,2:4]
dif.2a<-matrix(c(diag(tc),tc[1,1]+tc[2,2]-2*tc[2,1],tc[1,1]+tc[3,3]-2*tc[3,1],tc[2,2]+tc[3,3]-2*tc[3,2]),nrow=1)
colnames(dif.2a)<-c("A-B","A-C","A-O","B-C","B-O","C-O")
rownames(dif.2a)<-"Unscaled Variance"
print(xtable(cov.2a,digits=d),type="latex",file="hw2/hw2_2a_cov.tex")
print(xtable(dif.2a,digits=d),type="latex",file="hw2/hw2_2a_dif_var.tex")

#b
trt<-as.factor(c(o,a,b,c, o,o,a,b, o,a,c, o,b,c, o,a, o,b, o,c))
trt.2b<-trt
mm.2b<-model.matrix(~trt+blk)
cov.2b<-solve(t(mm.2b)%*%mm.2b)
tc<-cov.2b[2:4,2:4]
dif.2b<-matrix(c(diag(tc),tc[1,1]+tc[2,2]-2*tc[2,1],tc[1,1]+tc[3,3]-2*tc[3,1],tc[2,2]+tc[3,3]-2*tc[3,2]),nrow=1)
colnames(dif.2b)<-c("A-B","A-C","A-O","B-C","B-O","C-O")
rownames(dif.2a)<-"Unscaled Variance"
print(xtable(cov.2b,digits=d),type="latex",file="hw2/hw2_2b_cov.tex")
print(xtable(as.matrix(dif.2b),digits=d),type="latex",file="hw2/hw2_2b_dif_var.tex")

### 3 
wht<-read.table("hw2/wheat56.dat",sep="",header=T)
wht$Block <- as.factor(rep(1:4,each=56))
wht$Variety <- as.factor(wht$Variety)
attr(wht$Variety,"contrasts")<-"contr.sum"
n<-nrow(wht)

# a
wht.lm<-lm(Yield~Variety + Block,data=wht)
print(xtable(summary(wht.lm),digits=d),type="latex",file="hw2/hw2_3a_lin_reg.tex")
print(xtable(anova(wht.lm),digits=d),type="latex",file="hw2/hw2_3a_anova.tex")

pdf("hw2/hw2_3a_coef_plot.pdf")
coefplot(wht.lm,outerCI=1.96,innerCI=1.96,color=1,decreasing=T,predictors="Variety",title="Variety Coef and 95% CIs")
dev.off()

# b
wht.cf<-coef(wht.lm)[2:56]
wht.vcov<-vcov(wht.lm)[2:56,2:56]
a<-as.matrix(c(rep(1/20+1/36,20),rep(-1/36+1/36,35)))
var<-t(a) %*% wht.vcov %*% a
est<- t(a) %*% wht.cf
ci<-est+c(qt(.025,165),qt(.975,165))*sqrt(var)

# c

pdf('hw2/hw2_3c_plot.pdf')
rtob<-function(x) rgb(x,0,1-x)
plot(0,ylim=c(4,48),xlim=c(1,27),xlab="Longitude",ylab="Latitude",main="Residual Plot")
mtext("Variety Labeled")
text(wht$Longitude,wht$Latitude,labels=wht$Variety,col=rtob(residuals(wht.lm)/850+.5))
ulat<-sort(unique(wht$Latitude))
ulon<-sort(unique(wht$Longitude))
dlat<-4.3/2
dlon<-1.2/6
x0<-     c(0      ,ulon[6],ulon[6],0       ,ulon[19],ulon[19],0       ,ulon[10],ulon[10])
x1<-     c(ulon[6],ulon[6],28     ,ulon[19],ulon[19],28      ,ulon[10],ulon[10],28      )
y0<-dlat+c(ulat[4],ulat[4],ulat[3],ulat[6] ,ulat[6] ,ulat[5] ,ulat[9] ,ulat[9] ,ulat[8] )
y1<-dlat+c(ulat[4],ulat[3],ulat[3],ulat[6] ,ulat[5] ,ulat[5] ,ulat[9] ,ulat[8] ,ulat[8] )
segments(x0,y0,x1,y1,lty=2)
text(rep(27.5,4),ulat[c(2,5,7,10)]-c(0,dlat,0,0),paste("Block",1:4),srt=270,cex=.8)
color.legend(xl=1,yb=3,xr=12,yt=4,legend=c(-500,-250,0,250,500),rect.col=rtob(0:64/64),gradient="x",cex=.65)
segments(x0=c(0,13),y0=c(dlat+ulat[1],dlat+ulat[1]),x1=c(13,13),y1=c(dlat+ulat[1],0))
text(6.5,5.5,"Residual",cex=.8)
dev.off()

#d

wht2.lm<-lm(Yield~Variety + Block+Latitude+Longitude,data=wht)

pdf("hw2/hw2_3d_coef_plot.pdf")
coefplot(wht2.lm,outerCI=1.96,innerCI=1.96,color=1,decreasing=T,predictors="Variety",title="Variety Coef and 95% CIs, Including Linear Function Coordinates")
dev.off()

pdf('hw2/hw2_3d_plot.pdf')
rtob<-function(x) rgb(x,0,1-x)
plot(0,ylim=c(4,48),xlim=c(1,27),xlab="Longitude",ylab="Latitude",main="Yield Residual Plot, Linear Function Coordinates Included")
mtext("Variety Labeled")
text(wht$Longitude,wht$Latitude,labels=wht$Variety,col=rtob(residuals(wht2.lm)/850+.5))
ulat<-sort(unique(wht$Latitude))
ulon<-sort(unique(wht$Longitude))
dlat<-4.3/2
dlon<-1.2/6
x0<-     c(0      ,ulon[6],ulon[6],0       ,ulon[19],ulon[19],0       ,ulon[10],ulon[10])
x1<-     c(ulon[6],ulon[6],28     ,ulon[19],ulon[19],28      ,ulon[10],ulon[10],28      )
y0<-dlat+c(ulat[4],ulat[4],ulat[3],ulat[6] ,ulat[6] ,ulat[5] ,ulat[9] ,ulat[9] ,ulat[8] )
y1<-dlat+c(ulat[4],ulat[3],ulat[3],ulat[6] ,ulat[5] ,ulat[5] ,ulat[9] ,ulat[8] ,ulat[8] )
segments(x0,y0,x1,y1,lty=2)
text(rep(27.5,4),ulat[c(2,5,7,10)]-c(0,dlat,0,0),paste("Block",1:4),srt=270,cex=.8)
color.legend(xl=1,yb=3,xr=12,yt=4,legend=c(-500,-250,0,250,500),rect.col=rtob(0:64/64),gradient="x",cex=.65)
segments(x0=c(0,13),y0=c(dlat+ulat[1],dlat+ulat[1]),x1=c(13,13),y1=c(dlat+ulat[1],0))
text(6.5,5.5,"Residual",cex=.8)
dev.off()

#e

wht3.lm<-lm(Yield~Variety + Block+poly(Latitude,Longitude,degree=2),data=wht)

pdf("hw2/hw2_3e_coef_plot.pdf")
coefplot(wht3.lm,outerCI=1.96,innerCI=1.96,color=1,decreasing=T,predictors="Variety",title="Variety Coef and 95% CIs, Including Linear Function Coordinates")
dev.off()

pdf('hw2/hw2_3e_plot.pdf')
rtob<-function(x) rgb(x,0,1-x)
plot(0,ylim=c(4,48),xlim=c(1,27),xlab="Longitude",ylab="Latitude",main="Residual Plot, Quadratic Function Coordinates Included")
mtext("Variety Labeled")
text(wht$Longitude,wht$Latitude,labels=wht$Variety,col=rtob(residuals(wht3.lm)/850+.5))
ulat<-sort(unique(wht$Latitude))
ulon<-sort(unique(wht$Longitude))
dlat<-4.3/2
dlon<-1.2/6
x0<-     c(0      ,ulon[6],ulon[6],0       ,ulon[19],ulon[19],0       ,ulon[10],ulon[10])
x1<-     c(ulon[6],ulon[6],28     ,ulon[19],ulon[19],28      ,ulon[10],ulon[10],28      )
y0<-dlat+c(ulat[4],ulat[4],ulat[3],ulat[6] ,ulat[6] ,ulat[5] ,ulat[9] ,ulat[9] ,ulat[8] )
y1<-dlat+c(ulat[4],ulat[3],ulat[3],ulat[6] ,ulat[5] ,ulat[5] ,ulat[9] ,ulat[8] ,ulat[8] )
segments(x0,y0,x1,y1,lty=2)
text(rep(27.5,4),ulat[c(2,5,7,10)]-c(0,dlat,0,0),paste("Block",1:4),srt=270,cex=.8)
color.legend(xl=1,yb=3,xr=12,yt=4,legend=c(-500,-250,0,250,500),rect.col=rtob(0:64/64),gradient="x",cex=.65)
segments(x0=c(0,13),y0=c(dlat+ulat[1],dlat+ulat[1]),x1=c(13,13),y1=c(dlat+ulat[1],0))
text(6.5,5.5,"Residual",cex=.8)
dev.off()




