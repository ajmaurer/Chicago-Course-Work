library('xtable')
library('stats4')
library('locfit')
library('coefplot')
library('SDMTools')
library('plotrix')
library('combinat')
library('reshape2')

### Problem 1
# b
#I'm lazy, and this only invovles only hundreds of checks, so I'm just brute force going through all combinations to see which works best
all.comb<-combn(LETTERS[1:6],4)
first10<-combn(1:10,6) # choices of the 6 which inclue "a"
bvar<-10000000
bcomb<-NULL
bP<-NULL
bSadj<-NULL
for (n1 in 1:ncol(first10)) {
    Bct<-sum(all.comb[,first10[,n1]]=="B") # number of bs in first 6 chosen
    if (!(4>= Bct & Bct>=3)) next
    remain<-combn(11:14,6-Bct)
    if (Bct==4) remain <-rbind(remain,rep(15,ncol(remain))) # if we only pick 2 of 11-14, we must include 15
    for (n2 in 1:ncol(remain)) {
        comb<-all.comb[,c(first10[,n1],remain[,n2])]
        N<-matrix(0,nrow=9,ncol=6)
        for (i in 1:9) for (j in 1:6) N[i,j]<-sum(comb[,i]==LETTERS[j]) # occurances in each block 
        P<-t(N)%*%N
        Padj<-P-diag(diag(P))
        Sadj<-t(Padj)%*%Padj
        diag(P)<-NA
        diag(Sadj)<-NA
        rangeP<-range(as.vector(P),na.rm=T)
        if (rangeP[2]-rangeP[1]>1) next
        var<-var(as.vector(Sadj),na.rm=T)
        if (var<bvar) {
            bvar<-var
            bcomb<-comb
            bP<-P
            bSadj<-Sadj
        }
    }
}

rownames(bcomb)<-1:4
colnames(bcomb)<-paste('B',1:9,sep="")
rownames(bP)<-LETTERS[1:6]
colnames(bP)<-LETTERS[1:6]
rownames(bSadj)<-LETTERS[1:6]
colnames(bSadj)<-LETTERS[1:6]

print(xtable(bcomb,digits=0),type="latex",file="hw3/hw3_1a_comb.tex")
print(xtable(bP,digits=0),type="latex",file="hw3/hw3_1a_P.tex")
print(xtable(bSadj,digits=0),type="latex",file="hw3/hw3_1a_S.tex")

### problem 2
slp<-read.table("hw3/sleep.txt",sep="",header=F,col.names=c("per","trt","slp"))
# a
lm.all<-lm(slp~as.factor(trt)+as.factor(per),data=slp)
sig.all<-summary(lm.all)$sigma
lm.CD<-lm(slp~as.factor(trt)+as.factor(per),data=slp,subset=trt %in% c("C","D"))
sig.CD<-summary(lm.CD)$sigma
slp.rs<-dcast(melt(slp,id=c("per","trt")),per~trt)

pdf('hw3/hw3_2a_plot.pdf')
matplot(x=slp.rs[1]-.15,y=slp.rs[2],xaxp=c(1,10,9),pch="A",cex=.8,main="Outcome By Person and Treatment",xlab="Person",ylab="Sleep",ylim=c(0,90),col=1)
for (i in 2:4) {
    matplot(x=slp.rs[1]+i*.075-.225,y=slp.rs[i+1],pch=LETTERS[i],cex=.8,col=i,add=T)
}
dev.off()


