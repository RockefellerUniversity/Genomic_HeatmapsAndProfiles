## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,eval=T)


## ----eval=TRUE,echo=FALSE,include=FALSE-------------------------------------------------
timeCourse <- read.delim("./data/timecourse.csv",sep=",",header = T) 
summary(timeCourse)


## ----eval=FALSE-------------------------------------------------------------------------
## timeCourse <- read.delim("./data/timecourse.csv",sep=",",header = T)
## summary(timeCourse)


## ---------------------------------------------------------------------------------------
hist(as.matrix(timeCourse)[,2:3],main = "Control",xlab = "Control")
hist(as.matrix(timeCourse)[,4:5],main = "Treatment",xlab = "Treatment")


## ---------------------------------------------------------------------------------------
par(mfrow=c(2,1))
hist(as.matrix(timeCourse)[,2:3],xlim = c(0,20),main = "Control",xlab = "Control")
hist(as.matrix(timeCourse)[,4:5],xlim = c(0,20),main = "Treatment",xlab = "Treatment")


## ---------------------------------------------------------------------------------------
par(mfrow=c(1,1))

control1 <- timeCourse[,2]
barplot(control1,names.arg = timeCourse[,2],xlab = "Time")



## ---------------------------------------------------------------------------------------
par(mfrow=c(1,1))

controlMeans <- rowMeans(timeCourse[,c(2,3)])
treatmentMeans <- rowMeans(timeCourse[,c(4,5)])
meanTable <- rbind(controlMeans,treatmentMeans)
barplot(meanTable,names.arg = timeCourse[,1],xlab = "Time",beside=TRUE,col=c("blue","red"))

legend("topleft",legend = c("Control","Treatment"),fill=c("blue","red"))



## ---------------------------------------------------------------------------------------
topOfY <- max(timeCourse[,-1])+3
 
plot(timeCourse[,2],type="o",xlab="Time",ylab="Score",pch=1,ylim=c(0,topOfY),col="blue") 
lines(timeCourse[,3],pch=1,col="blue",type="o") 
lines(timeCourse[,4],pch=1,col="red",type="o") 
lines(timeCourse[,5],pch=1,col="red",type="o") 
legend("topleft",legend = c("Control","Treatment"),col=c("blue","red"),lwd = 2)

