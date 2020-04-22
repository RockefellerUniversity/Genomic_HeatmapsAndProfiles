params <-
list(isSlides = "no")

## ----setup_basePlotting, include=FALSE---------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ---- results='asis',include=TRUE,echo=FALSE---------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Plotting in R

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 


"    
  )
}else{
  cat("# Plotting in R


"    
  )
  
}



## ----setwd_introtoR,eval=F---------------------------------------
## setwd("/PathToMyDownload/Plotting_In_R/r_course")
## # e.g. setwd("~/Downloads/Intro_To_R_1Day/r_course")


## ----treatmentVector_basePlotting--------------------------------
treatment <- c(0.02,1.8, 17.5, 55,75.7, 80)


## ----defaultPlotEval_basePlotting,fig.width=8,fig.height=4.5,dpi=300,out.width="1920px",height="1080px"----
plot(treatment)


## ----plotType_basePlotting,eval=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o")


## ----plotTypeSEval_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l")


## ----plotTypeS_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="p")


## ----plotLabels_basePlotting,eval=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, main="My Plot", sub="a plot")


## ----plotAsisLabels_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, xlab="Position", ylab="score")


## ----plotAsisLasL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, las=1)


## ----plotAsisLasR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, las=2)


## ----plotAsisCexL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, cex=2)


## ----plotAsisCexR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, cex=0.5)


## ----plotAsisPchL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, pch=1)


## ----plotAsisPchR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, pch=20)


## ----plotAsisLwdL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l",lwd=10)


## ----plotAsisLwdR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l",lwd=0.5)


## ----plotAsisLtyL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l",lty=1)


## ----plotAsisLtyR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l",lty=2)


## ----plotColL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l", col="red")


## ----plotColR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="l", col="dodgerblue")


## ----plotTwoVectors_basePlotting,eval=TRUE,collapse=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1000px",height="1080px"----
control <- c(0, 20, 40, 60, 80,100)
plot(treatment,control)


## ----lines_basePlotting,eval=TRUE,collapse=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o", col="blue")
lines(control, type="o", pch=22, lty=2, col="red")


## ----yxlim_basePlotting,eval=TRUE,collapse=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1000px",height="1080px"----
control <- c(0, 20, 40, 60, 80,100)
plot(treatment, type="o", col="blue",ylim=c(0,100))
lines(control, type="o", pch=22, lty=2, col="red")


## ----range_basePlotting------------------------------------------
g_range <- range(0, treatment, control)
g_range


## ----plotNoAxesNoAnnotation_basePlotting,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o", col="blue", 
     ylim=g_range, axes=FALSE, ann=FALSE)


## ----axis_basePlotting,eval=FALSE--------------------------------
## axis(side=1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))


## ----axisEval_basePlotting,eval=TRUE,echo=FALSE,fig.width=8,fig.height=4,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o", col="blue", 
     ylim=g_range,axes=FALSE, ann=FALSE)
axis(1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))



## ----axisSides_basePlotting,eval=FALSE---------------------------
## axis(2, las=1, at=rep(0,g_range[2],by=20))


## ----axisSidesEval_basePlotting,eval=TRUE,echo=FALSE,fig.height=3.5,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o", col="blue", 
     ylim=g_range,axes=FALSE, ann=FALSE)
axis(1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))
axis(2, las=1, at=20*0:g_range[2])



## ----box_basePlotting,eval=FALSE---------------------------------
## box()


## ----boxEval_basePlotting,eval=TRUE,echo=FALSE,fig.height=4,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o", col="blue", 
     ylim=g_range,axes=FALSE, ann=FALSE)
axis(1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))
axis(2, las=1, at=20*0:g_range[2])
box()


## ----addLines_basePlotting,eval=FALSE,echo=TRUE,collapse=TRUE----
## lines(control, type="o", pch=22, lty=2, col="red")
## 


## ----addLinesEval_basePlotting,eval=TRUE,echo=FALSE,fig.width=8,fig.height=4,dpi=300,out.width="1000px",height="1080px"----
plot(treatment, type="o", col="blue", 
     ylim=g_range,axes=FALSE, ann=FALSE)
axis(1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))
axis(2, las=1, at=20*0:g_range[2])
box()
lines(control, type="o", pch=22, lty=2, col="red")


## ----legend_basePlotting,eval=FALSE------------------------------
## legend("topleft",legend=c("treatment","control"),
##        col=c("blue","red"), pch=21:22, lty=1:2);


## ----legendEval_basePlotting,echo=FALSE,fig.width=8,fig.height=3.5,dpi=300,out.width="1920px",height="1080px"----
plot(treatment, type="o", col="blue", ylim=g_range,axes=FALSE, ann=FALSE)
axis(1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))
axis(2, las=1, at=20*0:g_range[2])
box()

lines(control, type="o", pch=22, lty=2, col="red")
legend("topleft",legend=c("treatment","control"),
       col=c("blue","red"), pch=21:22, lty=1:2);  


## ----readability_basePlotting,echo=FALSE,fig.width=8,fig.height=3.5,dpi=300,out.width="1920px",height="1080px"----

plot(treatment, type="o", col="blue", lwd=1, ylim=g_range,axes=FALSE, ann=FALSE)
axis(1, at=1:6, lab=c("Mon","Tue","Wed","Thu","Fri","Sat"))
axis(2, las=1, at=20*0:g_range[2])
box()

lines(control, type="o", pch=22, lty=2, col="red", lwd=2.5)
legend("topleft",legend=c("treatment","control"),col=c("blue","red"), pch=21:22, lty=1:2, lwd=c(1,2.5))


## ---- echo=F,warning=F,message=F---------------------------------
library(viridis)
library(scales)


## ---- eval=F-----------------------------------------------------
## install.packages('viridis')
## library(viridis)
## viridis(5)

## ---- echo=F-----------------------------------------------------
viridis(5)


## ---- echo=F-----------------------------------------------------
show_col(viridis(5))


## ----barplotEval_basePlotting,fig.width=8,fig.height=4,dpi=300,out.width="1920px",height="1080px"----
barplot(treatment)


## ----barplotNames_basePlotting,eval=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1920px",height="1080px"----
barplot(treatment,
        names.arg=c("Mon","Tue","Wed","Thu","Fri","Sat"))


## ----barplotPreNames_basePlotting,eval=TRUE,fig.width=8,fig.height=3.5,dpi=300,out.width="1920px",height="1080px"----
names(treatment) <- c("Mon","Tue","Wed","Thu","Fri","Sat")
barplot(treatment)


## ----readTable_basePlotting--------------------------------------
data <- read.table("data/example_plot.txt", header=T, row.names=1, sep=",")


## ----barplotMatrix_basePlotting,eval=TRUE,fig.width=8,fig.height=4,dpi=300,out.width="1920px",height="1080px"----
barplot(as.matrix(data))


## ----barplotMatrixBeside_basePlotting,eval=TRUE,fig.width=8,fig.height=4,dpi=300,out.width="1920px",height="1080px"----
barplot(as.matrix(data),beside=TRUE)


## ----barplotComplexEval_basePlotting, fig.width=8,fig.height=3,dpi=300,out.width="1920px",height="1080px"----
barplot(as.matrix(data), main="Daily progression of X in\nControl and Treatment", ylab= "Total", beside=TRUE, 
        col= viridis(6))
legend("topleft", c("Mon","Tue","Wed","Thu","Fri","Sat"), cex=0.8,
        fill=  viridis(6))


## ----hist_basePlotting,echo=TRUE,fig.width=8,fig.height=3,dpi=300,out.width="1920px",height="1080px"----
hist(treatment)  


## ----histcolorYlim_basePlotting,echo=TRUE,fig.width=8,fig.height=3,dpi=300,out.width="1920px",height="1080px"----
hist(treatment, col="lightblue", ylim=c(0,5),cex.main=0.8)


## ----histBreaksEval_basePlottingL_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
hist(treatment, col="lightblue", 
     ylim=c(0,5), cex.main=0.8, 
     breaks = 2)


## ----histBreaksEval_basePlottingR_basePlotting,eval=TRUE,fig.width=8,fig.height=4.5,dpi=300,out.width="1000px",height="1080px"----
hist(treatment, col="lightblue", 
     ylim=c(0,5), cex.main=0.8, 
     breaks = 10)


## ----dotchartEval_basePlotting,fig.width=8,fig.height=3,dpi=300,out.width="1500px",height="750px"----
dotchart(t(data))  


## ----dotchartColEval_basePlotting,fig.width=6,fig.height=3,dpi=300,out.width="1920px",height="1080px"----

dotchart(t(data), color=c("red","blue"),main="Dotchart", cex=0.5)



## ----readForBoxplot_basePlotting,eval=T,fig.width=6,fig.height=3,dpi=300,out.width="1920px",height="1080px"----
exprs <- read.delim("data/gene_data.txt",sep="\t",h=T,row.names = 1)
head(exprs)


## ----Boxplot_basePlotting,eval=T,fig.width=6,fig.height=3,dpi=300,out.width="1920px",height="1080px"----
boxplot(exprs)


## ----BoxplotEval_basePlotting,eval=T,fig.width=6,fig.height=3,dpi=300,out.width="1920px",height="1080px"----
boxplot(log2(exprs),ylab="log2 Expression",
        col=c("red","red","blue","blue"))


## ----parMfrow_basePlotting,eval=FALSE----------------------------
## par(mfrow=c(2,2))


## ----setParMfrowEval_basePlotting,fig.width=8,fig.height=5,dpi=300,out.width="490px",height="270px"----
par(mfrow=c(2,2))
hist(exprs$Untreated1)
hist(exprs$Untreated2)
hist(exprs$Treated1)
hist(exprs$Treated2)


## ----setParMfrowEval_for_basePlotting,fig.width=8,fig.height=5,dpi=300,out.width="490px",height="270px"----
par(mfrow=c(2,2))
for (i in 1:4){
hist(exprs[,i])
}


## ---- echo=F-----------------------------------------------------
par(mfrow=c(1,1))



## ----fig.width=8,fig.height=5,dpi=300,out.width="490px",height="270px"----
plot(control, treatment)
text(20,60, 'THIS IS MY PLOT', col='red')
text(control, treatment, letters[1:6], adj=c(0,-1), col='blue')



## ----fig.width=8,fig.height=5,dpi=300,out.width="490px",height="270px"----
plot(control, treatment)
abline(h=10, col='blue')
abline(v=50, col='red', lwd=2)
abline(a=0, b=1, lty=2)



## ----fig.width=8,fig.height=5,dpi=300,out.width="490px",height="270px"----
plot(control, treatment)
polygon(c(50,50,100,100),c(50,80,80,50), col='gray', density=5)


## ----save_basePlotting,eval=FALSE--------------------------------
## bmp(filename, width = 480, height = 480, units = "px",
##     point-size = 12)
## jpeg(filename, width = 480, height = 480, units = "px",
##      point-size = 12, quality = 75)


## ----save_bmp_basePlotting,eval=FALSE----------------------------
## bmp(file = "control.bmp")
## plot(control)
## dev.off()


## ----save_jpeg_basePlotting,eval=FALSE---------------------------
## jpeg(file = "control.jpg", quality = 20)
## plot(control)
## dev.off()


## ----save_ps_basePlotting,eval=FALSE-----------------------------
## postscript(file = "control.ps")
## plot(control)
## dev.off()


## ----save_pdf_basePlotting,eval=FALSE----------------------------
## 
## pdf(file = "control.pdf", paper = "A4")
## plot(control)
## dev.off()

