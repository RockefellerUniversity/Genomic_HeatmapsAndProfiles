## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = T, eval=T)


## ----exerciseReadin, echo=T,eval=T------------------------------------------------------
library(ggplot2)

patients_clean <- read.delim("./data/patients_clean_ggplot2.txt",sep="\t")

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight))+geom_point()
plot


## ----exercise1--------------------------------------------------------------------------
plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight))+geom_point()
plot+scale_x_continuous(breaks=c(20,30,40),label=c(20,30,40),limits=c(20,40))+
  scale_y_continuous(breaks=seq(60,100,by=5),label=seq(60,100,by=5),
                     name="Weight (kilos)")


## ----exercise2--------------------------------------------------------------------------
plot <- ggplot(data=patients_clean,
               mapping=aes(x=factor(Age),y=BMI))+geom_violin(aes(fill=factor(Age)))+
               scale_fill_brewer(palette="Blues", na.value="black")


## ----exercise3--------------------------------------------------------------------------
library(ggplot2)

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()+
  scale_colour_gradient2(low="blue",high="yellow",mid="grey",midpoint=mean(patients_clean$Height))
plot


## ----exercise4--------------------------------------------------------------------------
library(ggplot2)

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()+
  scale_colour_gradient2(low="blue",high="yellow",
                         mid="grey",midpoint=mean(patients_clean$Height),
                         limits=c(min(patients_clean$Height),180),
                         na.value=NA)
plot


## ----exercise5--------------------------------------------------------------------------
library(ggplot2)

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()+
  scale_colour_gradient2(low="blue",high="yellow",
                         mid="grey",midpoint=mean(patients_clean$Height),
                         breaks=c(min(patients_clean$Height),
                                  median(patients_clean$Height),
                                  quantile(patients_clean$Height)[4]),
                         labels=c(signif(min(patients_clean$Height),3),
                                  signif(median(patients_clean$Height),3),
                                  signif(quantile(patients_clean$Height)[4],3)))
plot


## ----exercise6--------------------------------------------------------------------------

library(ggplot2)

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height,alpha=Count,size=Overweight))+geom_point()
plot


## ----exercise7--------------------------------------------------------------------------

library(ggplot2)

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=factor(Age)))+geom_point()+stat_smooth(method="lm",se=F)

plot


## ----exercise8--------------------------------------------------------------------------


plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=factor(Age)))+geom_point()+stat_smooth(method="lm",se=F)

plot <- plot+theme(legend.title=element_blank(),legend.background=element_rect(fill="white"),legend.key=element_rect(fill="white"),legend.position="bottom")

plot


## ----exercise9--------------------------------------------------------------------------


plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=factor(Age)))+geom_point()+stat_smooth(method="lm",se=F)

plot <- plot+theme(legend.title=element_blank(),legend.background=element_rect(fill="white"),legend.key=element_rect(fill="white"),legend.position="bottom")

plot <- plot+ggtitle("BMI vs Weight")+theme(panel.grid.minor=element_blank())

plot

ggsave(plot,file="BMIvsWeight.png",units = "in",height = 7,width = 7)

