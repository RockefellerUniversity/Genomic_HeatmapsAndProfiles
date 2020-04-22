## ----exerciseReadin, echo=T-------------------------------------------------------------
getwd()
patients_clean <- read.delim("./data/patients_clean_ggplot2.txt",sep="\t")


## ----exercise1--------------------------------------------------------------------------
library(ggplot2)

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight))+geom_point()
plot


## ----exercise2--------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()
plot


## ----exercise3--------------------------------------------------------------------------
plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()
plot+facet_grid(Sex~Smokes)


## ----exercise4--------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()+
  geom_smooth()
plot


## ----exercise5--------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(x=BMI,y=Weight,colour=Height))+geom_point()+
  geom_smooth(method="lm",se=F)
plot


## ----exercise6--------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(x=Smokes,y=BMI))+geom_boxplot()
plot


## ----exercise7--------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(x=Smokes,y=BMI,colour=Sex))+geom_boxplot()

plot


## ----exercise8--------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(x=Smokes,y=BMI,colour=Sex))+
  geom_boxplot()+
  facet_wrap(~Age)
plot


## ----exercise9--------------------------------------------------------------------------
plot <- ggplot(data=patients_clean,
               mapping=aes(x=Sex,y=BMI,colour=factor(Age)))+
  geom_boxplot()+
  facet_wrap(~Smokes)
plot


## ----exercise10-------------------------------------------------------------------------
plot <- ggplot(data=patients_clean,
               mapping=aes(x=Sex,y=BMI,colour=factor(Age)))+
  geom_violin()+
  facet_wrap(~Smokes)
plot


## ----exercise11-------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(BMI))+
  geom_histogram(fill="blue")
plot


## ----exercise12-------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(BMI))+ geom_density(aes(fill=Sex),alpha=0.5)
plot


## ----exercise13-------------------------------------------------------------------------

plot <- ggplot(data=patients_clean,
               mapping=aes(BMI))+ geom_density(aes(fill=Sex),alpha=0.5)
plot+facet_wrap(~Grade)

