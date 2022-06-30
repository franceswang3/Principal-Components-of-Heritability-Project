rm(list=ls()); gc()
library(readxl)
library(WriteXLS)
library(data.table)
library(haven)
library(nlme)
library(lme4)
library(lmerTest)
library(pscl)
library(tidyverse)
library(ordinal)
library(sjPlot)

install.packages("Matrix")
library(glmmTMB)

data<-read_sav("C:/Users/wangfl/OneDrive - University of Pittsburgh/_BoxMigration/PCH analyses/Merge with w9 sub data.sav")
#descriptives
hist(data$alcon, breaks=seq(-1,15,1))
hist(data$hvydrkOLD, breaks=seq(-1,6,1))

#Predict Alc Probs ZIP Model
summary(modelzip<-glmmTMB(alcon ~PCH.1+PCH.2+PCH.3+PCH.4+PCH.5+PCH.6+PCH.7+
                            agealcC+GENDER+(1|familyid), 
                          zi=~PCH.1+PCH.2+PCH.3+PCH.4+PCH.5+PCH.6+PCH.7+
                            agealcC+GENDER+(1|familyid),family=poisson, data=data))

#HEAVY DRINK
summary(modelzipHD<-glmmTMB(hvydrkOLD ~PCH.1+PCH.2+PCH.3+PCH.4+PCH.5+PCH.6+PCH.7+
                              agealcC+GENDER+(1|familyid), 
                          zi=~PCH.1+PCH.2+PCH.3+PCH.4+PCH.5+PCH.6+PCH.7+
                            agealcC+GENDER+(1|familyid),family=poisson, data=data))

tab_model(modelzip, modelzipHD, transform=NULL)

#SuperPCH correlations with Items

APCHn<-data[, c(7:19, 25:63)]
APCHn$disobeypar.pch<-APCHn$disobeypar.pch*-1
APCHn$disobeyschool.pch<-APCHn$disobeyschool.pch*-1
APCHn$badfriend.pch<-APCHn$badfriend.pch*-1
APCHn$temper.pch<-APCHn$temper.pch*-1
APCHn$secret.pch<-APCHn$secret.pch*-1
APCHn$guilty.pch<-APCHn$guilty.pch*-1
APCHn$rcadnofun.pch<-APCHn$rcadnofun.pch*-1
APCHn$atten.pch<-APCHn$atten.pch*-1
APCHn$g6act.1<-APCHn$g6act.1*-1
APCHn$g6atten.1<-APCHn$g6atten.1*-1

R=round(cor(APCHn[,-c(1:13)], APCHn[,1]),3)
items=rownames(R)
i=order(-abs(R))
superpch1items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch1items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch1items.csv")

R=round(cor(APCHn[,-c(1:13)], APCHn[,2]),3)
items=rownames(R)
i=order(-abs(R))
superpch2items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch2items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch2items.csv")

R=round(cor(APCHn[,-c(1:13)], APCHn[,3]),3)
items=rownames(R)
i=order(-abs(R))
superpch3items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch3items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch3items.csv")

R=round(cor(APCHn[,-c(1:13)], APCHn[,4]),3)
items=rownames(R)
i=order(-abs(R))
superpch4items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch4items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch4items.csv")

R=round(cor(APCHn[,-c(1:13)], APCHn[,5]),3)
items=rownames(R)
i=order(-abs(R))
superpch5items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch5items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch5items.csv")

R=round(cor(APCHn[,-c(1:13)], APCHn[,6]),3)
items=rownames(R)
i=order(-abs(R))
superpch6items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch6items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch6items.csv")

R=round(cor(APCHn[,-c(1:13)], APCHn[,7]),3)
items=rownames(R)
i=order(-abs(R))
superpch7items<-data.frame(item=items[i],correlation=R[i])[1:40,]
write.csv(superpch7items, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/superpch7items.csv")

#linear regression predicting each pch
rawdata<-read.table("C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/Paper analysis documents/PCH-ResultsPackage-08-13-2021/Step-1/Output/std_phenotypes_07-27-2021.txt", header=TRUE)
superpchdata<-data[,c(1, 7:19)]
regdata<-merge(superpchdata, rawdata, by="ISN")

trait<-c("mean.1","disobeypar.1","disobeyschool.1","noguilt.1","rules.1","badfriend.1","lie.1","olderkids.1","swear.1",
         "tease.1","temper.1","argue.1","atten.1","stubborn.1","mood.1","susp.1","loud.1","eat.1","perfect.1","alone.1",
         "guilty.1","overtired.1","secret.1","sleepless.1","sleepmore.1","troublesleep.1","noenergy.1","sad.1","noinvolve.1",
         "rcadsad.1","rcadeat.1","rcadmove.1","rcadnofun.1","rcadtired.1", "g6frus.1", "g6fear.1",  "g6atten.1",  "g6act.1",  "g6surg.1")

for (trt in trait){
  formula2=as.formula(paste(trt," ~ AGE.1+GENDER"))
  model2<-lm(formula2, data=regdata)
  name<-paste(trt, "rz", sep="")
  regdata[[name]]<-scale(model2$residuals)
}
summary(lm(PCH.1~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))
summary(lm(PCH.2~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))
summary(lm(PCH.3~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))
summary(lm(PCH.4~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))
summary(lm(PCH.5~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))
summary(lm(PCH.6~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))
summary(lm(PCH.7~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.1rz+disobeypar.1rz+disobeyschool.1rz+noguilt.1rz+rules.1rz+badfriend.1rz+lie.1rz+olderkids.1rz+swear.1rz+tease.1rz+temper.1rz+
             argue.1rz+atten.1rz+stubborn.1rz+mood.1rz+susp.1rz+loud.1rz+eat.1rz+perfect.1rz+alone.1rz+guilty.1rz+overtired.1rz+secret.1rz+sleepless.1rz+sleepmore.1rz+troublesleep.1rz+noenergy.1rz+sad.1rz+noinvolve.1rz+
             rcadsad.1rz+rcadeat.1rz+rcadmove.1rz+rcadnofun.1rz+rcadtired.1rz, data=regdata))


trait<-c("mean.2","disobeypar.2","disobeyschool.2","noguilt.2","rules.2","badfriend.2","lie.2","olderkids.2","swear.2",
         "tease.2","temper.2","argue.2","atten.2","stubborn.2","mood.2","susp.2","loud.2","eat.2","perfect.2","alone.2",
         "guilty.2","overtired.2","secret.2","sleepless.2","sleepmore.2","troublesleep.2","noenergy.2","sad.2","noinvolve.2",
         "rcadsad.2","rcadeat.2","rcadmove.2","rcadnofun.2","rcadtired.2", "g6frus.1", "g6fear.1",  "g6atten.1",  "g6act.1",  "g6surg.1")

rsquared <- data.frame(matrix(ncol = 0, nrow = 1))
for (trt in trait){
  formula2=as.formula(paste(trt," ~ AGE.2+GENDER"))
  model2<-lm(formula2, data=regdata)
  name<-paste(trt, "rz", sep="")
  regdata[[name]]<-scale(model2$residuals)
}
rsquared$r21<-summary(lm(PCH.1~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquared$r22<-summary(lm(PCH.2~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquared$r23<-summary(lm(PCH.3~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquared$r24<-summary(lm(PCH.4~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquared$r25<-summary(lm(PCH.5~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquared$r26<-summary(lm(PCH.6~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquared$r27<-summary(lm(PCH.7~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                           argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                           rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))$r.squared
rsquaredT2PredPCH<-write.csv(rsquared, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/rsquaredT2PredPCH2-25-22.csv")

coefs=list()
coef1<-summary(lm(PCH.1~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))

coef2<-summary(lm(PCH.2~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))
coef3<-summary(lm(PCH.3~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))
coef4<-summary(lm(PCH.4~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))
coef5<-summary(lm(PCH.5~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))
coef6<-summary(lm(PCH.6~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))
coef7<-summary(lm(PCH.7~g6frus.1rz+g6fear.1rz+g6atten.1rz+g6act.1rz+g6surg.1rz+mean.2rz+disobeypar.2rz+disobeyschool.2rz+noguilt.2rz+rules.2rz+badfriend.2rz+lie.2rz+olderkids.2rz+swear.2rz+tease.2rz+temper.2rz+
                    argue.2rz+atten.2rz+stubborn.2rz+mood.2rz+susp.2rz+loud.2rz+eat.2rz+perfect.2rz+alone.2rz+guilty.2rz+overtired.2rz+secret.2rz+sleepless.2rz+sleepmore.2rz+troublesleep.2rz+noenergy.2rz+sad.2rz+noinvolve.2rz+
                    rcadsad.2rz+rcadeat.2rz+rcadmove.2rz+rcadnofun.2rz+rcadtired.2rz, data=regdata))
coefpch1<-(as.data.frame(coef1[["coefficients"]]))
write.csv(coefpch1, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch1 2-25-22.csv")
coefpch2<-(as.data.frame(coef2[["coefficients"]]))
write.csv(coefpch2, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch2 2-25-22.csv")
coefpch3<-(as.data.frame(coef3[["coefficients"]]))
write.csv(coefpch3, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch3 2-25-22.csv")
coefpch4<-(as.data.frame(coef4[["coefficients"]]))
write.csv(coefpch4, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch4 2-25-22.csv")
coefpch5<-(as.data.frame(coef5[["coefficients"]]))
write.csv(coefpch5, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch5 2-25-22.csv")
coefpch6<-(as.data.frame(coef6[["coefficients"]]))
write.csv(coefpch6, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch6 2-25-22.csv")
coefpch7<-(as.data.frame(coef7[["coefficients"]]))
write.csv(coefpch7, "C:/Users/DR.  Frances Wang/Box Sync/PCH analyses/coefpch7 2-25-22.csv")

