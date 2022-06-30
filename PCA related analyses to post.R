rm(list=ls()); gc()
library(readxl)
library(WriteXLS)
library(data.table)
library(haven)
library(lmerTest)
library(sjPlot)
library(glmmTMB)
library(nlme)
library(lme4)
library(lmerTest)
pcadata <- read_excel("C:/Users/wangfl/Box Sync/PCH analyses/Paper analysis documents/PCH-ResultsPackage-08-13-2021/Step-3/Output/VCE-SOLN-PCH-T1-T2.xlsx")
pcadata <- data.frame(pcadata, row.names = 1)
pcadata$mean.rij<-as.numeric(pcadata$mean.rij)
data<- read.table("C:/Users/wangfl/Box Sync/PCH analyses/Paper analysis documents/PCH-ResultsPackage-08-13-2021/Step-1/Output/std_phenotypes_07-27-2021.txt", header=TRUE)

#Residualize T1 and T2 CP and DEP on age and gender and scale and put into data
trait<-c("mean","disobeypar","disobeyschool","noguilt","rules","badfriend","lie","olderkids","swear",
         "tease","temper","argue","atten","stubborn","mood","susp","loud","eat","perfect","alone",
         "guilty","overtired","secret","sleepless","sleepmore","troublesleep","noenergy","sad","noinvolve",
         "rcadsad","rcadeat","rcadmove","rcadnofun","rcadtired")

for (trt in trait){
  trt1<-paste(trt,".1", sep="")
  trt2<-paste(trt,".2", sep="")
  formula=as.formula(paste(trt1," ~ AGE.1+GENDER"))
  model1<-lm(formula, data=data)
  formula2=as.formula(paste(trt2," ~ AGE.2+GENDER"))
  model2<-lm(formula2, data=data)
  name1<-paste(trt1, "rz", sep="")
  name2<-paste(trt2, "rz", sep="")
  data[[name1]]<-scale(model1$residuals)
  data[[name2]]<-scale(model2$residuals)
}


#create PCA of each item (combine T1 and T2), using standardized raw data and G+R matrix created using standardized data, load it into datafile, 'data'
for (trt in trait){
  trt1<-paste(trt,".1rz", sep="")
  trt2<-paste(trt,".2rz", sep="")
  cortrt<-cor(data[,c(trt1, trt2)])
  pca<-eigen(cortrt)
  vector<-as.data.frame(pca[["vectors"]])
  vector<-as.matrix(vector)
  rawdata<-data[, c(trt1, trt2)]
  rawdata<-as.matrix(rawdata)
  name<-paste(trt, "pca", sep="")
  data[[name]]<-scale(c(rawdata%*%vector[,1]))
}

#Residualize temperament on age and gender and scale 
temp<-c("g6frus.1","g6fear.1", "g6atten.1", "g6act.1", "g6surg.1")
for (trt in temp){
  formula=as.formula(paste(trt," ~ AVE.AGE+GENDER"))
  model<-lm(formula, data=data)
  namet<-paste(trt, "rz", sep="")
  data[[namet]]<-scale(model$residuals)
}

#Super PCA, calculate and load into datafile, 'data'
rawsuperpca<-data[,c(148:186)]
covsuper<-cor(rawsuperpca)
superpca<-eigen(covsuper)
vectorsuper<-as.data.frame(superpca[["vectors"]])
vectorsuper<-as.matrix(vectorsuper)
rawsuperpca<-as.matrix(rawsuperpca)
data$SuperPCA1<-scale(c(rawsuperpca%*%vectorsuper[,1]))
data$SuperPCA2<-scale(c(rawsuperpca%*%vectorsuper[,2]))
data$SuperPCA3<-scale(c(rawsuperpca%*%vectorsuper[,3]))
data$SuperPCA4<-scale(c(rawsuperpca%*%vectorsuper[,4]))
data$SuperPCA5<-scale(c(rawsuperpca%*%vectorsuper[,5]))
data$SuperPCA6<-scale(c(rawsuperpca%*%vectorsuper[,6]))
data$SuperPCA7<-scale(c(rawsuperpca%*%vectorsuper[,7]))

#Code SuperPCA in direction of more risk for CP, DEP, or TEMP
data$SuperPCA1<-data$SuperPCA1*-1
data$SuperPCA2<-data$SuperPCA2*-1



#Predicting alcohol use
#ALC PROB
sub<-read_sav("C:/Users/wangfl/Box Sync/PCH analyses/Merge with w9 sub data.sav")
alc<-merge(data, sub, by="ISN")
alc$agealcC<-scale(alc$agealcon)
summary(modelzip<-glmmTMB(alcon ~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+
                            agealcC+GENDER.x+(1|familyid), 
                          zi=~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+
                            agealcC+GENDER.x+(1|familyid),family=poisson, data=alc))

#HEAVY DRINK
alc$hvydrkOLD2<-.bincode(alc$hvydrkOLD,breaks = c(-1,.5,100))-1

summary(modelHD2<-glmer(hvydrkOLD2 ~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+
                          agealcC+GENDER.x+(1|familyid) 
                         ,family=binomial, data=alc,  control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e5))))

summary(modelzipHD<-glmmTMB(hvydrkOLD ~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+
                              agealcC+GENDER.x+(1|familyid), 
                            zi=~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+
                              agealcC+GENDER.x+(1|familyid),family=poisson, data=alc))

tab_model(modelzip, modelzipHD,transform=NULL)

#AUC/ROC
library(ROCR)
alc$alconB<-.bincode(alc$alcon,breaks = c(-1,.5,100))-1
alc$hvydrkOLDB<-.bincode(alc$hvydrkOLD,breaks = c(-1,.5,100))-1


modelAP<-glmer(alconB~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+agealcon+GENDER.x+(1|familyid),
               family="binomial", data=alc, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e5)))
modelHD2<-glmer(hvydrkOLDB~SuperPCA1+SuperPCA2+SuperPCA3+SuperPCA4+SuperPCA5+SuperPCA6+SuperPCA7+agehvydrkOLD+GENDER.x+(1|familyid),
               family="binomial", data=alc, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun=2e5)))

alcna<-na.omit(alc)
alcna$predprobAP<-data.frame(predict(modelzip, newdata=alcna, type="response"))
pred_ROCR <- prediction(alcna$predprobAP, alcna$alcon)
roc_ROCR <- performance(pred_ROCR, measure = "tpr", x.measure = "fpr")
plot(roc_ROCR, main = "ROC curve", colorize = T)
abline(a = 0, b = 1)
auc_ROCR <- performance(pred_ROCR, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]

alcna$predprobHD<-data.frame(predict(modelHD, newdata=alcna, type="response"))
pred_ROCRHD <- prediction(alcna$predprobHD, alcna$g6SUUAUSE31B)
roc_ROCRHD <- performance(pred_ROCRHD, measure = "tpr", x.measure = "fpr")
plot(roc_ROCRHD, main = "ROC curve", colorize = T)
abline(a = 0, b = 1)
auc_ROCRHD <- performance(pred_ROCRHD, measure = "auc")
auc_ROCRHD <- auc_ROCRHD@y.values[[1]]


