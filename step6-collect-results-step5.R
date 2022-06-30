rm(list=ls()); gc()
options(stringsAsFactors = F)

# collect the results from step 2

library(data.table)
library(doParallel)
library(boa)
library(WriteXLS)

# define whether working on standardized results
standardized=T

# specify the work directory
if(standardized){
  workdir="/data3/Bert/Frances/July2021/ALL-PAIRS-VCE-STD"
}else{
  workdir="/data3/Bert/Frances/July2021/ALL-PAIRS-VCE"
}


# read the item.combn
item.combinations=read.table(paste(workdir,"/item_combinations.txt",sep=""),header=T)


#### collect the bivariate results
collectBivariateVCE<-function(label,item.combinations=item.combinations,workdir=workdir){
  estimates=fread(paste(workdir,"/",label,"/postgibbs_samples",sep=""),header=F,data.table=F) # read the postgibbs samples
  estimates=estimates[,-c(1:3)] # remove useless columns
  colnames(estimates)=c("gii","gij","gjj","rii","rij","rjj")  # give column names

  # determine derived parameters
  estimates[,"h2i"]=round(estimates[,"gii"]/rowSums(estimates[,c("gii","rii")]),6)
  estimates[,"h2j"]=round(estimates[,"gjj"]/rowSums(estimates[,c("gjj","rjj")]),6)
  estimates[,"rg"]=round(estimates[,"gij"]/sqrt(estimates[,c("gii")]*estimates[,c("gjj")]),6)
  estimates[,"rr"]=round(estimates[,"rij"]/sqrt(estimates[,c("rii")]*estimates[,c("rjj")]),6)
  
  # calculate the average for each of the estimates
  parameters=colnames(estimates)
  MN=NULL
  for(parameter in parameters){
    MN[paste("mean",parameter,sep=".")]=sprintf("%7.6f",mean(estimates[,parameter]))
  } 

  # and the highest posterior density
  HPD=NULL
  for(parameter in parameters){
    HPD[paste("hpd",parameter,sep=".")]=sprintf("%4.3f_%4.3f",boa.hpd(estimates[,parameter],0.05)["Lower Bound"],boa.hpd(estimates[,parameter],0.05)["Upper Bound"])
  } 
  
  # set up a data frame with the reults
  return(data.frame(label=label,trti=item.combinations[label,"item1"],trtj=item.combinations[label,"item2"],t(MN),t(HPD)))
}

re.VCE=as.data.frame(rbindlist(lapply(rownames(item.combinations),collectBivariateVCE,item.combinations=item.combinations,workdir=workdir)))
rownames(re.VCE)=re.VCE$label

# create matrices with variance components out of this
items=unique(c(re.VCE$trti,re.VCE$trtj))
# genetic
Cg=diag(length(items)); rownames(Cg)=items; colnames(Cg)=items
Cg[cbind(re.VCE$trti,re.VCE$trtj)]=as.numeric(re.VCE$mean.rg)
Cg[cbind(re.VCE$trtj,re.VCE$trti)]=as.numeric(re.VCE$mean.rg)
# residual
Cr=diag(length(items)); rownames(Cr)=items; colnames(Cr)=items
Cr[cbind(re.VCE$trti,re.VCE$trtj)]=as.numeric(re.VCE$mean.rr)
Cr[cbind(re.VCE$trtj,re.VCE$trti)]=as.numeric(re.VCE$mean.rr)

# collect the variances
Vg=rep(0,length(items)); names(Vg)=items
Vr=rep(0,length(items)); names(Vr)=items
for(item in items){
  Vg[item]=(sum(as.numeric(re.VCE[re.VCE$trti == item,"mean.gii"]))+sum(as.numeric(re.VCE[re.VCE$trtj == item,"mean.gjj"])))/(length(items)-1)
  Vr[item]=(sum(as.numeric(re.VCE[re.VCE$trti == item,"mean.rii"]))+sum(as.numeric(re.VCE[re.VCE$trtj == item,"mean.rjj"])))/(length(items)-1)
}

H2=data.frame(item=items,h2=Vg/(Vg+Vr))

# determine the covariance matrices, calculated fromt he average variance estimate and the estimate of the correlation
G=diag(sqrt(Vg))%*%Cg%*%diag(sqrt(Vg)); rownames(G)=items; colnames(G)=items
R=diag(sqrt(Vr))%*%Cr%*%diag(sqrt(Vr)); rownames(R)=items; colnames(R)=items



# collect the covariate estimates
cov.labels=c("INT","ISN","GENDER","AGE.1","AGE.2","AVE.AGE")

re.SOLN=NULL
re.SOLN=data.frame(T=NA,COV=cov.labels,matrix(0,ncol=2*length(items),nrow=length(cov.labels)))
colnames(re.SOLN)[-c(1,2)]=paste(expand.grid(c("B","SE"),items)[,2],expand.grid(c("B","SE"),items)[,1],sep=".")
rownames(re.SOLN)=cov.labels
for(label in rownames(item.combinations)){
  soln=fread(paste(workdir,"/",label,"/final_solutions",sep=""),header=F,data.table=F)[,-3] # don't need the third column
  colnames(soln)=c("T","COV",paste("B",sep=""),paste("SE",sep=""))
  soln$COV=cov.labels[soln$COV]  
  # for item1
  item1=item.combinations[label,"item1"]
  re.SOLN[soln[soln[,"T"] == 1,"COV"],paste(item1,"B",sep=".")]=re.SOLN[soln[soln[,"T"] == 1,"COV"],paste(item1,"B",sep=".")]+soln[soln[,"T"] == 1,"B"] # B
  re.SOLN[soln[soln[,"T"] == 1,"COV"],paste(item1,"SE",sep=".")]=re.SOLN[soln[soln[,"T"] == 1,"COV"],paste(item1,"SE",sep=".")]+soln[soln[,"T"] == 1,"SE"] # SE
  # for item1
  item2=item.combinations[label,"item2"]
  re.SOLN[soln[soln[,"T"] == 2,"COV"],paste(item2,"B",sep=".")]=re.SOLN[soln[soln[,"T"] == 2,"COV"],paste(item2,"B",sep=".")]+soln[soln[,"T"] == 2,"B"] # B
  re.SOLN[soln[soln[,"T"] == 2,"COV"],paste(item2,"SE",sep=".")]=re.SOLN[soln[soln[,"T"] == 2,"COV"],paste(item2,"SE",sep=".")]+soln[soln[,"T"] == 2,"SE"] # SE
  
}
# take the average for the B and SE
re.SOLN[,-c(1,2)]=re.SOLN[,-c(1,2)]/(length(items)-1)
  

# do some rounding
re.SOLN[,-c(1,2)]=round(re.SOLN[,-c(1,2)],6)
G=round(G,6)
R=round(R,6)
Cg=round(Cg,6)
Cr=round(Cr,6)
H2$item=NULL
H2=round(H2,6)

# save all the results to a excel spreadsheet
OUT=list(VCE=re.VCE,SOLN=re.SOLN,G=as.data.frame(G),R=as.data.frame(R),Cg=as.data.frame(Cg),Cr=as.data.frame(Cr),H2=H2)
WriteXLS(OUT,paste(workdir,"/VCE-SOLN-ALL-PAIRS.xlsx",sep=""),row.names = T)



