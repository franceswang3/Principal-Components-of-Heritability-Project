rm(list=ls()); gc()
options(stringsAsFactors = F)

# collect the results from step 11

library(data.table)
library(doParallel)
library(boa)
library(WriteXLS)

# define whether working on standardized results
standardized=T 

# specify the work directory
if(standardized){
  workdir="/data3/Bert/Frances/July2021/SUPER-PCH-VCE-STD/"
}else{
  workdir="/data3/Bert/Frances/July2021/PCH-VCE-STD"
}


# read the item.combn
item.combinations=read.table(paste(workdir,"/item_combinations.txt",sep=""),header=T)


#### collect the bivariate results
collectUnivariateVCE<-function(label,item.combinations=item.combinations,workdir=workdir){
  estimates=fread(paste(workdir,"/",label,"/postgibbs_samples",sep=""),header=F,data.table=F) # read the postgibbs samples
  estimates=estimates[,-c(1:3)] # remove useless columns
  colnames(estimates)=c("gii","rii")  # give column names

  # determine derived parameters
  estimates[,"h2i"]=round(estimates[,"gii"]/rowSums(estimates[,c("gii","rii")]),6)

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
  return(data.frame(label=label,trti=item.combinations[label,"item"],t(MN),t(HPD)))
}

re.VCE=as.data.frame(rbindlist(lapply(rownames(item.combinations),collectUnivariateVCE,item.combinations=item.combinations,workdir=workdir)))
rownames(re.VCE)=re.VCE$label


# collect the covariate estimates
cov.labels=c("INT","ISN","GENDER","AGE.1","AGE.2","AGE.AVE")

re.SOLN=NULL
for(label in rownames(item.combinations)){
 soln=fread(paste(workdir,"/",label,"/final_solutions",sep=""),header=F,data.table=F)[,-3] # don't need the third column
  colnames(soln)=c("T","COV",paste(label,".B",sep=""),paste(label,".SE",sep=""))
  soln$T=paste("T",soln$T,sep="")
  soln$COV=cov.labels[soln$COV]
  # round to 3 digits
  soln[,paste(label,".B",sep="")]=round(soln[,paste(label,".B",sep="")],6)
  soln[,paste(label,".SE",sep="")]=round(soln[,paste(label,".SE",sep="")],6)
  
  if(is.null(re.SOLN)){
    # set up the beginning of the data frame to collect all results
    re.SOLN=soln[,c("T","COV")]
  }
  # now add the remaining estimates
  re.SOLN=cbind.data.frame(re.SOLN,data.frame(soln[,c(paste(label,".B",sep=""),paste(label,".SE",sep=""))]))
}
  
# and calculating the primary PCH
# function to perform PCH and extract the primary PCH
primaryPchWeights<-function(vg,vr){
  # calculate the weights for the first PCH
  L=t(chol(vr))                    # determine the cholesky of the residual variance
  LGL=solve(L)%*%vg%*%t(solve(L))  # transform the genetic variance to account for the decomposition of the residual variance
  P=eigen(LGL)$vectors             # eigenvectors
  D=eigen(LGL)$values              # eigenvalues
  
  TM=t(P)%*%solve(L)              # transformation matrix
  rownames(TM)=paste("PCH",1:nrow(TM),sep=".")
  colnames(TM)=c("item1","item2")

  out=data.frame(t(TM["PCH.1",]),exp.h2=D[1]/(1+D[1]))
  return(out)
}


re.PCH=NULL
for(label in rownames(item.combinations)){
  G=matrix(NA,2,2)
  R=matrix(NA,2,2)
  
  # fill in the matrices
  G[1,1]=as.numeric(re.VCE[label,"mean.gii"]); G[1,2]=as.numeric(re.VCE[label,"mean.gij"]); G[2,1]=as.numeric(re.VCE[label,"mean.gij"]);G[2,2]=as.numeric(re.VCE[label,"mean.gjj"])
  rownames(G)=item.combinations[label,]; colnames(G)=item.combinations[label,]
  R[1,1]=as.numeric(re.VCE[label,"mean.rii"]); R[1,2]=as.numeric(re.VCE[label,"mean.rij"]); R[2,1]=as.numeric(re.VCE[label,"mean.rij"]);R[2,2]=as.numeric(re.VCE[label,"mean.rjj"])
  rownames(R)=item.combinations[label,]; colnames(R)=item.combinations[label,]
  
  # calculate the primary PCH
  re.PCH=rbind.data.frame(re.PCH,data.frame(label,primaryPchWeights(vg=G,vr=R)))
}
rownames(re.PCH)=rownames(item.combinations)
re.PCH[,-1]=round(re.PCH[,-1],6)


# save all the results to a excel spreadsheet
OUT=list(VCE=re.VCE,SOLN=re.SOLN)
WriteXLS(OUT,paste(workdir,"/VCE-SOLN-SUPER-PCH.xlsx",sep=""))



