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
  workdir="/data3/Bert/Frances/July2021/PCH-VCE-STD"
}else{
  workdir="/data3/Bert/Frances/July2021/PCH-VCE"
}


# read the item.combn
item.combinations=read.table(paste(workdir,"/item_combinations.txt",sep=""),header=T)



# collect the results from the multiple chains
ESTIMATES=NULL
n.chain=32
for(i.chain in 1:n.chain){
  estimates=fread(paste(workdir,"/MUV.chain",i.chain,"/postgibbs_samples",sep=""),header=F,data.table=F) # read the postgibbs samples
  estimates=estimates[,-c(1:3)] # remove useless columns
  ESTIMATES=rbind(ESTIMATES,estimates)
}
idx=which(lower.tri(matrix(0,ncol(item.combinations),ncol(item.combinations)),diag=T),arr.ind=T)
est.names=c(paste("g",idx[,2],idx[,1],sep="_"),paste("r",idx[,2],idx[,1],sep="_"))
colnames(ESTIMATES)=est.names

### calculate the correlations and heritabilities
for(i.pch in 1:ncol(item.combinations)){
  ESTIMATES[,paste("h2",i.pch,sep="_")]=ESTIMATES[,paste("g",i.pch,i.pch,sep="_")]/
    (ESTIMATES[,paste("g",i.pch,i.pch,sep="_")]+ESTIMATES[,paste("r",i.pch,i.pch,sep="_")])
}

idx.nd=which(lower.tri(matrix(0,ncol(item.combinations),ncol(item.combinations)),diag=F),arr.ind=T)
# genetic correlation
for(i in 1:nrow(idx.nd)){
  i.pch=idx.nd[i,2]; j.pch=idx.nd[i,1]
  ESTIMATES[,paste("rg",i.pch,j.pch,sep="_")]=ESTIMATES[,paste("g",i.pch,j.pch,sep="_")]/
    sqrt(ESTIMATES[,paste("g",i.pch,i.pch,sep="_")]*ESTIMATES[,paste("g",j.pch,j.pch,sep="_")])
}
# residual correlation
for(i in 1:nrow(idx.nd)){
  i.pch=idx.nd[i,2]; j.pch=idx.nd[i,1]
  ESTIMATES[,paste("rr",i.pch,j.pch,sep="_")]=ESTIMATES[,paste("r",i.pch,j.pch,sep="_")]/
    sqrt(ESTIMATES[,paste("r",i.pch,i.pch,sep="_")]*ESTIMATES[,paste("r",j.pch,j.pch,sep="_")])
}

# now for each pair of items determine the mean of the parameters and the highest posterior density
re.VCE=NULL
for(i in 1:nrow(idx.nd)){
  pch1=idx.nd[i,2]; pch2=idx.nd[i,1]
  re.VCE=rbind.data.frame(re.VCE,data.frame(trti=item.combinations[1,pch1],
                                      trtj=item.combinations[1,pch2],
                                      mean.gii=sprintf("%7.6f",mean(ESTIMATES[,paste("g",pch1,pch1,sep="_")])),
                                      mean.gij=sprintf("%7.6f",mean(ESTIMATES[,paste("g",pch1,pch2,sep="_")])),
                                      mean.gjj=sprintf("%7.6f",mean(ESTIMATES[,paste("g",pch2,pch2,sep="_")])),
                                      mean.rii=sprintf("%7.6f",mean(ESTIMATES[,paste("r",pch1,pch1,sep="_")])),
                                      mean.rij=sprintf("%7.6f",mean(ESTIMATES[,paste("r",pch1,pch2,sep="_")])),
                                      mean.rjj=sprintf("%7.6f",mean(ESTIMATES[,paste("r",pch2,pch2,sep="_")])),
                                      mean.h2i=sprintf("%7.6f",mean(ESTIMATES[,paste("h2",pch1,sep="_")])),
                                      mean.h2j=sprintf("%7.6f",mean(ESTIMATES[,paste("h2",pch2,sep="_")])),
                                      mean.rg=sprintf("%7.6f",mean(ESTIMATES[,paste("rg",pch1,pch2,sep="_")])),
                                      mean.rr=sprintf("%7.6f",mean(ESTIMATES[,paste("rr",pch1,pch2,sep="_")])),
                                      hpd.gii=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("g",pch1,pch1,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("g",pch1,pch1,sep="_")],0.05)["Upper Bound"]),
                                      hpd.gij=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("g",pch1,pch2,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("g",pch1,pch2,sep="_")],0.05)["Upper Bound"]),
                                      hpd.gjj=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("g",pch2,pch2,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("g",pch2,pch2,sep="_")],0.05)["Upper Bound"]),
                                      hpd.rii=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("r",pch1,pch1,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("r",pch1,pch1,sep="_")],0.05)["Upper Bound"]),
                                      hpd.rij=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("r",pch1,pch2,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("r",pch1,pch2,sep="_")],0.05)["Upper Bound"]),
                                      hpd.rjj=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("r",pch2,pch2,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("r",pch2,pch2,sep="_")],0.05)["Upper Bound"]),
                                      hpd.h2i=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("h2",pch1,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("h2",pch1,sep="_")],0.05)["Upper Bound"]),
                                      hpd.h2j=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("h2",pch2,sep="_")],0.05)["Lower Bound"],
                                                      boa.hpd(ESTIMATES[,paste("h2",pch2,sep="_")],0.05)["Upper Bound"]),
                                      hpd.rg=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("rg",pch1,pch2,sep="_")],0.05)["Lower Bound"],
                                                     boa.hpd(ESTIMATES[,paste("rg",pch1,pch2,sep="_")],0.05)["Upper Bound"]),
                                      hpd.rr=sprintf("%4.3f_%4.3f",boa.hpd(ESTIMATES[,paste("rr",pch1,pch2,sep="_")],0.05)["Lower Bound"],
                                                     boa.hpd(ESTIMATES[,paste("rr",pch1,pch2,sep="_")],0.05)["Upper Bound"])
  ))
}


# create matrices with variance components out of this
items=unlist(item.combinations)
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
re.SOLN=matrix(0,ncol=2*length(items),nrow=length(cov.labels))
colnames(re.SOLN)=paste(expand.grid(c("B","SE"),items)[,2],expand.grid(c("B","SE"),items)[,1],sep=".")
rownames(re.SOLN)=cov.labels

for(i.chain in 1:n.chain){
  soln=fread(paste(workdir,"/MUV.chain",i.chain,"/final_solutions",sep=""),header=F,data.table=F)[,-3] # don't need the third column
  colnames(soln)=c("ITEM","COV",paste("B",sep=""),paste("SE",sep=""))
  soln$ITEM=items[soln$ITEM]
  soln$COV=cov.labels[soln$COV]
  # collect the solutiopns
  re.SOLN[cbind(soln[,"COV"],paste(soln[,"ITEM"],"B",sep="."))]=re.SOLN[cbind(soln[,"COV"],paste(soln[,"ITEM"],"B",sep="."))]+soln[,"B"]
  re.SOLN[cbind(soln[,"COV"],paste(soln[,"ITEM"],"SE",sep="."))]=re.SOLN[cbind(soln[,"COV"],paste(soln[,"ITEM"],"SE",sep="."))]+soln[,"SE"]
}
# take the average for the B and SE
re.SOLN=re.SOLN/n.chain
re.SOLN=cbind.data.frame(T=NA,COV=cov.labels,re.SOLN)
  

# calculate the principal components
pchWeights<-function(vg,vr){
  # calculate the weights for the first PCH
  L=t(chol(vr))                    # determine the cholesky of the residual variance
  LGL=solve(L)%*%vg%*%t(solve(L))  # transform the genetic variance to account for the decomposition of the residual variance
  P=eigen(LGL)$vectors             # eigenvectors
  D=eigen(LGL)$values              # eigenvalues
  
  TM=t(P)%*%solve(L)              # transformation matrix
  rownames(TM)=paste("PCH",1:nrow(TM),sep=".")
  colnames(TM)=colnames(vr)
  
  out=list(PCH=t(TM),exp.h2=D/(1+D))
  return(out)
}
re.PCH=pchWeights(vg=G,vr=R)
PCH=re.PCH$PCH
exp.h2.PCH=data.frame(exp.h2=re.PCH$exp.h2)
rownames(exp.h2.PCH)=colnames(PCH)

# do some rounding
re.SOLN[,-c(1,2)]=round(re.SOLN[,-c(1,2)],6)
G=round(G,6)
R=round(R,6)
Cg=round(Cg,6)
Cr=round(Cr,6)
rownames(H2)=H2$item
H2$item=NULL
H2=round(H2,6)
PCH=round(PCH,6)
exp.h2.PCH=round(exp.h2.PCH,6)

# save all the results to a excel spreadsheet
OUT=list(VCE=re.VCE,SOLN=re.SOLN,G=as.data.frame(G),R=as.data.frame(R),Cg=as.data.frame(Cg),Cr=as.data.frame(Cr),H2=H2,PCH=as.data.frame(PCH),exp.h2.PCH=as.data.frame(exp.h2.PCH))
WriteXLS(OUT,paste(workdir,"/VCE-SOLN-MUV-PCH.xlsx",sep=""),row.names = T)



