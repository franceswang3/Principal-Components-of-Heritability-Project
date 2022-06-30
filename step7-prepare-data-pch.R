rm(list=ls()); gc()
options(sstringsAsFactors = F)

library(data.table)
library(readxl)

# this sets up the data for step in which VCE are calculated for all the pairwise pch
# first the input items are adjusted for the results from the previous step
# next the items are clustered based on the genetic correlation matrix
# finally the first for pch for each cluster is determined
# this include calculate PCH for the conduct problem and depression items

# specify if you want to analyze the standardized data (T) or not (F)
standardized=T

# read the data files
ped=fread("Data/pedigree_07-27-2021.txt",header=T,data.table=F)
if(standardized){
  phn=fread("ALL-PAIRS-VCE-STD/std_phenotypes_pch.txt",header=T,data.table=F)
  SOLN=as.data.frame(read_excel("ALL-PAIRS-VCE-STD/VCE-SOLN-ALL-PAIRS.xlsx",sheet = "SOLN"))
  G=as.data.frame(read_excel("ALL-PAIRS-VCE-STD/VCE-SOLN-ALL-PAIRS.xlsx",sheet = "G"))
  R=as.data.frame(read_excel("ALL-PAIRS-VCE-STD/VCE-SOLN-ALL-PAIRS.xlsx",sheet = "R"))
}else{
  phn=fread("ALL-PAIRS-VCE/phenotypes_pch.txt",header=T,data.table=F)
  SOLN=as.data.frame(read_excel("ALL-PAIRS-VCE/VCE-SOLN-ALL-PAIRS.xlsx",sheet = "SOLN"))
  G=as.data.frame(read_excel("ALL-PAIRS-VCE/VCE-SOLN-ALL-PAIRS.xlsx",sheet = "G"))
  R=as.data.frame(read_excel("ALL-PAIRS-VCE/VCE-SOLN-ALL-PAIRS.xlsx",sheet = "R"))
}
rownames(SOLN)=SOLN$...1
rownames(G)=G$...1; G=as.matrix(G[,-1])
rownames(R)=R$...1; R=as.matrix(R[,-1])

### adjust the items for the solutions of the covariates
for(item in colnames(phn)[!(colnames(phn) %in% rownames(SOLN))]){
  phn[,item]=phn[,item]-as.matrix(phn[,c("INT","GENDER","AGE.1","AGE.2","AVE.AGE")])%*%SOLN[c("INT","GENDER","AGE.1","AGE.2","AVE.AGE"),paste(item,"B",sep=".")]
} 

res=phn
### only needed when using them to find the important items for each of the SUPER-PCH
if(standardized){
  for(item in colnames(res)[-c(1:6)]){
    res[,item]=round(scale(res[,item]),6)
  }
  fwrite(res,paste("ALL-PAIRS-VCE-STD/std_residual_phenotypes.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
}else{
  fwrite(res,paste("ALL-PAIRS-VCE/residual_pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
}


### determine the clustering based on the genetic correlation matrix
H=hclust(as.dist(1-round(cov2cor(G),6)))
pdf("item-clusters-h-0.70.pdf",height=6,width=8)
plot(H)
hcut=0.70
rect.hclust(H,h=hcut)
dev.off()
Htree=cutree(H,h=hcut)  # find the items clustering together
# collect the clusters
clusters=NULL
for(i in sort(unique(Htree))){
  clusters[[paste("pch",i,sep="")]]=names(Htree)[Htree == i]
}
cluster.list=data.frame(item=names(Htree),cluster=Htree)
cluster.list=cluster.list[order(cluster.list$cluster),]

#### calculate the pch

primaryPchWeights<-function(vg,vr){
  # calculate the weights for the first PCH
  L=t(chol(vr))                    # determine the cholesky of the residual variance
  LGL=solve(L)%*%vg%*%t(solve(L))  # transform the genetic variance to account for the decomposition of the residual variance
  P=eigen(LGL)$vectors             # eigenvectors
  D=eigen(LGL)$values              # eigenvalues
  
  TM=t(P)%*%solve(L)              # transformation matrix
  rownames(TM)=paste("PCH",1:nrow(TM),sep=".")
  colnames(TM)=rownames(vg)
  
  out=list(PCH=t(TM["PCH.1",]),exp.h2=D[1]/(1+D[1]))
  print(D);flush.console()
  return(out)
}

# calculate the pch weightd antransform the data
PCH.results=matrix(NA,nrow=length(clusters),ncol=ncol(G)+1)
rownames(PCH.results)=paste("pch",1:length(clusters),sep="")
colnames(PCH.results)=c(colnames(G),"exp.h2")
for(pch in names(clusters)){
  items=clusters[[pch]]
  if(length(items) > 1){
    re.pch=primaryPchWeights(vg=G[items,items],vr=R[items,items])
    pch.weights=re.pch$PCH
    PCH.results[pch,colnames(pch.weights)]=pch.weights
    PCH.results[pch,"exp.h2"]=re.pch$exp.h2
    phn[,pch]=as.matrix(phn[,items])%*%as.matrix(pch.weights[,items])
  }else{
    phn[,pch]=phn[,items]
    PCH.results[pch,items]=1
    PCH.results[pch,"exp.h2"]=G[items,items]/(G[items,items]+R[items,items])
  }
  if(standardized){
    phn[,pch]=scale(phn[,pch])
  }
  phn[,pch]=round(phn[,pch],6)
}
PCH.results=round(PCH.results,6)
# remove all previous items. only retain the pch
phn=phn[,!(colnames(phn) %in% colnames(G))]




# write the files
if(standardized){
  outdir="/data3/Bert/Frances/July2021/PCH-VCE-STD"
}else{
  outdir="/data3/Bert/Frances/July2021/PCH-VCE"
}
dir.create(outdir,showWarnings = F)

if(standardized){
  fwrite(phn,paste(outdir,"/std_pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
  fwrite(cluster.list,paste(outdir,"/std_pch_clusters.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
  WriteXLS::WriteXLS(as.data.frame(PCH.results),paste(outdir,"/std-pch-weight-exp-h2.xlsx",sep=""),SheetNames = "pch",row.names = T)
}else{
  fwrite(phn,paste(outdir,"/pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
  fwrite(cluster.list,paste(outdir,"/pch_clusters.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
  WriteXLS::WriteXLS(as.data.frame(PCH.results),paste(outdir,"/pch-weight-exp-h2.xlsx",sep=""),SheetNames = "pch",row.names = T)
}




