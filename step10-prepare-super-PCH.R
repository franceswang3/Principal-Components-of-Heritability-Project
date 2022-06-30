rm(list=ls()); gc()
options(sstringsAsFactors = F)

library(data.table)
library(readxl)

# this sets up the data for super PCH
# first the input items are adjusted for the results from the previous step
# next the items are clustered based on the genetic correlation matrix
# finally the first for pch for each cluster is determined
# this include calculate PCH for the conduct problem and depression items

# specify if you want to analyze the standardized data (T) or not (F)
standardized=T

# read the data files
ped=fread("Data/pedigree_07-27-2021.txt",header=T,data.table=F)
if(standardized){
  phn=fread("PCH-VCE-STD/std_pch.txt",header=T,data.table=F)
  SOLN=as.data.frame(read_excel("PCH-VCE-STD/VCE-SOLN-MUV-PCH.xlsx",sheet = "SOLN"))
  G=as.data.frame(read_excel("PCH-VCE-STD/VCE-SOLN-MUV-PCH.xlsx",sheet = "G"))
  R=as.data.frame(read_excel("PCH-VCE-STD/VCE-SOLN-MUV-PCH.xlsx",sheet = "R"))
  PCH=as.data.frame(read_excel("PCH-VCE-STD/VCE-SOLN-MUV-PCH.xlsx",sheet = "PCH"))
}else{
  phn=fread("PCH-VCE/pch.txt",header=T,data.table=F)
  SOLN=as.data.frame(read_excel("PCH-VCE/VCE-SOLN-MUV-PCH.xlsx",sheet = "SOLN"))
  G=as.data.frame(read_excel("PCH-VCE/VCE-SOLN-MUV-PCH.xlsx",sheet = "G"))
  R=as.data.frame(read_excel("PCH-VCE/VCE-SOLN-MUV-PCH.xlsx",sheet = "R"))
  PCH=as.data.frame(read_excel("PCH-VCE/VCE-SOLN-MUV-PCH.xlsx",sheet = "PCH"))
}
rownames(SOLN)=SOLN$...1
rownames(G)=G$...1; G=as.matrix(G[,-1])
rownames(R)=R$...1; R=as.matrix(R[,-1])
rownames(PCH)=PCH$...1; PCH=as.matrix(PCH[,-1])

### adjust the items for the solutions of the covariates
for(item in colnames(phn)[!(colnames(phn) %in% rownames(SOLN))]){
  phn[,item]=phn[,item]-as.matrix(phn[,c("INT","GENDER","AGE.1","AGE.2","AVE.AGE")])%*%SOLN[c("INT","GENDER","AGE.1","AGE.2","AVE.AGE"),paste(item,"B",sep=".")]
} 


# calculate the PCH
for(pch in colnames(PCH)){
  phn[,pch]=as.matrix(phn[,rownames(PCH)])%*%PCH[,pch]
}

# remove all previous items. only retain the pch
phn=phn[,!(colnames(phn) %in% colnames(G))]

# scale the super PCH when standardized is chosen
if(standardized){
  for(pch in colnames(PCH)){
    phn[,pch]=scale(phn[,pch])
  }
}

# round the pch
for(pch in colnames(PCH)){
  phn[,pch]=round(phn[,pch],6)
}


# write the files
if(standardized){
  outdir="/data3/Bert/Frances/July2021/SUPER-PCH-VCE-STD"
}else{
  outdir="/data3/Bert/Frances/July2021/SUPER-PCH-VCE"
}
dir.create(outdir,showWarnings = F)

if(standardized){
  fwrite(phn,paste(outdir,"/std_super_pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
}else{
  fwrite(phn,paste(outdir,"/super_pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
}



