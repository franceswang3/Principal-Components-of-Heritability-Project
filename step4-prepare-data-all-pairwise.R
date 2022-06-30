rm(list=ls()); gc()
options(sstringsAsFactors = F)

library(data.table)
library(readxl)

# this sets up the data for step in which VCE are calcualted for all pairwise combinations
# this include calculate PCH for the conduct problem and depression items

# specify if you want to analyze the standardized data (T) or not (F)
standardized=T

# read the data files
items=fread("Data/items_07-27-2021.txt",header=T,data.table=F)
ped=fread("Data/pedigree_07-27-2021.txt",header=T,data.table=F)
if(standardized){
  phn=fread("Data/std_phenotypes_07-27-2021.txt",header=T,data.table=F)
}else{
  phn=fread("Data/phenotypes_07-27-2021.txt",header=T,data.table=F)
}

# read the items that are combined in T1 and T2 using PCH
if(standardized){
  workdir="/data3/Bert/Frances/July2021/T1-T2-VCE-STD"
}else{
  workdir="/data3/Bert/Frances/July2021/T1-T2-VCE"
}

t12.items=read.table(paste(workdir,"/item_combinations.txt",sep=""),header=T)

# read the solutions to the covariates
SOLN=as.data.frame(read_excel(paste(workdir,"/VCE-SOLN-PCH-T1-T2.xlsx",sep=""),sheet = "SOLN"))
rownames(SOLN)=paste(SOLN$T,SOLN$COV,sep="_")
PCH=as.data.frame(read_excel(paste(workdir,"/VCE-SOLN-PCH-T1-T2.xlsx",sep=""),sheet = "PCH"))
rownames(PCH)=PCH$label

# adjust the phenotypes for the T1, T2 covariates
for(item in colnames(phn)[colnames(phn) %in% unlist(t12.items)]){
  if(substr(item,nchar(item),nchar(item)) == 1){
    phn[,item]=phn[,item]-as.matrix(phn[,c("INT","GENDER","AGE.1","AGE.2","AVE.AGE")])%*%SOLN[SOLN$T == "T1",paste(unlist(strsplit(item,"\\."))[1],"B",sep=".")]
  }else{
    phn[,item]=phn[,item]-as.matrix(phn[,c("INT","GENDER","AGE.1","AGE.2","AVE.AGE")])%*%SOLN[SOLN$T == "T2",paste(unlist(strsplit(item,"\\."))[1],"B",sep=".")]
  }
}

# perform PCH on the T1-T2 items
for(item in rownames(PCH)){
  phn[,paste(item,".pch",sep="")]=as.vector(t(as.matrix(PCH[item,c("item1","item2")])%*%t(as.matrix(phn[,c(paste(item,".1",sep=""),paste(item,".2",sep=""))]))))
  if(standardized){
    phn[,paste(item,".pch",sep="")]=scale(phn[,paste(item,".pch",sep="")])
  }
  phn[,paste(item,".pch",sep="")]=round(phn[,paste(item,".pch",sep="")],6)
}

# remove the columns that are for the items with pch
i.keep=which(!(colnames(phn) %in% unlist(t12.items)))
phn=phn[,i.keep]

# write the files
if(standardized){
  outdir="/data3/Bert/Frances/July2021/ALL-PAIRS-VCE-STD"
}else{
  outdir="/data3/Bert/Frances/July2021/ALL-PAIRS-VCE"
}
dir.create(outdir,showWarnings = F)

# now write the phenotype file
if(standardized){
  fwrite(phn,paste(outdir,"/std_phenotypes_pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
}else{
  fwrite(phn,paste(outdir,"/phenotypes_pch.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t",na=NA)
}




