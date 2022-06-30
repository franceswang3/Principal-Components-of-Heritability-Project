rm(list=ls()); gc()
options(stringsAsFactors = F)

# run GIBBS sampler for the items at the two time points
# you will get the following error
# *** error on final_solutions
# this can be ignored

library(data.table)
library(doParallel)

setwd("/data3/Bert/Frances/July2021")

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

# since this is estimate the VC between T1 and T2 for each item we do not use the temperament items
items=items[items$category %in% c("conduct problem","depression"),]
# set up a data.frame with the combination to analyze
item.combn=data.frame(item1=paste(items$item,".1",sep=""),item2=paste(items$item,".2",sep=""))
rownames(item.combn)=items$item

# specify the covariance model to use
models=list()
models[["AGE.1"]]=c("INT","ISN","GENDER","AGE.1")
models[["AGE.2"]]=c("INT","ISN","GENDER","AGE.2")

# create the the working directory for this analysis
if(standardized){
  workdir="/data3/Bert/Frances/July2021/T1-T2-VCE-STD"
}else{
  workdir="/data3/Bert/Frances/July2021/T1-T2-VCE"
}
dir.create(workdir,showWarnings = F,recursive = T)
setwd(workdir) # go to the directory

fwrite(item.combn,"item_combinations.txt",row.names=T,col.names=T,quote=F,sep="\t",na=NA)

runGibbsSampler<-function(i.combn,item.combn=item.combn,models=models,ped=ped,phn=phn,workdir=workdir,n.gibbs,n.burnin,n.thin,n.backup){

  setwd(workdir)
  # set up the data
  items=unlist(item.combn[i.combn,])
  # create a directory in which the estimation occurs
  estm.dir=paste(workdir,"/",rownames(item.combn)[i.combn],sep="")
  dir.create(estm.dir,recursive = T,showWarnings = F)
  # go to this directory
  setwd(estm.dir)
  
  # determine which model to use for each item in the combination
  item.model=list()
  for(item in items){
    if(substr(item,nchar(item),nchar(item)) == "1"){
      # model for AGE.1
      item.model[[item]]="AGE.1"
    }else{
      # model for AGE.2
      item.model[[item]]="AGE.2"
    }
  }
  
  # set up the pedigree and phenotype files
  # pedigree
  fwrite(ped,"ped",row.names=F,col.names=F,quote=F,sep=" ")
  
  # collect information on model, it is easiest to carry all covariates that you have (even the ones that you might not use)
  effects=c("INT","ISN","GENDER","AGE.1","AGE.2","AVE.AGE")
  for(item in items){
    effects=unique(c(effects,models[[item.model[[item]]]]))
  }
  # add to that the ites of interest
  phn.columns=c(effects,items)
  sel.phn=phn[,phn.columns]
  if(length(items) > 1){
    sel.phn=sel.phn[rowSums(sel.phn[,items]) != 0,]
  }else{
    sel.phn=sel.phn[sel.phn[,items] != 0,]
  }
  fwrite(sel.phn,"phn",row.names=F,col.names=F,quote=F,sep=" ")

  # set up a model data.frame, this is input for the gibbs sampler
  model.df=data.frame(effect=effects)
  for(item in items){
    model.df[match(models[[item.model[[item]]]],model.df$effect),paste("effect.col.",item,sep="")]=match(models[[item.model[[item]]]],colnames(sel.phn),nomatch=0)
  }
  model.df[is.na(model.df)]=0
  model.df=cbind.data.frame(model.df,data.frame(effect.N=1,effect.type="cov"))
  model.df$effect.type[model.df$effect %in% "ISN"]="cross"
  model.df$effect.N[model.df$effect == "ISN"]=max(ped$ISN)

  # starting values of variance components
  R=matrix(.5,length(items),length(items)); diag(R)=1
  G=matrix(.5,length(items),length(items)); diag(G)=1
  

  #### set up the gibbs sampler input file
  gibbs.inp="gibbs-sampler.inp"
  write.table("DATAFILE",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=F)
  write.table("phn",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("NUMBER_OF_TRAITS",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(length(items),gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("NUMBER_OF_EFFECTS",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(length(effects),gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("OBSERVATION(S)",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(paste(match(items,colnames(sel.phn)),collapse=" "),gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("WEIGHT(S)",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT [EFFECT NESTED]",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(model.df[,-1],gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  
  ### random residual
  write.table("RANDOM_RESIDUAL",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(R,gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  
  ### random additive effect
  write.table("RANDOM_GROUP",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(which(model.df$effect == "ISN"),gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("RANDOM_TYPE",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("add_animal",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("FILE",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("ped",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table("(CO)VARIANCES",gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(G,gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)

  #### add OPTIONS
  write.table(paste("OPTION","seed",paste(t(sample(10000000,2)),collapse=" ")),gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(paste("OPTION save_halfway_samples",n.backup,sep=" "),gibbs.inp,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(paste("OPTION solution mean",paste(which(model.df$effect != "ISN"),collapse = " "),sep=" "),gibbs.inp,row.names=F,col.names=F,quote=F,append=T)
  write.table(paste("OPTION missing -999",paste(which(model.df$effect != "ISN"),collapse = " "),sep=" "),gibbs.inp,row.names=F,col.names=F,quote=F,append=T)
  
  ### now create the gibbs sampler control file
  gibbs.control="gibbs-sampler.cntrl"
  write.table(gibbs.inp,gibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=F)
  write.table(paste(n.gibbs,n.burnin,sep=" "),gibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(n.thin,gibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  
  # gibbs sampler program
  gibbs.program="/data3/Software/BLUPF90/gibbs3f90"

  gibbs.command=paste(gibbs.program," < ",gibbs.control," > ","gibbs.log")
  system(gibbs.command)  

  ### run postgibbs
  postgibbs.control="postgibbs.cntrl"
  write.table(gibbs.inp,postgibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=F)
  write.table(0,postgibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(n.thin,postgibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  write.table(0,postgibbs.control,row.names=F,col.names=F,quote=F,sep=" ",append=T)
  
  postgibbs.program="/data3/Software/BLUPF90/postgibbsf90"
  pg.command=paste(postgibbs.program," < ",postgibbs.control," > ","postgibbs.log")
  system(pg.command)

  
  setwd(workdir)
  return(items)
}

system.time(re.gibbs<-mclapply(seq(1,nrow(item.combn),1),runGibbsSampler,item.combn=item.combn,
                             models=models,ped=ped,phn=phn,workdir=workdir,
                             n.gibbs=3*550000,n.burnin=3*50000,n.thin=3*100,n.backup=10000,
                             mc.cores=min(nrow(item.combn),100),mc.preschedule = F,mc.silent = T))

setwd("/data/Bert/Frances/July2021")
