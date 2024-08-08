rm(list=ls()); gc()
options(stringAsFactors = F)

tissue="SYNAPTOSOME-PEPTIDE"

########################################################
#### THIS VERSION SKIPS THE SURROGATE VARIABLES AND 
#### TAKES ITS INPUT FROM THE MODEL SELECTION STEP
########################################################

# this part performs differential intensity analysis using a paired t-test
# the paaorch is to ranomly shuffle the ASD/NT status in each of n.perm permutations
# this will give a distribution for the test statistic and this can then be used to calculate 
# a p-value



#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(limma)
require(sva)
require(qvalue)

########################
# FUNCTIONS
########################
permutationFunction<-function(i,df){
  # this function performs each of the permutations
  # permute within pairs
  iperm=sample(2,nrow(df),replace=T) # random determine whether cases and controls should be swapped or not
  fwd.perm=c("ctl.int","case.int")[iperm]  # used to select which sample to be treated as case and which as control
  rev.perm=c("case.int","ctl.int")[iperm]  # this selecs the other sample
  # determine the paired-t-test statistic
  t.stat=t.test(as.numeric(df[cbind(rownames(df),fwd.perm)]),as.numeric(df[cbind(rownames(df),rev.perm)]),paired=T)$statistic
  return(t.stat)
}

permutationTest <- function(feature,n.perm=10000,log2.int,df){
  # set up data by adding the intensity for the feature to the data frame
  df$ctl.int=log2.int[feature,df$ctl.id]
  df$case.int=log2.int[feature,df$case.id]
  #calcualte the observed data paired-t-test statistic
  re.ttest=t.test(df$case.int,df$ctl.int,paired=T)
  # collect useful information
  b=re.ttest$estimate
  deg.fr=re.ttest$parameter
  se=mean(abs(re.ttest$conf.int-b)/qt(.975,df=deg.fr))
  t.obs=re.ttest$statistic
  
  # perform permutations
  t.perm=unlist(lapply(1:n.perm,permutationFunction,df=df))
  
  
  # collect results mn and sd can be used for an approximate test statistic, this is currently not used
  df.out=data.frame(feature,effect=b,std.error=se,statistic=t.obs,p.emperical=sum(abs(t.perm) >= abs(t.obs))/n.perm,mn.perm=mean(t.perm),sd.perm=sd(t.perm))
  return(df.out)
}



### load the normalized data
load(paste0("WorkData/NORMALIZED-",tissue,".RData"))



# set up a paired data data frame with identification and diagnosis
df=data.frame(pair=unique(meta.samples$PAIR))
rownames(df)=df$pair   # 28 pairs
df$ctl.id=rownames(meta.samples)[match(paste(df$pair,"CTL",sep="-"),paste(meta.samples$PAIR,meta.samples$DX,sep="-"))]
df$case.id=rownames(meta.samples)[match(paste(df$pair,"ASD",sep="-"),paste(meta.samples$PAIR,meta.samples$DX,sep="-"))]

# run the base analysis for all peptides use 1e3 permutations
system.time(re.perm.all<-rbindlist(mclapply(rownames(norm.pep),permutationTest,n.perm=1e3,log2.int = log2(norm.pep),df=df,mc.cores=min(190,nrow(norm.pep)),mc.preschedule = T,mc.silent=T)))
re.perm.all=as.data.frame(re.perm.all); rownames(re.perm.all)=re.perm.all$feature
re.perm.all$nperm="1e3"

# select proteins with p < 0.05 use 1e4 permutations
peptides=re.perm.all$feature[re.perm.all$p.emperical < 0.05] # 3370
system.time(re.perm.05<-rbindlist(mclapply(peptides,permutationTest,n.perm=1e4,log2.int = log2(norm.pep[peptides,]),df=df,mc.cores=min(190,length(peptides)),mc.preschedule = T,mc.silent=T)))
re.perm.05=as.data.frame(re.perm.05); rownames(re.perm.05)=re.perm.05$feature
re.perm.05$nperm="1e4"

# select proteins with p < 0.005 use 1e5 permutations
peptides=re.perm.05$feature[re.perm.05$p.emperical < 0.005] # 538
system.time(re.perm.005<-rbindlist(mclapply(peptides,permutationTest,n.perm=1e5,log2.int = log2(norm.pep[peptides,]),df=df,mc.cores=min(190,length(peptides)),mc.preschedule = T,mc.silent=T)))
re.perm.005=as.data.frame(re.perm.005); rownames(re.perm.005)=re.perm.005$feature
re.perm.005$nperm="1e5"

# select proteins with p < 0.0005 use 1e6 permutations
peptides=re.perm.005$feature[re.perm.005$p.emperical < 0.0005] # 87
system.time(re.perm.0005<-rbindlist(mclapply(peptides,permutationTest,n.perm=1e6,log2.int = log2(norm.pep[peptides,]),df=df,mc.cores=min(190,length(peptides)),mc.preschedule = T,mc.silent=T)))
re.perm.0005=as.data.frame(re.perm.0005); rownames(re.perm.0005)=re.perm.0005$feature
re.perm.0005$nperm="1e6"

re.perm.all[rownames(re.perm.05),]=re.perm.05
re.perm.all[rownames(re.perm.005),]=re.perm.005
re.perm.all[rownames(re.perm.0005),]=re.perm.0005

pdf(paste0("Plots/emperical-p-values-peptides-",tissue,".pdf"),height=8,width=12)
hist(re.perm.all$p.emperical,breaks=seq(0,1,.01),las=1,xlab="emperical P",main="",col=viridis(3)[2])
dev.off()


colnames(re.perm.all)[1]="peptide"

1-qvalue(re.perm.all$p.emperical)$pi0 # 0.225
re.perm.all$q=qvalue(re.perm.all$p.emperical)$qvalues
sum(re.perm.all$q < 0.05) # 2
sum(re.perm.all$q < 0.10) # 69
sum(re.perm.all$q < 0.20) # 792
max(re.perm.all$p.emperical[re.perm.all$q < 0.05])  # 0.000003
max(re.perm.all$p.emperical[re.perm.all$q < 0.10])  # 0.00035
max(re.perm.all$p.emperical[re.perm.all$q < 0.20])  # 0.0080


save(re.perm.all,re.perm.05,re.perm.005,re.perm.0005,file=paste0("WorkData/PAIRED-T-TEST-",tissue,".RData"))


# finish up and create a clean file.
re.paired.t.test=re.perm.all
save(re.paired.t.test,re.perm.all,re.perm.05,re.perm.005,re.perm.0005,file=paste0("WorkData/PAIRED-T-TEST-",tissue,".RData"))

WriteXLS::WriteXLS(re.paired.t.test,paste0("Results/PAIRED-T-TEST-",tissue,".xlsx"),SheetNames = tissue)


