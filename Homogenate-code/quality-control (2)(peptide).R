rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOMOGENATE-PEPTIDE"
# this code performs basic QC steps. Since these are dataset dependent the user will need to make the appropriate changes to this code

#### libraries used
library(data.table)           # fast version of read.table
library(readxl)               # used for reading Excel format files
library(viridis)              # colorblind friendly color palette
library(matrixStats)          # to get statistics on matrices

load(paste0("WorkData/RAW-",tissue,".RData"))
dim(raw.pep)  # 79,952 x 80

# remove peptides without intensity data
sum(rowSums(!is.na(raw.pep)) == 0)  # 6652
raw.pep=raw.pep[rowSums(!is.na(raw.pep)) > 0,] 
dim(raw.pep)  # 73,300 x 74
meta.pep=meta.pep[rownames(raw.pep),]  # 73,300


# remove peptides with missing master protein (ROLLUP) information   - 121
sum(is.na(meta.pep$ROLLUP)) # - 185
meta.pep=meta.pep[!is.na(meta.pep$ROLLUP),]  
dim(meta.pep)                # 73,115


length(unique(unlist(strsplit(meta.pep$ROLLUP,";")))) # 7514
# select the corresponding peptide intensities
raw.pep=raw.pep[rownames(meta.pep),]  
dim(raw.pep)                     # 73,115 x 80


###### also retain peptides with non-missing values for half the samples and half the pools
raw.pep=raw.pep[rowSums(is.na(raw.pep[,samples])) <= length(samples)/2 & 
                  rowSums(is.na(raw.pep[,pools])) <= length(pools)/2,]
dim(raw.pep)  # 35,103 x 80
meta.pep=meta.pep[rownames(raw.pep),]
length(unique(unlist(strsplit(meta.pep$ROLLUP,";")))) # 5746


###### THIS STEP IS UP FOR DEBATE #############################
# remove outliers
# determine outliers by plex and peptide
outlierDetection<-function(feature,XX,PLEX,range=1.5){
  tmp=XX[feature,]
  # standardize to a common mean
  mn=mean(XX[feature,],na.rm=T)
  mns=aggregate(XX[feature,]~factor(PLEX,levels=unique(PLEX)),FUN="mean") # mean by plex
  plex.means=unlist(mns[,2]); names(plex.means)=mns[,1]   
  
  tmp=tmp*(mn/plex.means[PLEX])              # adjusteds observation
  outliers=names(boxplot(log2(tmp),range=range,plot=F)$out)     # outliers based in log2 transformed data using range (quantile box +/- range*IQR)
  if(length(outliers) > 0){                                         # make a list of outliers
    re.out=data.frame(feature=feature,sample=outliers)
  }else{
    re.out=data.frame(feature=feature,sample=NA)
  }
  return(re.out)
}

outliers=rbindlist(mclapply(rownames(raw.pep),outlierDetection,XX=log2(raw.pep),PLEX=meta.samples$PLEX,range=3,mc.cores=64,mc.preschedule=T))
outliers=as.data.frame(outliers[!is.na(outliers$sample),])
nrow(outliers)
nrow(outliers)/sum(!is.na(raw.pep))
### set outliers to missing
qc.pep=raw.pep
# SKIP THE NEXT LINE WHEN YOU DO NOT WANT TO REMOVE OUTLIERS
qc.pep[cbind(outliers$feature,outliers$sample)]=NA


save(meta.samples,meta.pep,qc.pep,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/QC-",tissue,".RData"))


pdf(paste0("Plots/total-intensity-QC-",tissue,".pdf"),height=8,width=12)
plot(colSums(qc.pep,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total~intensity (millions)",line = 2.5)
title(main="Total intensity for samples and pools by plex")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
legend("bottomleft",legend = c("ASD","CONTROL","POOL"),pch=pchDX,col="grey75",bty="n")
dev.off()

pdf(paste0("Plots/fraction-called-QC-",tissue,".pdf"),height=8,width=12)
plot(colSums(!is.na(qc.pep))/nrow(qc.pep),col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
title(ylab="Fraction",line = 3)
title(main="Fraction peptides with intensities")
legend("bottomleft",legend = c("ASD","CONTROL","POOL"),pch=pchDX,col="grey75",bty="n")
dev.off()

pdf(paste0("Plots/MDS-QC-",tissue,".pdf"),height = 8, width =8)
plotMDS(log2(qc.pep),col=colPLEX[meta.samples$PLEX],las=1,main="MDS of QC log2(intensity) data",dim.plot = c(1,2))
legend("topleft",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
dev.off()






