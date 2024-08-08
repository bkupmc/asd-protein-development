rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="SYNAPTOSOME-PEPTIDE"

############################################
### This versions is designed for protein level wort
### An initial step ius performed in which peptides are rolled up into proteins
### IN THIS VERSION ONLY THE SAMPLE LOADING NORMALIZATION IS PERFORMED
#####################################


#### this code performs the three  normalization steps as used by Ying
#### Normalizations are performed on the observed scale

#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(limma)


### load the VAE imptued data
load(paste0("WorkData/VAE-",tissue,".RData"))

##### ROLL the peptides into proteins
# only use the peptides that map uniquely to a single protein
i.multi=grep(";",meta.pep$ROLLUP) # 1820
if(length(i.multi) > 0){meta.pep=meta.pep[-i.multi,]} # 27,491
length(unique(meta.pep$ROLLUP)) # 4287

vae.pep=vae.pep[rownames(meta.pep),]
vae.prt=aggregate(vae.pep~meta.pep$ROLLUP,FUN="sum")
rownames(vae.prt)=vae.prt[,1]
vae.prt=as.matrix(vae.prt[,-1])
dim(vae.prt)

tissue="SYNAPTOSOME-PROTEIN"

##### PLOTS
pdf(paste0("Plots/total-intensity-ROLLEDUP-",tissue,".pdf"),height=8,width=12)
plot(colSums(vae.prt,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total-intensity (millions)",line = 2.5)
title(main="Total intensity for samples by plex")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
legend("bottomleft",legend = levels(meta.samples$DX),pch=pchDX,col="grey75",bty="n")
dev.off()


pdf(paste0("Plots/MDS-ROLLEDUP-",tissue,".pdf"),height = 8, width =12)
par(mfrow=c(1,2))
par(oma=c(0,0,2,0))
par(mar=c(5,5,1,2))
plotMDS(log2(vae.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(1,2))
legend("topright",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
plotMDS(log2(vae.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(3,4))
mtext(text="MDS of normalized log2(intensity) data",side=3,outer=T,cex=1.5,line=0)

dev.off()






##########
### Sample loading normalization
average.sum=sum(vae.prt)/ncol(vae.prt)
channel.sum=colSums(vae.prt)
sl.prt=t(t(vae.prt)*average.sum/channel.sum)

pdf(paste0("Plots/MDS-SL-NORM-",tissue,".pdf"),height = 8, width =8)
plotMDS(log2(sl.prt),col=colPLEX[meta.samples$PLEX],las=1,main="MDS of SL normalized log2(intensity) data",dim.plot = c(1,2))
legend("topleft",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
dev.off()

####### Normalized peptides
norm.prt=sl.prt   # 342,763,333


#### save the resulting data
save(meta.samples,meta.pep,norm.prt,sl.prt,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/NORMALIZED-",tissue,".RData"))



######PLOTS






par(mfrow=c(1,1))


pdf(paste0("Plots/total-intensity-NORMALIZED-",tissue,".pdf"),height=8,width=12)
plot(colSums(norm.prt,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total~intensity (millions)",line = 2.5)
title(main="Total intensity for samples by plex")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
legend("bottomleft",legend = levels(meta.samples$DX),pch=pchDX,col="grey75",bty="n")
dev.off()


pdf(paste0("Plots/MDS-NORMALIZED-",tissue,".pdf"),height = 8, width =12)
par(mfrow=c(1,2))
par(oma=c(0,0,2,0))
par(mar=c(5,5,1,2))
plotMDS(log2(norm.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(1,2))
legend("topright",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
plotMDS(log2(norm.prt),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(3,4))
mtext(text="MDS of normalized log2(intensity) data",side=3,outer=T,cex=1.5,line=0)

dev.off()


# boxplot(log2(norm.pep),range=3,las=2,xaxt="n",col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
# abline(v=plex.border,lty=3,col="grey50")
# axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))






