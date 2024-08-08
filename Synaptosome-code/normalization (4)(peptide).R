rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="SYNAPTOSOME-PEPTIDE"

############################################
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

##########
### Sample loading normalization
average.sum=sum(vae.pep)/ncol(vae.pep)
channel.sum=colSums(vae.pep)
sl.pep=t(t(vae.pep)*average.sum/channel.sum)

pdf(paste0("Plots/MDS-SL-NORM-",tissue,".pdf"),height = 8, width =8)
plotMDS(log2(sl.pep),col=colPLEX[meta.samples$PLEX],las=1,main="MDS of SL normalized log2(intensity) data",dim.plot = c(1,2))
legend("topleft",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
dev.off()



####### Normalized peptides
norm.pep=sl.pep    # 402,708,213


#### save the resulting data
save(meta.samples,meta.pep,norm.pep,sl.pep,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/NORMALIZED-",tissue,".RData"))



######PLOTS


pdf(paste0("Plots/total-intensity-NORMALIZED-",tissue,".pdf"),height=8,width=12)
plot(colSums(norm.pep,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
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
plotMDS(log2(norm.pep),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(1,2))
legend("bottomleft",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
plotMDS(log2(norm.pep),col=colPLEX[meta.samples$PLEX],las=1,main="",dim.plot = c(3,4))
mtext(text="MDS of normalized log2(intensity) data",side=3,outer=T,cex=1.5,line=0)

dev.off()







