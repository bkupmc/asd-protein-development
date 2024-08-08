rm(list=ls()); gc()
options(stringsAsFactors = F)
tissue="HOMOGENATE-PEPTIDE"
#####################################################################################################
#### THIS VERSION OF THE CODE ONLY PERFORMS THE IMPUTATION FOR THE SAMPLES
#### POOLS WILL BE REMOVED FROM THE ANALYSES


### this code performs the VAE imputation using softImpute to come up with the initialization values
### the code should not have to be changed other than when one wants to use another initialization process
### the user needs to specify the covariates to use in the VAE
#### libraries use
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(softImpute)
require(limma)

#### copies of VAEIT and vae-python.R need to be in this directory


##### load the data
load(paste0("WorkData/QC-",tissue,".RData"))



#1) find starting values using softImpute
mseImputation <- function(lambda,r,sample.col,pool.col,Cna,thresh=1e-4,maxit=100){
  # function that finds the lambda for a given r using the minimum distance between the observed
  # observed average intensity and the imputed average intensity
  # lambda - lambda to use in the imputation
  # rank - rank to use for the svd decomposition
  # Cna - matrix with missing values    
  re.softimpute<-softImpute(Cna,rank.max = r,lambda=lambda,type="als",trace.it  = F,thresh = thresh, maxit = maxit)
  if(r > 1){
    Cudv=(re.softimpute$u%*%diag(re.softimpute$d))%*%t(re.softimpute$v)
  }else{
    Cudv=(re.softimpute$u*re.softimpute$d)%*%t(re.softimpute$v)
  }
  colnames(Cudv)=colnames(Cna)
  rownames(Cudv)=rownames(Cna)
  
  i.complete=which(rowSums(is.na(Cna)) == 0)
  mse=mean((rowMeans(Cudv[i.complete,sample.col],na.rm=T)-rowMeans(Cna[i.complete,pool.col],na.rm=T))^2)
  
  return(mse)
}

# find the missing values in qc.pep
ina <- which(is.na(qc.pep))  # 540,972

# find the maximum value for lambda. This defines the searched range
(lam_max = lambda0(log2(qc.pep), lambda = 0, maxit = 100, trace.it  = FALSE, thresh = 1e-05))  # 17747.96

# find the optimal lambda, using a rank = 1 (r=1) softImpute, needs to be written to a function to increase the speed.
lambda.sv=NULL
for(i in 1:10){
  print(i); flush.console()  
  re.optimize=optimize(mseImputation,interval = c(0,lam_max),tol = 1e-8,
                       r=1, sample.col=samples, pool.col=pools, Cna=log2(qc.pep), thresh=1e-5, maxit=100)
  lambda.sv=c(lambda.sv,re.optimize$minimum)
}
mean(lambda.sv)
lambda=round(mean(lambda.sv),2)  # 0.53

# perform softimpute given lambda and r=1
re.softimpute=softImpute(log2(qc.pep), rank.max = 1, lambda=lambda, type="als", trace.it = F, thresh = 1e-5, maxit=100, final.svd=T)
# all imputed/smoothed values
si.pep=complete(log2(qc.pep), re.softimpute)
si.pep=2^si.pep  # put back on original scale

#################################
###### select samples only
#################################
meta.samples=meta.samples[samples,]
meta.samples$DX=factor(meta.samples$DX,levels=levels(meta.samples$DX)[1:2])
qc.pep=qc.pep[,samples]
si.pep=si.pep[,samples]

# now some things for for plotting
plex.border=cumsum(table(meta.samples$PLEX))
plex.border=plex.border[-length(plex.border)]+.5
plex.mid=c(0,cumsum(table(meta.samples$PLEX)))
plex.mid=.5+(plex.mid[-length(plex.mid)]+plex.mid[-1])/2
names(plex.mid)=plexes

pchDX=c(15,17);names(pchDX)=unique(meta.samples$DX)
colPLEX=viridis(length(unique(meta.samples$PLEX))); names(colPLEX)=unique(meta.samples$PLEX)


#####################################
# now perform the VAE step
####################################

##### Instructions to perform setting up an environment
# conda install mamba -n base -c conda-forge
# conda config --add channels conda-forge
# mamba create -n py_vae python==3.8 -y   # run in shell
# conda activate py_vae
# mamba install -n py_vae -c conda-forge -c esri tensorflow==2.6 tensorflow-probability==0.14 tensorflow-addons scikit-learn numpy -y

# start conda environment
# make sure a copy of vae-python.R and the directory VAEIT are in the working directory
if(!dir.exists("VAEIT")){
  print("move VAE software to the current directory")
  system("cp -r /data3/Software/VAEIT/VAEIT VAEIT" )
  system("cp /data3/Software/VAEIT/vae_python.r vae_python.r")
}else{
  print("copy of VAE software located in the current directory")
}

### library reticulate is needed to run the python code
require(reticulate)
use_condaenv(condaenv = "py_vae") 
# load the vae code
source("vae_python.r")

# set up the model
# categorical covariates
batches_cate <- meta.samples[,c("PLEX","DX","RACE","SEX")]
batches_cate <- as.matrix(batches_cate)
# continuous covariates
batches_cont <- meta.samples[,c("AGEyr","AGEyr2","PMIhr")]
batches_cont <- as.matrix(batches_cont)

# create the input matrices
Xinit=as.matrix(log2(si.pep))
Xna=as.matrix(log2(qc.pep))

re.vae<-vae(Xna=Xna, batches_cate=batches_cate, batches_cont=batches_cont, X_init=Xinit, 
            num_epoch=20L, print_every_epoch=1L, X=NULL,
            beta_kl=1., beta_unobs=0.5, p_feat=.2)


save(re.vae,qc.pep,si.pep,batches_cate,batches_cont,file=paste0("WorkData/VAE-output-",tissue,".RData"))
vae.pep=2^re.vae$X_imp                  # put back onto observed scale
rownames(vae.pep)=rownames(qc.pep)
colnames(vae.pep)=colnames(qc.pep)

save(meta.samples,meta.pep,vae.pep,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/VAE-",tissue,".RData"))

mean((log2(qc.pep)-log2(vae.pep))^2,na.rm=T)   # 0.0005   0.0018
median((log2(qc.pep)-log2(vae.pep))^2,na.rm=T) # 0.0001   0.0004

######PLOTS

pdf(paste0("Plots/total-intensity-VAE-",tissue,".pdf"),height=8,width=12)
plot(colSums(vae.pep,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total~intensity (millions)",line = 2.5)
title(main="Total intensity for samples by plex")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
legend("bottomleft",legend = levels(meta.samples$DX),pch=pchDX,col="grey75",bty="n")
dev.off()

pdf(paste0("Plots/fraction-called-VAE-",tissue,".pdf"),height=8,width=12)
plot(colSums(!is.na(vae.pep))/nrow(vae.pep),col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
title(ylab="Fraction",line = 3)
title(main="Fraction peptides with intensities")
legend("bottomleft",legend = levels(meta.samples$DX),pch=pchDX,col="grey75",bty="n")
dev.off()

pdf(paste0("Plots/MDS-VAE-",tissue,".pdf"),height = 8, width =8)
plotMDS(log2(vae.pep),col=colPLEX[meta.samples$PLEX],las=1,main="MDS of imputed log2(intensity) data",dim.plot = c(1,2))
legend("bottomleft",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
dev.off()








