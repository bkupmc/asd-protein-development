rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="HOMOGENATE-PEPTIDE"

### this code performs the VAE imputation using softImpute to come up with the initialization values
### the code should not have to be changed other than when one wants to use another initialization process
### the user needs to specify the covariates to use in the VAE
#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(softImpute)
require(limma)

#### copies of VAEIT and vae-python.R need to be in this directory

############## function
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


##### load the data
load(paste0("WorkData/QC-",tissue,".RData"))

# simulate a dataset with missingness
# order the peptides by average expression
seed=0
set.seed(seed)
qc.pep=qc.pep[order(rowMeans(qc.pep,na.rm=T)),]
Xall=qc.pep      # all data
Xcomplete=qc.pep[rowSums(is.na(qc.pep)) == 0,] # complete data

# select the pattern to simulate from Xall
i.xall=sort(sample(nrow(Xall),nrow(Xcomplete)))
# find the locations with missing values in Xcomplete for the selected peptides
i.na=which(is.na(Xall[i.xall,]),arr.ind = F)
# now create a matrix with missing vales based on this
Xna=Xcomplete                               # start with the complete data
Xna[i.na]=NA                                # add the missingnes




#1) find starting values using softImpute

# find the maximum value for lambda. This defines the searched range
(lam_max = lambda0(log2(Xna), lambda = 0, maxit = 100, trace.it  = FALSE, thresh = 1e-05))  # 11052.99

# find the optimal lambda, using a rank = 1 (r=1) softImpute
lambda.sv=NULL
for(i in 1:10){
  print(i); flush.console()
  re.optimize=optimize(mseImputation,interval = c(0,lam_max),tol = 1e-8,
                       r=1, sample.col=samples, pool.col=pools, Cna=log2(Xna), thresh=1e-5, maxit=100)
  lambda.sv=c(lambda.sv,re.optimize$minimum)
}
mean(lambda.sv)
lambda=round(mean(lambda.sv),2) # 0.35



# perform softimpute given lambda and r=1
re.softimpute=softImpute(log2(Xna), rank.max = 1, lambda=lambda, type="als", trace.it = F, thresh = 1e-5, maxit=100, final.svd=T)
# all imputed/smoothed values
Xsi=complete(log2(Xna), re.softimpute)
Xsi=2^Xsi  # put back on original scale

mean((log2(Xsi)-log2(Xna))^2,na.rm=T)
mean((log2(Xsi)-log2(Xcomplete))^2)
mean((log2(Xsi[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 1.844

###### select samples only
meta.samples=meta.samples[samples,]
meta.samples$DX=factor(meta.samples$DX,levels=levels(meta.samples$DX)[1:2])
Xcomplete=Xcomplete[,samples]
Xna=Xna[,samples]
Xsi=Xsi[,samples]


#####################################
# now perform the VAE step
####################################

##### Instructions to perform setting up an environment
# conda install mamba -n base -c conda-forge
# conda config --add channels conda-forge
# mamba create -n py_aev python==3.8 -y   # run in shell
# conda activate py_aev
# mamba install -n py_env -c conda-forge -c esri tensorflow==2.6 tensorflow-probability==0.14 tensorflow-addons scikit-learn numpy -y




# start conda environment
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
Yinit=as.matrix(log2(Xsi))
Yna=as.matrix(log2(Xna))

re.vae<-vae(Xna=Yna, batches_cate=batches_cate, batches_cont=batches_cont, X_init=Yinit, 
            num_epoch=20L, print_every_epoch=1L, X=NULL,
            beta_kl=1., beta_unobs=0.5, p_feat=.2)

Xvae=2^re.vae$X_imp

mean((log2(Xsi)-log2(Xna))^2,na.rm=T)
mean((log2(Xvae)-log2(Xcomplete))^2)
mean((log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 1.7828
mean((log2(Xvae[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)]))^2)  # 0.0009

median((log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]))^2)  # 0.5994
median((log2(Xvae[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)]))^2)  # 0.0001

save(seed,Xall,Xcomplete,Xsi,Xvae,Xna,re.softimpute,re.vae,file="WorkData/VAE-TEST-HOMOGENATE.RData")





plot(log2(Xcomplete[300,]),col=viridis(3)[3],pch=15)
points(log2(Xna[300,]),col=viridis(3)[2],pch=17)
points(log2(Xvae[300,]),col=viridis(3)[1],pch=19)

hist(abs(log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)])),breaks=seq(0,9,.1))
max(abs(log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)])))
hist(log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)]),breaks=seq(-9,9,.1))

mean(abs(log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)])))
median(abs(log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)])))

mean(abs(log2(Xvae[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)])))
median(abs(log2(Xvae[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)])))



mean((log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)])))
median((log2(Xvae[is.na(Xna)])-log2(Xcomplete[is.na(Xna)])))

mean((log2(Xvae[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)])))
median((log2(Xvae[!is.na(Xna)])-log2(Xcomplete[!is.na(Xna)])))

save(seed,Xall,Xcomplete,Xsi,Xvae,Xna,re.softimpute,re.vae,file="WorkData/VAE-TEST-HOMOGENATE.RData")

