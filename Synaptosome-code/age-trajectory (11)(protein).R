rm(list=ls()); gc()
options(stringsAsFactors = F)

# this portion performs a tracjectory analysis for the age effects

require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(limma)
require(sva)
require(flexmix)

#lcm1=flexmix(.~.|veh, model=FLXMRlmer(formula = chargedFlag ~ lastTripDummy + endsOnWeekend + endHome, random=~1, weighted = FALSE), data=trips4,
#             + control = list(minprior = 0.2), k=6)

# load the data to use
load("WorkData/SURROGATE-VARIABLES-HOMOGENATE.RData")

#1) fit the full model including the surrogate variables
X=model.matrix(as.formula(paste("~",paste(selected.model.sv,collapse="+"),sep = "")),data=meta.samples.sv)   # design matrix X
Y=log2(norm.pep[,samples])                                                                                   # observation matrix Y

BETA=t(solve(t(X)%*%X,t(Y%*%X)))                                                                   # calcualte solutions for all features
RES=Y-BETA%*%t(X)                                                                                  # determine the residuals
RESage=RES+BETA[,grep("AGE",colnames(BETA))]%*%t(X[,grep("AGE",colnames(X))])                      # add the age effects back in 

#### Run flexmix on the residuals using age and age^2 as covariates
# create an effect for the interaction term
meta.samples.sv$`AGEyr:DX`=as.data.frame(X)$`AGEyr:DXASD`

df=data.frame(peptide=rep(rownames(RESage),each=ncol(RESage)),sample=rep(colnames(RESage),times=length(rownames(RESage))))
df[,c("AGEyr","AGEyr2")]=meta.samples.sv[df$sample,c("AGEyr","AGEyr2")]
df[,"RESage"]=RESage[cbind(df$peptide,df$sample)]

# run flexmix to find the mixture distribution, this can take some time
re.flexmix <- flexmix(RESage ~ AGEyr + AGEyr2 | peptide, k = 3, 
                      data = df, control = list(tolerance = 10^-3),model=FLXMRglmfix(varFix = T))

table(clusters(re.flexmix))/length(samples)
#    1    2    3
# 7432  6351  19099 

df.param=parameters(re.flexmix)
df.age=df.param[grep("coef.AGE",rownames(df.param)),]
x.age=poly(seq(4,32),2,raw=T)
y.age=t(t(x.age%*%df.age)+df.param[1,])  # intercept plus age effect

plot(NULL,NULL,xlim=c(4,32),ylim=range(y.age),las=1,xlab="Age",ylab="Effect")
for(i in 1:ncol(y.age)){
  points(x.age[,1],y.age[,i],type="l",col=viridis(3)[i],lwd=2,lty=i)
}
abline(h=0,col="grey50",lty=3)
legend("bottomleft",legend=paste("Trajectory",1:ncol(y.age),sep="-"),lty=1:ncol(y.age),lwd = 2,pch=NULL,col=viridis(3)[1:ncol(y.age)],bty="n")

save(re.flexmix,meta.samples,meta.samples.sv,meta.pep,norm.pep,selected.model,selected.model.sv,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file="WorkData/TRAJECTORY-HOMOGENATE.RData")
