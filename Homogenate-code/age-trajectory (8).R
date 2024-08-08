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

tissue="HOMOGENATE-PROTEIN"

### load thenormalized data and the informaation on the model
load(paste0("WorkData/MODEL-SELECTION-",tissue,".RData"))

# functions needed
lmerResiduals <- function(feature,Y,XX,covariates,random.effect,add.back){
  # function performs a linear regression given a set of covariates and returns the residuals
  # add.back are the effect as specified in the lmer coefficient that you want to backin (use NULL if you don't want to anything back)
  model.formula<-formula(paste("Y['",feature,"',] ~ ",paste(c(covariates),collapse=" + "),paste0(" + (1|",random.effect,")"),sep=""))
  re.lmer=lmer(model.formula,data=XX,REML = T,control = lmerControl(optimizer="bobyqa"))
  re.residuals=residuals(re.lmer)
  re.coef=summary(re.lmer)$coefficients
  
  # add back the effect of interest
  if(!is.null(add.back)){
    # set up the matrix for the covariates
    covMatrix=model.matrix(formula(paste("~",paste(selected.model,collapse=" + "))),data=XX)
    # find the regression coefficients of interest
    beta=as.matrix(re.coef[add.back,"Estimate"])
    # determine the adjustment
    cov.adjust=covMatrix[,add.back]%*%beta
    # add this effect back in
    re.residuals=re.residuals+cov.adjust
  }
  
  # add back the add.back effect
  return(data.frame(t(re.residuals)))
}


### calcualte the residuals adding back the effect of age and age:dx interaction
system.time(RESage <- mclapply(rownames(norm.prt),lmerResiduals,Y=log2(norm.prt[,samples]),XX=meta.samples[samples,],
                                covariates=selected.model,random.effect="PAIR",add.back=c("AGEyr"),mc.cores=64,mc.preschedule=T))
RESage=as.data.frame(rbindlist(RESage))
colnames(RESage)=colnames(norm.prt)
rownames(RESage)=rownames(norm.prt)


# load the data to 
#load("WorkData/SURROGATE-VARIABLES-HOMOGENATE.RData")


# set up the dataframe to use for flexmix.
asd=rownames(meta.samples)[meta.samples$DX == "ASD"]
ctl=rownames(meta.samples)[meta.samples$DX == "CTL"]
proteins=rownames(norm.prt)

# use only the the controls
df=data.frame(protein=rep(proteins,each=length(ctl)),sample=rep(ctl,times=length(proteins)))
df[,c("AGEyr")]=meta.samples[df$sample,c("AGEyr")]
df[,"RESage"]=RESage[cbind(df$protein,df$sample)]

# run flexmix to find the mixture distribution, this can take some time
re.flexmix <- flexmix(RESage ~ AGEyr  | protein, k = 3, 
                      data = df, control = list(tolerance = 10^-5),model=FLXMRglmfix(varFix = T))

(n.tbl=table(clusters(re.flexmix))/length(ctl))
#1    2    3 
#3230 1036 461

df$cluster=clusters(re.flexmix)
prt.clusters=df[!duplicated(df[,c("protein","cluster")]),c("protein","cluster")]


df.param=parameters(re.flexmix)
df.param=df.param[rownames(df.param) != "sigma",]
df.age=df.param[grep("coef.AGE",rownames(df.param)),]
cov.matrix=data.frame(intercept=1,age=seq(min(df$AGEyr),max(df$AGEyr),1))  
y.age=as.matrix(cov.matrix)%*%as.matrix(df.param)  # intercept plus age effect

pdf("HOMOGENATE-PROTEIN-AGE-mixture-modeling_11-29-2022.pdf",height=8,width=8)
plot(NULL,NULL,xlim=range(cov.matrix$age),ylim=range(y.age),las=1,xlab="Age",ylab="Effect")
abline(h=0,col="grey50",lty=3)
for(i in 1:ncol(y.age)){
  points(cov.matrix$age,y.age[,i],type="l",col=viridis(9)[c(2,4,6)][i],lwd=2,lty=i)
}

plus.minus=c("-","+","+");names(plus.minus)=c("-1","0","1")
legend("topleft",legend=c(sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param[1,1],2),plus.minus[as.character(sign(df.param[2,1]))],abs(df.param[2,1]),n.tbl[1]),
                             sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param[1,2],2),plus.minus[as.character(sign(df.param[2,2]))],abs(df.param[2,2]),n.tbl[2]),
                             sprintf("abnd = %4.2f %s %4.3f * age; n.prt: %4i",round(df.param[1,3],2),plus.minus[as.character(sign(df.param[2,3]))],abs(df.param[2,3]),n.tbl[3])),
                              lty=1:ncol(y.age),lwd = 2,pch=NULL,col=viridis(9)[c(2,4,6)][1:ncol(y.age)],bty="n")
dev.off()


RE=list(parameters=as.data.frame(df.param),clustering=prt.clusters)
WriteXLS::WriteXLS(RE,"HOMOGENATE-PROTEIN-AGE-mixture-modeling_11-29-2022.xlsx")


save(re.flexmix,meta.samples,meta.pep,norm.prt,selected.model,RESage,RE,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file="WorkData/TRAJECTORY-HOMOGENATE.RData")

load("WorkData/TRAJECTORY-HOMOGENATE.RData")

##### calculate the slope for both asd and ctl

linearRegression<-function(feature,yy,xx){
  lm.formula=as.formula(paste0("yy['",feature,"',]","~","xx"))
  re.lm=summary(lm(lm.formula))
  a=re.lm$coef[1,1]
  b=re.lm$coef[2,1]
  se=re.lm$coef[2,2]
  return(data.frame(feature,a,b,se))
}
lrCTL=as.data.frame(rbindlist(lapply(prt.clusters$protein,linearRegression,yy=as.matrix(RESage[,ctl]),xx=meta.samples[ctl,"AGEyr"])))
lrASD=as.data.frame(rbindlist(lapply(prt.clusters$protein,linearRegression,yy=as.matrix(RESage[,asd]),xx=meta.samples[asd,"AGEyr"])))

prt.clusters=RE$clustering
rownames(prt.clusters)=prt.clusters$protein
prt.clusters[lrASD$feature,c("a.asd","b.asd","se.asd")]=lrASD[,c("a","b","se")]
prt.clusters[lrCTL$feature,c("a.ctl","b.ctl","se.ctl")]=lrCTL[,c("a","b","se")]


plot(prt.clusters[prt.clusters$cluster == 3,"b.ctl"],prt.clusters[prt.clusters$cluster == 3,"b.asd"])
re.lm=lm(prt.clusters[prt.clusters$cluster == 3,"b.asd"]~prt.clusters[prt.clusters$cluster == 3,"b.ctl"])
abline(a=0,b=1)
abline(re.lm$coef)

t.test(prt.clusters[prt.clusters$cluster == 1,"b.asd"],prt.clusters[prt.clusters$cluster == 1,"b.ctl"],paired = T)
t.test(prt.clusters[prt.clusters$cluster == 2,"b.asd"],prt.clusters[prt.clusters$cluster == 2,"b.ctl"],paired = T)
t.test(prt.clusters[prt.clusters$cluster == 3,"b.asd"],prt.clusters[prt.clusters$cluster == 3,"b.ctl"],paired = T)

hist(prt.clusters[prt.clusters$cluster == 2,"b.asd"]-prt.clusters[prt.clusters$cluster == 2,"b.ctl"])


regr.coef=aggregate(cbind(a.asd,b.asd,a.ctl,b.ctl)~cluster,FUN="mean",data=as.data.frame(prt.clusters))

pdf("HOMOGENATE-AGE-TRAJ-ASD-NT-2022-11-29.pdf",height=8,width=8)
plot(NULL,NULL,xlim=c(4,32),ylim=c(-.5,.75),las=1,xlab="Age",ylab="Effect")
abline(h=0,col="grey50")
for(i in 1:3){
  abline(a=regr.coef[i,"a.ctl"],b=regr.coef[i,"b.ctl"],lty=2,col=viridis(9)[c(2,5,8)][i],lwd=2)
  abline(a=regr.coef[i,"a.asd"],b=regr.coef[i,"b.asd"],lty=1,col=viridis(9)[c(2,5,8)][i],lwd=2)
}
legend("topleft",legend=c("NT","ASD"),lty=c(2,1),lwd=2,bty="n")
legend("bottomleft",legend=paste0(c("stable","decreasing","increasing"),"; n = ",n.tbl),col=viridis(9)[c(2,5,8)],lty=1,lwd=2,bty="n")
dev.off()


MLM=as.data.frame(read_excel("Results/HOMOGENATE-PROTEIN.xlsx",sheet="MLM.DX"))
rownames(MLM)=MLM$protein


posterior=re.flexmix@posterior$scaled[seq(1,nrow(re.flexmix@posterior$scaled),31),]
prt.clusters[,c("pp1","pp2","pp3")]=posterior
t.test(prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .75 & prt.clusters$pp2 <= .90,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .75 & prt.clusters$pp2 <= .90,"b.asd"],paired = T)

t.test(prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .90 & prt.clusters$pp2 <= 1,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 2 & prt.clusters$pp2 > .90 & prt.clusters$pp2 <= 1,"b.asd"],paired = T)


t.test(prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .75 & prt.clusters$pp3 <= .90,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .75 & prt.clusters$pp3 <= .90,"b.asd"],paired = T)

t.test(prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .90 & prt.clusters$pp3 <= 1,"b.ctl"],
       prt.clusters[prt.clusters$cluster == 3 & prt.clusters$pp3 > .90 & prt.clusters$pp3 <= 1,"b.asd"],paired = T)



prt.clusters$delta=prt.clusters$b.asd-prt.clusters$b.ctl
prt.clusters$se.delta=sqrt(prt.clusters$se.asd^2+prt.clusters$se.ctl^2)
prt.clusters$p.delta=2*pt(-abs(prt.clusters$delta/prt.clusters$se.delta),df=27)
prt.clusters$q.delta=qvalue(prt.clusters$p.delta)$qvalue

table(prt.clusters$cluster[prt.clusters$q.delta < 0.05])

prt.clusters[MLM$protein,"dx.fc"]=MLM$Estimate

WriteXLS::WriteXLS(prt.clusters,paste0("HOMOGENATE-AGE-TRAJ-ASD-NT-",Sys.Date(),".xls"))

#### examples
which(prt.clusters$cluster == 2 & prt.clusters$b.ctl < -0.02 & abs(prt.clusters$b.asd) < 0.005)

k=1867
prt.clusters[k,]
plot(NULL,NULL,xlim=c(4,32),ylim=c(-1,1),las=1,xlab="Age",ylab="Effect")
abline(h=0,col="grey50")
abline(a=prt.clusters[k,"a.ctl"],b=prt.clusters[k,"b.ctl"],lty=2,col=viridis(9)[c(2,5,8)][2],lwd=2)
abline(a=prt.clusters[k,"dx.fc"],b=prt.clusters[k,"b.asd"],lty=1,col=viridis(9)[c(2,5,8)][2],lwd=2)
