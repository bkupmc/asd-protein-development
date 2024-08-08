rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="SYNAPTOSOME-PEPTIDE"

#### Here the data is read from the PSM files and manipulated to create a set of standard input data.frames that are used throughout the analysis
#### since this tends to be different for each analysis, this code needs tp be changed to meet the project's needs
#### at the end there should be a series of data.frames (names can be customized to individual needs)
#### meta.samples - this contains
####                 - ID for each sample and pool 
####                 - ROLE identifier of sample or pool
####                 - PLEX and CHANNEL 
####                 - any covariates thought to be of interest for the analysis
####                       - continous   - factors
####                       - categorical - factors
####                 - rownames are the ID
####                 - colnames are ID, ROLE, PLEX, CHANNEL, names of covariates
#### meta.feat   - this contains as minimum information on the feature (peptide) and the rollup (protein)
####                  - FEATURE
####                  - ROLLUP
#### intensity     - matrix with the features as rows and samples and pools as columns
####               - rownames are the feature names
####               - colnames are the sample and pool ID from the meta.samples data.frame

#### libraries used throughout the pipeline
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(softImpute)
require(limma)
require(sva)

###################################
########## sample meta file
######################################
samples=fread("/data3/DownLoadedData/McDonald/ASD_Development_June2022/Synaptosome/Run1/Shelby_Prelim_Analysis/Peptide_Normalization_ProteinRollup/Covariates/ASD_Dev_limma.csv",
           header=T,data.table=F)
colnames(samples)
colnames(samples)=c("ID1","AGEyr","PMIhr","RACE","SEX","DX","COMORBIDITIES","ID","PAIR")
samples$ID=paste0("S",samples$ID)
samples$PAIR=paste0("MP",samples$PAIR)
rownames(samples)=samples$ID1  # this is the ID information in the plex information file

# read the plex information for the samples
plex=fread("/data3/DownLoadedData/McDonald/ASD_Development_June2022/Synaptosome/Run1/Shelby_Prelim_Analysis/Peptide_Normalization_ProteinRollup/Peptide_input_files/ASD_Dev_Clinic_PD.csv",
           header=T,data.table=F)
colnames(plex)=c("PD.ID","PLEX","TMT.PLEX","CHANNEL","RYAN.PLEX","RYAN.CHANNEL","ID","ID1")
rownames(plex)=plex$ID1
sum(rownames(plex) %in% rownames(samples))  # check to see if the IDs match between


#### create the meta file for the samples and pools
meta.samples=plex[,c("PD.ID","PLEX","CHANNEL")]
# add ROLE
meta.samples$ROLE="SAMPLE"                                           #56
meta.samples$ROLE[substr(rownames(meta.samples),1,1) == "P"]="POOL"  #18
# add the information from samples
i.samples=match(rownames(meta.samples),rownames(samples))
meta.samples=cbind.data.frame(meta.samples,samples[i.samples,c("ID","SEX","DX","AGEyr","PMIhr","RACE","COMORBIDITIES","PAIR")])

### turn categorical variables into factors
meta.samples$SEX=factor(c("Male","Female")[meta.samples$SEX], levels=c("Male","Female"))

meta.samples$DX[is.na(meta.samples$DX)]=3
meta.samples$DX=factor(c("CTL","ASD","POOL")[meta.samples$DX], levels=c("CTL","ASD","POOL"))

meta.samples$RACE=factor(c("White","Black or African-American","Unknown")[meta.samples$RACE],levels=c("White","Black or African-American","Unknown"),labels=c("CAU","AFR","UNK"))

meta.samples$PLEX=factor(meta.samples$PLEX, levels=unique(meta.samples$PLEX))

#meta.samples$COMORBIDITIES[is.na(meta.samples$COMORBIDITIES) & meta.samples$ROLE == "SAMPLE"]="NONE"
#meta.samples$COMORBIDITIES=factor(c("NONE","ADHD","CP","EPLSY","PW")[meta.samples$COMORBIDITIES],levels=c("NONE","ADHD","CP","EPLSY","PW"))

# assign the pools an ID
meta.samples$ID[is.na(meta.samples$ID)]=paste0("P",1000*as.numeric(meta.samples$PLEX[is.na(meta.samples$ID)])+(1:sum(is.na(meta.samples$ID))))

# add a quadratic age effect
meta.samples$AGEyr2=meta.samples$AGEyr^2

######################################################
###### meta data for the features
######################################################

# read the psm
psm=as.data.frame(read_excel("/data3/DownLoadedData/McDonald/ASD_Development_June2022/Synaptosome/Run1/PD_export_originalfiles/PSM/ASD_Syn_complete_noP5F6_Run1_PSM.xlsx")) # 291,424
head(psm)


# create a meta file for the peptides
meta.pep=data.frame(FEATURE=unique(psm$`Annotated Sequence`))
meta.pep[,c("ROLLUP")]=psm[match(meta.pep$FEATURE,psm$`Annotated Sequence`),c("Master Protein Accessions")]
rownames(meta.pep)=meta.pep$FEATURE  # 61,060



#### now set up the intensity data
# read the psm
psm=as.data.frame(read_excel("/data3/DownLoadedData/McDonald/ASD_Development_June2022/Synaptosome/Run1/PD_export_originalfiles/PSM/ASD_Syn_complete_noP5F6_Run1_PSM.xlsx")) # 291,424

# set up a matrix to store the peptide level intensities (roll-up)
raw.pep=matrix(NA,nrow=nrow(meta.pep),ncol=nrow(meta.samples))
rownames(raw.pep)=meta.pep$FEATURE
colnames(raw.pep)=meta.samples$ID


##### process by plex (PD.FRACTION)
plexes=unique(meta.samples$PLEX)
for(plex in plexes){
  print(plex); flush.console()
  i.plex=which(substr(psm$`File ID`,1,3) == plex)
  plex.psm=psm[i.plex,c(5,31:46)]
  # generate column id based on sample ID
  i.id=match(colnames(plex.psm),paste("Abundance:",meta.samples$CHANNEL[meta.samples$PLEX == plex]))
  colnames(plex.psm)=meta.samples$ID[meta.samples$PLEX== plex][i.id]
  colnames(plex.psm)[1]="AnnotatedSequence"
  
  # remove rows with all NA
  plex.psm=plex.psm[rowSums(!is.na(plex.psm[,-1])) > 0,]
  
  # remove columns with NA name
  plex.psm=plex.psm[,!is.na(colnames(plex.psm))]
  
  # remove rows with all NA, this is needed for the aggregate
  plex.psm=plex.psm[rowSums(!is.na(plex.psm[,-1])) > 0,]
  
  
  # now perform the roll-up of psm values into peptide values
  ru.psm=aggregate(as.matrix(plex.psm[,-1])~plex.psm$AnnotatedSequence,FUN="sum",na.rm=T,na.action = "na.pass")
  rownames(ru.psm)=ru.psm[,1]
  ru.psm=as.matrix(ru.psm[,-1])
  ru.psm[ru.psm == 0]=NA  # aggregate introduces 0 during this operation. The original input does not contain 0
  
  # add this to the petide matrix
  raw.pep[rownames(ru.psm),colnames(ru.psm)]=ru.psm
}


### some additional information that will be used throughout the pipeline
rownames(meta.samples)=meta.samples$ID
meta.samples=meta.samples[,-1]             # remvoe the PD.ID
samples=rownames(meta.samples)[meta.samples$ROLE == "SAMPLE"]
pools=rownames(meta.samples)[meta.samples$ROLE == "POOL"]
plexes=unique(meta.samples$PLEX)

# now some things for for plotting
plex.border=cumsum(table(meta.samples$PLEX))
plex.border=plex.border[-length(plex.border)]+.5
plex.mid=c(0,cumsum(table(meta.samples$PLEX)))
plex.mid=.5+(plex.mid[-length(plex.mid)]+plex.mid[-1])/2
names(plex.mid)=plexes

pchDX=c(15,17,19);names(pchDX)=unique(meta.samples$DX)
colPLEX=viridis(length(unique(meta.samples$PLEX))); names(colPLEX)=unique(meta.samples$PLEX)

###################################
#### save the raw data
###################################
save(meta.samples,meta.pep,raw.pep,samples,pools,plexes,plex.border,plex.mid,pchDX,colPLEX,file=paste0("WorkData/RAW-",tissue,".RData"))


pdf(paste0("Plots/total-intensity-RAW-",tissue,".pdf"),height=8,width=12)
plot(colSums(raw.pep,na.rm=T)/1e6,col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
title(ylab="total~intensity (millions)",line = 2.5)
title(main="Total intensity for samples and pools by plex")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
legend("bottomleft",legend = c("ASD","CONTROL","POOL"),pch=pchDX,col="grey75",bty="n")
dev.off()

pdf(paste0("Plots/fraction-called-RAW-",tissue,".pdf"),height=8,width=12)
plot(colSums(!is.na(raw.pep))/nrow(raw.pep),col=colPLEX[meta.samples$PLEX],pch=pchDX[meta.samples$DX],
     las=1,xaxt="n",xlab="",ylab="")
abline(v=plex.border,lty=3,col="grey50")
axis(side = 1, at=plex.mid, labels = paste("Plex", 1:5))
title(ylab="Fraction",line = 3)
title(main="Fraction peptides with intensities")
legend("bottomleft",legend = c("ASD","CONTROL","POOL"),pch=pchDX,col="grey75",bty="n")
dev.off()

pdf(paste0("Plots/MDS-RAW-",tissue,".pdf"),height = 8, width =8)
plotMDS(log2(raw.pep),col=colPLEX[meta.samples$PLEX],las=1,main="MDS of raw log2(intensity) data",dim.plot = c(1,2))
legend("topleft",legend=paste("Plex",1:5),pch=15,pt.cex = 2,col=colPLEX,bty="n")
dev.off()


