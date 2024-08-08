rm(list=ls()); gc()
options(stringsAsFactors = F)

require(data.table)
require(WriteXLS)
require(doParallel)
require(viridis)

### load the peptide-protein-gene mapping information
load("WorkData/peptide-protein-gene-maps.RData")

### read the results information for the peptides
tissue="HOMOGENATE-PEPTIDE"
load(paste0("WorkData/PAIRED-T-TEST-",tissue,".RData"))

peptide.PTT=re.paired.t.test
peptide.PTT$protein=peptide.map[peptide.PTT$peptide,"protein"]
peptide.PTT$gene=peptide.map[peptide.PTT$peptide,"gene"]
peptide.PTT$ensg=peptide.map[peptide.PTT$peptide,"ensg"]

load(paste0("WorkData/LMER-REGRESSION-",tissue,".RData"))

peptide.MLM.DX=reDX
peptide.MLM.DX$protein=peptide.map[peptide.MLM.DX$peptide,"protein"]
peptide.MLM.DX$gene=peptide.map[peptide.MLM.DX$peptide,"gene"]
peptide.MLM.DX$ensg=peptide.map[peptide.MLM.DX$peptide,"ensg"]

peptide.MLM.AGE=reAGE
peptide.MLM.AGE$protein=peptide.map[peptide.MLM.AGE$peptide,"protein"]
peptide.MLM.AGE$gene=peptide.map[peptide.MLM.AGE$peptide,"gene"]
peptide.MLM.AGE$ensg=peptide.map[peptide.MLM.AGE$peptide,"ensg"]

RE.PEPTIDE=list()
RE.PEPTIDE[["PTT.DX"]]=peptide.PTT
RE.PEPTIDE[["MLM.DX"]]=peptide.MLM.DX
RE.PEPTIDE[["MLM.AGE"]]=peptide.MLM.AGE

WriteXLS(RE.PEPTIDE,paste0("Results/",tissue,".xlsx"))

### read the results information for the proteins
tissue="HOMOGENATE-PROTEIN"
load(paste0("WorkData/PAIRED-T-TEST-",tissue,".RData"))

protein.PTT=re.paired.t.test
protein.PTT$gene=protein.map[protein.PTT$protein,"gene"]
protein.PTT$ensg=protein.map[protein.PTT$protein,"ensg"]

load(paste0("WorkData/LMER-REGRESSION-",tissue,".RData"))

protein.MLM.DX=reDX
protein.MLM.DX$gene=protein.map[protein.MLM.DX$protein,"gene"]
protein.MLM.DX$ensg=protein.map[protein.MLM.DX$protein,"ensg"]

protein.MLM.AGE=reAGE
protein.MLM.AGE$gene=protein.map[protein.MLM.AGE$protein,"gene"]
protein.MLM.AGE$ensg=protein.map[protein.MLM.AGE$protein,"ensg"]

RE.PROTEIN=list()
RE.PROTEIN[["PTT.DX"]]=protein.PTT
RE.PROTEIN[["MLM.DX"]]=protein.MLM.DX
RE.PROTEIN[["MLM.AGE"]]=protein.MLM.AGE

WriteXLS(RE.PROTEIN,paste0("Results/",tissue,".xlsx"))

save(norm.pep,norm.prt,meta.samples,RE.PEPTIDE,RE.PROTEIN,
     colPLEX,pchDX,plex.border,plex.mid,plexes,pools,samples,selected.model,file="WorkData/RESULTS-HOMOGENATE.RData")




### make a plot compare results from PTT and MLM.DX
pdf("Plots/comparing-PTT-MLM-DX-HOMOGENATE.pdf",height=12,width=12)
par(mfrow=c(2,2))
par(mar=c(3,3,3,1))
par(oma=c(3,3,1,1))
plot(peptide.PTT$effect,peptide.MLM.DX$Estimate,las=1,pch=19,col=viridis(5)[2])
mtext("Peptide",side=2,line=3)
mtext("Fold Change",side=3,line=1)
abline(a=0,b=1,col="grey25")

plot(-log10(peptide.PTT$p.emperical),-log10(peptide.MLM.DX$p),las=1,pch=19,col=viridis(5)[2],xlim=range(-log10(peptide.PTT$p.emperical)),ylim=range(-log10(peptide.PTT$p.emperical)))
mtext("-log10(P)",side=3,line=1)
abline(a=0,b=1,col="grey25")

plot(protein.PTT$effect,protein.MLM.DX$Estimate,las=1,pch=19,col=viridis(5)[3])
mtext("Protein",side=2,line=3)
abline(a=0,b=1,col="grey25")

plot(-log10(protein.PTT$p.emperical),-log10(protein.MLM.DX$p),las=1,pch=19,col=viridis(5)[3],xlim=range(-log10(protein.PTT$p.emperical)),ylim=range(-log10(protein.PTT$p.emperical)),
    xlab="paired t-test",ylab="")
abline(a=0,b=1,col="grey25")

mtext("Paired T-Test",side=1,line=1,outer=T)
mtext("Mixed Linear Model",side=2,line=1,outer=T)
dev.off()


