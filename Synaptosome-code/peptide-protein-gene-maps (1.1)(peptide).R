rm(list=ls()); gc()
options(stringsAsFactors = F)

tissue="SYNAPTOSOME-PEPTIDE"
# script will create the annotations that might be needed later
library(doParallel)

load(paste0("WorkData/RAW-",tissue,".RData"))
colnames(meta.pep)=c("peptide","protein")

library(biomaRt)
mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl")


#### add the gene information to the peptide and protein information
# find gene names for the proteins on biomaRt
proteins=unique(gsub(" ","",unlist(strsplit(meta.pep$protein[!is.na(meta.pep$protein)],";"))))  # 7577

# find the ensembl and external gene id linked to the protein
reBM=getBM(c("uniprotswissprot","ensembl_gene_id","external_gene_name","chromosome_name"),filters="uniprotswissprot",
           values=proteins,mart=mart) # 8164
reBM=reBM[reBM$chromosome_name %in% c(1:23,"X","Y","MT"),]  # 7519

length(unique(reBM$uniprotswissprot))  # 7491
length(unique(reBM$ensembl_gene_id))    # 7516
length(unique(reBM$external_gene_name)) # 7513


# process them to allow for multiple genes for a protein
geneByProtein <- function(protein,reBM){
  i=which(reBM$uniprotswissprot == protein)
  genes=reBM$external_gene_name[i]; genes=genes[!is.na(genes)]
  genes=paste(genes,collapse="/")
  ensgs=reBM$ensembl_gene_id[i]; ensgs=ensgs[!is.na(ensgs)]
  ensgs=paste(ensgs,collapse="/")
  return(data.frame(protein,genes,ensgs))
}

re.protein=rbindlist(lapply(proteins,geneByProtein,reBM=reBM))
re.protein=as.data.frame(re.protein)
rownames(re.protein)=re.protein$protein

# prcess to allow for multiple proteins for a peptide
proteinsByPeptide <- function(peptide,meta.pep,re.protein){
  protein=unlist(strsplit(meta.pep[peptide,"protein"],"; "))
  gene=paste(re.protein[protein,"genes"],collapse = "; ")
  ensg=paste(re.protein[protein,"ensgs"],collapse = "; ")
  protein=paste(protein,collapse="; ")
  return(data.frame(peptide,protein,gene,ensg))
}

re.peptide=as.data.frame(rbindlist(mclapply(meta.pep$peptide,proteinsByPeptide,meta.pep=meta.pep,re.protein=re.protein,mc.cores=32,mc.preschedule = T)))
rownames(re.peptide)=re.peptide$peptide


protein.map=re.protein
peptide.map=re.peptide
colnames(protein.map)=c("protein","gene","ensg")
save(protein.map,peptide.map,file="WorkData/peptide-protein-gene-maps.RData")
