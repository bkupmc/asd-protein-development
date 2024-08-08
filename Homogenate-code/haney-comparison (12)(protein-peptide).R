rm(list=ls()); gc()
options(stringAsFactors=F)

require(read.table)
require(readxl)


### read the Haney et all data transcript level
haney=as.data.frame(read_excel("/data3/DownLoadedData/MacDonald/ASD_Development_June2022/Haney et al 2020/SupplementaryTable3.xlsx",sheet = 2)) # 99,819
ba17=haney[,c(1,2,3,grep("ASD_BA17",colnames(haney)))]  # 23,832 genes; 99,819 transcripts
# take the smallest q value per gene
ba17=ba17[order(ba17$ensembl_gene_id,ba17$ASD_BA17_FDR),]
ba17=ba17[!duplicated(ba17$ensembl_gene_id),] # 23,832



##### read our data peptide level
asd=as.data.frame(read_excel("Results/HOMOGENATE-PEPTIDE.xlsx"),sheet="PTT.DX")  # 35,103
# omit results from peptides mapping to multiple proteins
asd=asd[-grep(";",asd$ensg),]    # 32,882
# omit results for peptides mapping to multiple genes
asd=asd[-grep("/",asd$ensg),]    # 32,829
# select best results for each gene
asd=asd[order(asd$ensg,asd$q),]
asd=asd[!duplicated(asd$ensg),]  # 4700

all.genes=intersect(asd$ensg,ba17$ensembl_gene_id)  # 4610
ba17.genes=ba17$ensembl_gene_id[ba17$ASD_BA17_FDR < 0.05] # 1059
ba17.genes=ba17.genes[ba17.genes %in% all.genes] # 390

# q < 0.05
asd.genes=asd$ensg[asd$q < 0.05]                  # 63
asd.genes=asd.genes[asd.genes %in% all.genes]     # 62

ftable(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#       FALSE TRUE
#
#FALSE   4171  377
#TRUE      49   13

fisher.test(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#Fisher's Exact Test for Count Data
#
#data:  all.genes %in% asd.genes and all.genes %in% ba17.genes
#p-value = 0.001667
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.446846 5.553346
#sample estimates:
#odds ratio 
#  2.933976

#### perform the sign test
top.genes=intersect(ba17.genes,asd.genes)

top.asd=asd[asd$ensg %in% top.genes,]
top.ba17=ba17[ba17$ensembl_gene_id %in% top.genes,]

rownames(top.asd)=top.asd$ensg
rownames(top.ba17)=top.ba17$ensembl_gene_id

sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"]))  # 12
pbinom(sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"])),length(top.genes),p=.5,lower.tail = F) # 0.00012

####### 0.10
all.genes=intersect(asd$ensg,ba17$ensembl_gene_id)  # 4610
ba17.genes=ba17$ensembl_gene_id[ba17$ASD_BA17_FDR < 0.10] # 1059
ba17.genes=ba17.genes[ba17.genes %in% all.genes] # 644

# q < 0.10
asd.genes=asd$ensg[asd$q < 0.10]                  # 248
asd.genes=asd.genes[asd.genes %in% all.genes]     # 243

ftable(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#       FALSE TRUE
#
#FALSE   3792  575
#TRUE     174   69

fisher.test(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#Fisher's Exact Test for Count Data
#
#data:  all.genes %in% asd.genes and all.genes %in% ba17.genes
#p-value = 1.554e-09
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.921918 3.526143
#sample estimates:
#odds ratio 
# 2.614609

#### perform the sign test
top.genes=intersect(ba17.genes,asd.genes)

top.asd=asd[asd$ensg %in% top.genes,]
top.ba17=ba17[ba17$ensembl_gene_id %in% top.genes,]

rownames(top.asd)=top.asd$ensg
rownames(top.ba17)=top.ba17$ensembl_gene_id

sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"]))  # 61
pbinom(sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"])),length(top.genes),p=.5,lower.tail = F) # 2.0514e-12


ba17=ba17[ba17$ensembl_gene_id %in% all.genes,]
rownames(ba17)=ba17$ensembl_gene_id
asd=asd[asd$ensg %in% all.genes,]
rownames(asd)=asd$ensg

peptide.transcript=data.frame(ensg=all.genes)
peptide.transcript=cbind.data.frame(peptide.transcript,asd[peptide.transcript$ensg,],ba17[peptide.transcript$ensg,])
peptide.transcript[,c(4,5,6,7,8,9,12,13,15)]=NULL
colnames(peptide.transcript)=c("ensg","peptide","effectDI","q.DI","protein","transcript","gene","ASD_BA17_logFC.DE","ASD_BA17_FDR.DE")
peptide.transcript=peptide.transcript[,c("transcript","ensg","gene","peptide","protein","effectDI","q.DI","ASD_BA17_logFC.DE","ASD_BA17_FDR.DE")]

WriteXLS(peptide.transcript,"Results/HOMOGENATE-PEPTIDE-DI-vs-TRANSCRIPT-DE-ASD_10-19-2022.xlsx")






# now using the protein/gene results
haney=as.data.frame(read_excel("/data3/DownLoadedData/MacDonald/ASD_Development_June2022/Haney et al 2020/SupplementaryTable3.xlsx",sheet = 1)) # 24,836
ba17=haney[,c(1,2,3,grep("BA17",colnames(haney)))]  # 24,836 genes
# take the smallest q value per gene
ba17=ba17[order(ba17$ensembl_gene_id,ba17$ASD_BA17_FDR),]
ba17=ba17[!duplicated(ba17$ensembl_gene_id),] # 24,836



##### read our data peptide level
asd=as.data.frame(read_excel("Results/HOMOGENATE-PROTEIN.xlsx"),sheet="PTT.DX")  # 4727
# omit results for proteins mapping to multiple genes
asd=asd[-grep("/",asd$ensg),]    # 4715
# select best results for each gene, this removes proteins mapping to the same gene
asd=asd[order(asd$ensg,asd$q),]
asd=asd[!duplicated(asd$ensg),]  # 4700

all.genes=intersect(asd$ensg,ba17$ensembl_gene_id)  # 4618
ba17.genes=ba17$ensembl_gene_id[ba17$ASD_BA17_FDR <= 0.05] # 3264
ba17.genes=ba17.genes[ba17.genes %in% all.genes] # 1087

# q < 0.05
asd.genes=asd$ensg[asd$q <= 0.05]                  # 67
asd.genes=asd.genes[asd.genes %in% all.genes]     # 64



ftable(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#FALSE TRUE
#
#FALSE   3488  1066
#TRUE      43    21

fisher.test(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#Fisher's Exact Test for Count Data
#
#data:  all.genes %in% asd.genes and all.genes %in% ba17.genes
#p-value = 0.101
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.8964558 2.7663621
#sample estimates:
#odds ratio 
#  1.597794 

#### perform the sign test
top.genes=intersect(ba17.genes,asd.genes)

top.asd=asd[asd$ensg %in% top.genes,]
top.ba17=ba17[ba17$ensembl_gene_id %in% top.genes,]

rownames(top.asd)=top.asd$ensg
rownames(top.ba17)=top.ba17$ensembl_gene_id

sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"]))  # 21
pbinom(sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"])),length(top.genes),p=.5,lower.tail = F) # 0


####### 0.10
all.genes=intersect(asd$ensg,ba17$ensembl_gene_id)  # 4618
ba17.genes=ba17$ensembl_gene_id[ba17$ASD_BA17_FDR <= 0.10] # 3264
ba17.genes=ba17.genes[ba17.genes %in% all.genes] # 1478

# q < 0.10
asd.genes=asd$ensg[asd$q <= 0.10]                  # 237
asd.genes=asd.genes[asd.genes %in% all.genes]     # 218



ftable(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#FALSE TRUE
#
#FALSE   3013  1387
#TRUE      127  91

fisher.test(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#Fisher's Exact Test for Count Data
#
#data:  all.genes %in% asd.genes and all.genes %in% ba17.genes
#p-value = 0.002217
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.166672 2.069707
#sample estimates:
#odds ratio 
#  1.556401

#### perform the sign test
top.genes=intersect(ba17.genes,asd.genes)

top.asd=asd[asd$ensg %in% top.genes,]
top.ba17=ba17[ba17$ensembl_gene_id %in% top.genes,]

rownames(top.asd)=top.asd$ensg
rownames(top.ba17)=top.ba17$ensembl_gene_id

sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"]))  # 21
pbinom(sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"])),length(top.genes),p=.5,lower.tail = F) # 0

#### create a dataset
ba17=ba17[ba17$ensembl_gene_id %in% all.genes,]
rownames(ba17)=ba17$ensembl_gene_id
asd=asd[asd$ensg %in% all.genes,]
rownames(asd)=asd$ensg

protein.gene=data.frame(ensg=all.genes)
protein.gene=cbind.data.frame(protein.gene,asd[protein.gene$ensg,],ba17[protein.gene$ensg,])
protein.gene[,c(4,5,6,7,8,9,12,13,14,18,19)]=NULL
colnames(protein.gene)=c("ensg","protein","effectDI","q.DI","gene","gene_biotype","ASD_BA17_logFC.DE","ASD_BA17_FDR.DE")
protein.gene=protein.gene[,c("ensg","gene","gene_biotype","protein","effectDI","q.DI","ASD_BA17_logFC.DE","ASD_BA17_FDR.DE")]

WriteXLS(protein.gene,"Results/HOMOGENATE-PROTEIN-DI-vs-GENE-DE-ASD_10-19-2022.xlsx")