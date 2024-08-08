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
asd=as.data.frame(read_excel("Results/SYNAPTOSOME-PEPTIDE.xlsx"),sheet="PTT.DX")  # 29,311
# omit results from peptides mapping to multiple proteins
asd=asd[-grep(";",asd$ensg),]    # 27,491
# omit results for peptides mapping to multiple genes
asd=asd[-grep("/",asd$ensg),]    # 27,439
# select best results for each gene
asd=asd[order(asd$ensg,asd$q),]
asd=asd[!duplicated(asd$ensg),]  # 4260

all.genes=intersect(asd$ensg,ba17$ensembl_gene_id)  # 4188
ba17.genes=ba17$ensembl_gene_id[ba17$ASD_BA17_FDR < 0.05] # 1059
ba17.genes=ba17.genes[ba17.genes %in% all.genes] # 347

# q < 0.05
asd.genes=asd$ensg[asd$q < 0.05]                  # 2
asd.genes=asd.genes[asd.genes %in% all.genes]     # 2

ftable(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#       FALSE TRUE
#
#FALSE   3839  347
#TRUE       2    0


##### NO NEED TO DO THIS
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










# now using the protein/gene results
haney=as.data.frame(read_excel("/data3/DownLoadedData/MacDonald/ASD_Development_June2022/Haney et al 2020/SupplementaryTable3.xlsx",sheet = 1)) # 24,836
ba17=haney[,c(1,2,3,grep("BA17",colnames(haney)))]  # 24,836 genes
# take the smallest q value per gene
ba17=ba17[order(ba17$ensembl_gene_id,ba17$ASD_BA17_FDR),]
ba17=ba17[!duplicated(ba17$ensembl_gene_id),] # 24,836



##### read our data peptide level
asd=as.data.frame(read_excel("Results/SYNAPTOSOME-PROTEIN.xlsx"),sheet="PTT.DX")  # 4287
# omit results for proteins mapping to multiple genes
asd=asd[-grep("/",asd$ensg),]    # 4275
# select best results for each gene, this removes proteins mapping to the same gene
asd=asd[order(asd$ensg,asd$q),]
asd=asd[!duplicated(asd$ensg),]  # 4260

all.genes=intersect(asd$ensg,ba17$ensembl_gene_id)  # 4195
ba17.genes=ba17$ensembl_gene_id[ba17$ASD_BA17_FDR < 0.10] # 3264
ba17.genes=ba17.genes[ba17.genes %in% all.genes] # 1001

# q < 0.05
asd.genes=asd$ensg[asd$q < 0.10]                  # 5
asd.genes=asd.genes[asd.genes %in% all.genes]     # 5



ftable(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#FALSE TRUE
#
#FALSE   3191   999
#TRUE      3      2

fisher.test(all.genes %in% asd.genes,all.genes %in% ba17.genes)

#Fisher's Exact Test for Count Data
#
#data:  all.genes %in% asd.genes and all.genes %in% ba17.genes
#p-value = 0.3432
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.1775807 18.6079039
#sample estimates:
#odds ratio 
#  2.129062

##### NO NEED TO DO THIS
#### perform the sign test
top.genes=intersect(ba17.genes,asd.genes)

top.asd=asd[asd$ensg %in% top.genes,]
top.ba17=ba17[ba17$ensembl_gene_id %in% top.genes,]

rownames(top.asd)=top.asd$ensg
rownames(top.ba17)=top.ba17$ensembl_gene_id

sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"]))  # 21
pbinom(sum(sign(top.asd[top.genes,"effect"]) == sign(top.ba17[top.genes,"ASD_BA17_logFC"])),length(top.genes),p=.5,lower.tail = F) # 0
