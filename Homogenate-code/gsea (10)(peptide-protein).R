rm(list=ls()); gc()
options(stringsAsFactors = F)

# this part performs Gene Set Enrichment Analysis


tissue="HOMOGENATE"
#### libraries used
require(data.table)           # fast version of read.table
require(readxl)               # used for reading Excel format files
require(viridis)              # colorblind friendly color palette
require(matrixStats)          # to get statistics on matrices
require(doParallel)           # parallel computing
require(limma)
require(sva)


### load the results data
load(paste0("WorkData/RESULTS-",tissue,".RData"))


#### perform GSEA
# load the Msig database
msig=readRDS("/data3/DownLoadedData/GSEA-GeneSetEnrichment/7.5.1/MsigDB.v7.5.1_08-15-2022.rds")

# functions needed
unifiedGeneNames <- function(gene.names){
  # get the alias
  gene.alias=alias2SymbolTable(gene.names)
  # set the name for which their is no alias to the original gene name
  gene.alias[is.na(gene.alias)]=gene.names[is.na(gene.alias)]
  # replace gene.alias with original gene.names for the gene.alias that appears more than once
  i.dup=which(gene.alias %in% gene.alias[duplicated(gene.alias)])
  gene.alias[i.dup]=gene.names[i.dup]
  
  # return the list with updated names
  return(gene.alias)
}

geneSetEnrichment<-function(gs.term,genes=genes,top.genes=top.genes,gene.xref=gene.xref,GS){
  n.gene=length(genes)
  n.top=length(top.genes)
  # unify the gene names
  gs.genes=unifiedGeneNames(GS[[gs.term]])
  # determine the OR and p-value using Fisher exact test
  n.gs=sum(genes %in% gs.genes)
  n.gs.top=sum(top.genes %in% gs.genes)
  re<-fisher.test(matrix(c(n.gs.top,n.gs-n.gs.top,n.top-n.gs.top,n.gene-n.gs-n.top+n.gs.top),nrow=2,byrow = T))
  included.genes=top.genes[top.genes %in% gs.genes]
  included.genes=gene.xref[included.genes,"de.gene"]
  results=data.frame(GSterm=gs.term,GSgenes=n.gs,DEgenes=n.gs.top,expected=(n.gs/n.gene)*n.top,OR=re$estimate,p=re$p.value,genes=paste(included.genes,collapse=","))
  return(results)
}



geneSetEnrichmentAnalysis <- function(genes,sel.genes,pathways,msig){
  # function to perform the gene set enrichment given a set of genes
  # use the unified gene names in limma
  uni.genes=unifiedGeneNames(genes)
  # set up a cross-reference 
  gene.xref=data.frame(de.gene=genes,uni.gene=uni.genes)
  row.names(gene.xref)=gene.xref$uni.gene
  
  # go through each of the requested pathways and find the significant ones
  RE.sheet=NULL
  for(i.path in 1:nrow(pathways)){
    RE=rbindlist(mclapply(names(msig[[pathways[i.path,"cx"]]][[pathways[i.path,"set"]]]$gs),geneSetEnrichment,genes=genes,top.genes=sel.genes,
                        GS=msig[[pathways[i.path,"cx"]]][[pathways[i.path,"set"]]]$gs,gene.xref=gene.xref,mc.cores=16,mc.preschedule = T))
    RE=data.frame(RE)
    # apply fdr correction
    RE$adjP=p.adjust(RE$p,method="fdr")
    # sort by adjusted p-value
    RE=RE[order(RE$adjP,RE$p),]
    # select the pathways passing fdr limit
    RE=RE[RE$adjP < 0.05 & RE$OR > 1,]
    print(paste(pathways$set.name[i.path],nrow(RE))); flush.console()
    if(nrow(RE) > 0){
      RE$database=pathways$set.name[i.path]
      RE=RE[,c("database","GSterm","GSgenes","DEgenes","expected","OR","p","adjP","genes")]
      RE.sheet=rbind.data.frame(RE.sheet,RE)
    }
  }
  if(is.null(RE.sheet)){RE.sheet=data.frame(n=0)}
  return(RE.sheet)

}

# set up a set of pathways you want to test
# you can find the available ones using 
# names(msig)    # names of the major sets
# names(msig[["c2"]])  # names of the subsets within a set
# msig[["c2"]][["kegg"]]$description  # description of the subset
# msig[["c2"]][["kegg"]]$n  # number of pathways in the subset
# these are the common ones we use
test.pathways=data.frame(cx=c("c2","c2","c2","c2","c2","c5","c5","c5","c7","h"),
                         set=c("cgp","biocarta","kegg","pid","reactome","bp","cc","mf","immune","all"),
                         set.name=c("curated; chemical and genetic perturbations","curated; biocarta","curated; kegg","curated; pid","curated; reactome",
                                    "GO; biological processes","GO; cellular components","GO; molecular functions",
                                    "immunologic signatures",
                                    "hallmark gene sets"))

#test.pathways=test.pathways[c(2,7),] # this you can use to select a subset


# PEPTIDES
GSEA.PEPTIDE=list()
for(set in names(RE.PEPTIDE)){
  print(set); flush.console()
  # get the names for all the genes in the analysis
  genes=unique(unlist(strsplit(unlist(strsplit(RE.PEPTIDE[[set]]$gene,"/")),"; ")))
  genes=genes[genes != ""]
  # vector of selected genes
  if(length(RE.PEPTIDE[[set]]$gene[RE.PEPTIDE[[set]]$q < 0.05]) == 0){
    sel.genes=NULL
  }else{
    sel.genes=unique(unlist(strsplit(unlist(strsplit(RE.PEPTIDE[[set]]$gene[RE.PEPTIDE[[set]]$q < 0.05],"/")),"; ")))
    sel.genes=sel.genes[sel.genes != ""]
  }
  print(length(sel.genes)); flush.console()
  # run the GSEA
  if(length(sel.genes) > 0){
    GSEA.PEPTIDE[[set]]=geneSetEnrichmentAnalysis(genes,sel.genes=sel.genes,pathways=test.pathways,msig=msig)
 }else{
    GSEA.PEPTIDE[[set]]=data.frame(n=0)
  }
}

WriteXLS(GSEA.PEPTIDE,paste0("Results/GSEA",tissue,"-PEPTIDE.xlsx"))

# repeat for the PROTEINS
GSEA.PROTEIN=list()
for(set in names(RE.PROTEIN)){
  print(set); flush.console()
  # get the names for all the genes in the analysis
  genes=unique(unlist(strsplit(unlist(strsplit(RE.PROTEIN[[set]]$gene,"/")),"; ")))
  genes=genes[genes != ""]
  # vector of selected genes
  if(length(RE.PROTEIN[[set]]$gene[RE.PROTEIN[[set]]$q < 0.05]) == 0){
    sel.genes=NULL
  }else{
    sel.genes=unique(unlist(strsplit(unlist(strsplit(RE.PROTEIN[[set]]$gene[RE.PROTEIN[[set]]$q < 0.05],"/")),"; ")))
    sel.genes=sel.genes[sel.genes != ""]
  }
  print(length(sel.genes)); flush.console()
  # run GSEA
  if(length(sel.genes) > 0){
    GSEA.PROTEIN[[set]]=geneSetEnrichmentAnalysis(genes,sel.genes=sel.genes,pathways=test.pathways,msig=msig)
  }else{
    GSEA.PROTEIN[[set]]=data.frame(n=0)
  }
}
WriteXLS(GSEA.PROTEIN,paste0("Results/GSEA",tissue,"-PROTEIN.xlsx"))


save(GSEA.PEPTIDE,GSEA.PROTEIN,file="WorkData/GSEA-HOMOGENATE.RData")

