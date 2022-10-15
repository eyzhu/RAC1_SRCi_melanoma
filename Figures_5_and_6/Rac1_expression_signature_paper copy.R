library(viridis)
library(enrichR)
library(topGo)
library(Boruta)
library(DESeq2)
library(ggplot2)
library(limma) 
library(RColorBrewer)
library(pheatmap)
library(grid)
library(fmsb)
library(biomaRt)
library(singscore)
library(KEGGREST)
library(TCGAbiolinks)
library(ggfortify)
library(stringr)

query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")

mydata <- GDCprepare(query)
skcmMatrix <- assay(mydata)

skcmMatrix = skcmMatFinal
skcmMatFilt = skcmMatrix[rowSums(skcmMatrix>10)>10, ]

wantGenes = rownames(skcmMatFilt)

ensembl=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")  

hgncGenes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                             "hgnc_symbol"), values=wantGenes, mart= ensembl)

swapGenes = hgncGenes[match(wantGenes, hgncGenes[,1]) , 2]
badGeneInd = which(swapGenes=="")

rownames(skcmMatFinal) = swapGenes
