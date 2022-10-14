library(viridis)
library(DESeq2)
library(ggplot2)
library(limma) 
library(RColorBrewer)
library(pheatmap)
library(grid)
library(limma)
library(topGO)
library(biomaRt)

# txi object created using tximport package of kallisto abundance files
sampleLable = colnames(txi$abundance)
subtype = sub("_[1-9]$", "", sampleLable)

sampInfo = subtype
sampInfo = factor(sampInfo)

sampInfo = data.frame(group=sampInfo)
rownames(sampInfo) = sampleLable

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sampInfo,
                                   design = ~ group)

dds <- estimateSizeFactors(ddsTxi)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)
ntd_R1KD = normTransform(dds)
normExpression_R1KD =  assay(ntd_R1KD[,1:12])

mads = apply(normExpression, 1, mad)
qq = mads[order(-mads)][1:100]

A375_R1_KD = data.frame(results(dds, contrast=c("group","A375_KD", "A375_NT")))
a451_R1_KD = data.frame(results(dds, contrast=c("group","451_KD", "451_NT")))

A375_SRCi = data.frame(results(dds, contrast=c("group","A375_E_SRCi", "A375_E_VEH")))

sigA375 = A375_R1_KD[A375_R1_KD[,"padj"]<0.01,]
sigA375 = sigA375[!is.na(sigA375[,2]),]
sigA375[order(-sigA375[,"log2FoldChange"]),]

sig451 = a451_R1_KD[a451_R1_KD[,"padj"]<0.01,]
sig451 = sig451[!is.na(sig451[,2]),]
sig451[order(sig451[,"log2FoldChange"]),]

A375_SRCi = A375_SRCi[A375_SRCi[,"padj"]<0.01,]
A375_SRCi = A375_SRCi[!is.na(A375_SRCi[,"padj"]),]

write.csv(sigA375, file="a375_rac1KD.csv")
write.csv(sig451, file="451lu_rac1KD.csv")
write.csv(A375_SRCi, file="a375_src.csv")

a375_ranks = sigA375[,"log2FoldChange"]
names(a375_ranks) = rownames(sigA375)
a375_ranks = a375_ranks[order(a375_ranks)]

## invasive_genes2 are genes that define the Undifferentaited 
inv = read.csv("invasive_genes2.csv", stringsAsFactors = FALSE, header=FALSE)
invList = list(invasive=unlist(inv))

fgseaRes <- fgsea(pathways = invList,
                  stats    = a375_ranks)

fgseaRes[order(fgseaRes[,"NES"]),]
fgseaRes[order(-fgseaRes[,"NES"]),]

pdf(file="A375_enrich.pdf", width=3.5, height=3)
plotEnrichment(invList$invasive, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25))
dev.off()

## a375 enrichment plot
a375_ranks = sigA375[,"log2FoldChange"]
names(a375_ranks) = rownames(sigA375)
# a375_ranks = a375_ranks[order(a375_ranks)]
a375_ranks = a375_ranks[abs(a375_ranks)>0.5]

inv = read.csv("invasive_genes2.csv", stringsAsFactors = FALSE, header=FALSE)
invList = list(invasive=unlist(inv))

fgseaRes <- fgsea(pathways = invList,
                  stats    = a375_ranks)

fgseaRes <- fgsea(pathways = pathways.hallmark ,
                  stats = a375_ranks)

fgseaRes[order(fgseaRes[,"NES"]),]
fgseaRes[order(-fgseaRes[,"NES"]),]

fgseaResFilt = fgseaRes[which(abs(fgseaRes$padj) <0.05), ] 
fgseaResFilt = data.frame(pathway=fgseaResFilt$pathway, Padj=fgseaResFilt$padj, NES=fgseaResFilt$NES)

pdf("hallmark_pathways_A375.pdf", width=6, height=3)
ggplot(fgseaResFilt, aes(reorder(pathway, -NES), NES)) + geom_col() + 
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
  theme_minimal()
dev.off()


## 451Lu enrichment plot
a375_ranks = sig451[,"log2FoldChange"]
names(a375_ranks) = rownames(sig451)
# a375_ranks = a375_ranks[order(a375_ranks)]
a375_ranks = a375_ranks[abs(a375_ranks)>0.5]

fgseaRes <- fgsea(pathways = invList,
                  stats    = a375_ranks)

# fgseaRes <- fgsea(pathways = pathways.hallmark ,
#                   stats = a375_ranks[names(a375_ranks) %in% sigCommon])

fgseaRes <- fgsea(pathways = pathways.hallmark ,
                  stats = a375_ranks)

fgseaRes[order(fgseaRes[,"NES"]),]
fgseaRes[order(-fgseaRes[,"NES"]),]

fgseaResFilt = fgseaRes[which(abs(fgseaRes$padj) <0.05), ] 
fgseaResFilt = data.frame(pathway=fgseaResFilt$pathway, Padj=fgseaResFilt$padj, NES=fgseaResFilt$NES)

pdf("hallmark_pathways_451Lu.pdf", width=6, height=3)
  ggplot(fgseaResFilt, aes(reorder(pathway, -NES), NES)) + geom_col() + 
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score") + 
    theme_minimal()
dev.off()

fgseaResFilt[order(fgseaResFilt[,"NES"]),]
fgseaResFilt[order(-fgseaResFilt[,"NES"]),]

sigA375[order(sigA375[,2]),]

a3kdg = rownames(sigA375[sigA375[,2]< -2,])
l45kdg = rownames(sig451[sig451[,2]< -2,])

commonGenes = intersect(a3kdg, l45kdg)


pdf("RAC1_KD_genes_hm.pdf", width=6, height=3.2)
  pheatmap(scaled[c( "RAC1","AXL", "NGFR", "EGFR", "PDGFRB","CAV1", "ANKRD1" , "CTGF", "CYR61", "TGFB2"), 
                  c("A375_NT_1", "A375_NT_2", "A375_NT_3", "451_NT_1", "451_NT_2", "451_NT_3",
                    "A375_KD_1", "A375_KD_2", "A375_KD_3", "451_KD_1", "451_KD_2", "451_KD_3"
                  )],
           
           treeheight_row = 0, treeheight_col = 0, fontsize = 12, cluster_cols = FALSE, cluster_rows = FALSE, color = magma(10))
dev.off()

##
pdf(file="451Lu_enrich.pdf", width=3.5, height=3)
plotEnrichment(invList$invasive, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12), legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15))
dev.off()

qq = invList$invasive[invList$invasive %in% rownames(normExpression_R1KD)]

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

scaled = t(as.data.frame(apply(normExpression_R1KD, 1, normalize)))



pheatmap(scaled[fgseaRes$leadingEdge[[1]], c("451_NT_1", "451_NT_2", "451_NT_3",
                                             "451_KD_1", "451_KD_2", "451_KD_3")])

############### Src analysis
a375_ranks = A375_SRCi[,"log2FoldChange"]
names(a375_ranks) = rownames(A375_SRCi)
# a375_ranks = a375_ranks[order(a375_ranks)]
a375_ranks = a375_ranks[abs(a375_ranks) > 0.5]

fgseaRes <- fgsea(pathways = invList,
stats    = a375_ranks)

fgseaRes <- fgsea(pathways = list(Rac1Genes = commonGenes),
stats    = a375_ranks)

pdf(file="Src_inv_enrich.pdf", width=7, height=5)
plotEnrichment(invList$invasive, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25))
dev.off()

pdf(file="Src_rac1_enrich.pdf", width=7, height=5)
plotEnrichment(commonGenes, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25))
dev.off()

fgseaRes <- fgsea(pathways = pathways.hallmark ,
                  stats = a375_ranks)

fgseaRes <- fgsea(pathways = list(rac1=commonGenes), stats = a375_ranks)

plotEnrichment(pathways.hallmark$HALLMARK_TGF_BETA_SIGNALING, a375_ranks, gseaParam = 1, ticksSize = 0.3) +  
  theme(axis.title = element_text(size = 25), 
        axis.text = element_text(size = 15), legend.text = element_text(size = 20), 
        legend.title = element_text(size = 25))

fgseaRes[order(fgseaRes[,"NES"]),]
fgseaRes[order(-fgseaRes[,"NES"]),]

fgseaResFilt = fgseaRes[which(abs(fgseaRes$padj) <0.05), ] 
fgseaResFilt = data.frame(pathway=fgseaResFilt$pathway, Padj=fgseaResFilt$padj, NES=fgseaResFilt$NES)

pdf(file="hallmarks_SRCi.pdf", width=3.5, height=3)
  ggplot(fgseaResFilt, aes(reorder(pathway, -NES), NES)) + geom_col() + 
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score") + 
    theme_minimal()

  dev.off()

##
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

hgncGenes = getBM(filters= "hgnc_symbol", attributes= c( "entrezgene_id", "hgnc_symbol"), 
                  values=sigCommon, mart=human)

temp = goana(hgncGenes[,1], FDR=0.01)
temp = temp[temp[,2]=="BP", ]
temp = temp[temp[,3]<100, ]
temp[order(temp[,5]),]

erPathways = topGO(temp, "BP")
erPathways = erPathways[order(erPathways[,"P.DE"]),]

head(erPathways[order(erPathways[,"P.DE"]),],100)

fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats = a375_ranks)

library("ggVennDiagram")
library("ggvenn")

sigA3752 = sigA375[abs(sigA375[,"log2FoldChange"]) > 0.5,]

A375_SRCi2 = A375_SRCi[ abs(A375_SRCi[,"log2FoldChange"]) > 0.5, ]

qq = list(RAC1_KD=rownames(sigA3752), SRCi=rownames(A375_SRCi2))

pdf(file="venn.pdf", width = 7, height = 4)
  ggvenn(qq, columns = c("RAC1_KD", "SRCi"), fill_color = c("#EFC000FF", "#868686FF") )
dev.off()

