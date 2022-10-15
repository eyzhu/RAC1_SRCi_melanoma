library(viridis)
library(pheatmap)
library(M3C)
library(Rtsne)
library(DEseq2)
library(GEOquery)
library(limma)
library(Biobase)
library(e1071)
library(caret)
library(reshape2)

## load Verfaillie et al data
siMITF_sk28 = read.csv("data/Verfaillie/sk28MITF_si_DEGs.csv", header=TRUE, stringsAsFactors = FALSE, skip = 1)
siMITF_sk28_down = unique(siMITF_sk28[siMITF_sk28[,"log2fc"]<0, "gene_symbol"])
siMITF_sk28_up = unique(siMITF_sk28[siMITF_sk28[,"log2fc"]>0, "gene_symbol"])

siMITF_501 = read.csv("data/Verfaillie/501MITF_si_DEGs.csv", header=TRUE, stringsAsFactors = FALSE, skip = 1)
siMITF_501_down = unique(siMITF_501[siMITF_501[,"log2"]<0, "gene_symbol"])
siMITF_501_up = unique(siMITF_501[siMITF_501[,"log2"]>0, "gene_symbol"])

oeMITF_a375 = read.csv("data/Verfaillie/a375MITF_OE_DEGs.csv", header=TRUE, stringsAsFactors = FALSE, skip = 1)
oeMITF_a375_up = unique(oeMITF_a375[oeMITF_a375[,"log2"]>0, "gene_symbol"])
oeMITF_a375_down = unique(oeMITF_a375[oeMITF_a375[,"log2"]<0, "gene_symbol"])

mitfChangeGenes = as.character(c(siMITF_sk28_down,siMITF_501_down,oeMITF_a375_up))
mitfGenesFreq = sort(table(mitfChangeGenes),decreasing = TRUE)

MITFregGenes = names(mitfGenesFreq[mitfGenesFreq>1])

mitfChangeGenesUp = as.character(c(siMITF_sk28_up, siMITF_501_up, oeMITF_a375_down))
mitfGenesFreqUp = sort(table(mitfChangeGenesUp),decreasing = TRUE)

## load MITF chip
sk28CutNrun = unique(read.csv("data/Verfaillie/sk28MITF_cutandrun.csv", stringsAsFactors = FALSE, skip = 1)[,"SYMBOL"])
mel501_chip = unique(read.csv("data/Verfaillie/501MITF_chip.csv", stringsAsFactors = FALSE, skip = 1)[,"SYMBOL"])
colo829chip = unique(read.csv("data/Verfaillie/colo829MITF_chip.csv", stringsAsFactors = FALSE, skip = 1)[,"SYMBOL"])

mitfchipTars = unique(c(sk28CutNrun, mel501_chip, colo829chip))
MITFregGenes_filt = MITFregGenes[MITFregGenes %in% mitfchipTars]
MITFregGenes_filt = unique(MITFregGenes_filt)

## load TEADs si datasets
hacData = read.csv("data/Verfaillie/hac_farie_data.csv", header = TRUE, stringsAsFactors = FALSE, skip=2)
hacAp1 = hacData[,"H3K27ac_AP1_target"] =="AP1"
hacTEAD = hacData[,"H3K27ac_TEAD_target"] =="TEAD"
farAp1 = hacData[,"FAIRE_AP1_target"] =="AP1"

teadGenes = unique(hacData[ hacTEAD, "gene_ID"])
ap1Genes =  unique(hacData[hacAp1 | farAp1 , "gene_ID"])

combGenes = unique(c(teadGenes, ap1Genes))

teadSidata = read.csv("data/Verfaillie/tead_targets.csv", header=TRUE, skip = 3, stringsAsFactors = FALSE)
teadDownGenes = unique(teadSidata[teadSidata[,"log2FC.1"]<0, 1])

teadRegGenes = combGenes[combGenes %in% teadDownGenes]

allgenes = unique(c(MITFregGenes_filt, teadRegGenes))

proInvGenes = read.csv("data/Verfaillie/prolif_invasive_genes.csv", header=TRUE, stringsAsFactors = FALSE)
actualProGenes = proInvGenes[proInvGenes[,"signature"] == "proliferative", 1]
actualInvGenes = proInvGenes[proInvGenes[,"signature"] == "invasive", 1]

proInvGenes = proInvGenes[ proInvGenes[,"signature"] != "", 1]

normExpressionTCGA_filt = 
    normExpressionTCGA_ALL[proInvGenes[proInvGenes %in% rownames(normExpressionTCGA_ALL)],]

normExpressionTCGA_filt_t = t(normExpressionTCGA_filt)

corMat = cor(normExpressionTCGA_filt_t)
diag(corMat) = 0

minS = apply(corMat, 2, min)
maxS = apply(corMat, 2, max)

minFilt = minS[minS < -0.5]

## remove Y genes
maxFilt = maxS[maxS<0.96 & maxS>0.65]

wantMaxs = unique(c(names(minFilt), names(maxFilt)))

out = pheatmap(corMat[wantMaxs,wantMaxs])
clusts = cutree(out$tree_row, k=3)
proGenes = names(clusts[clusts==1])
invGenes = names(clusts[clusts==3])

proGenes = proGenes[proGenes %in% actualProGenes]
invGenes = invGenes[invGenes %in% actualInvGenes]

clustGenes = c(proGenes, invGenes)

## load TCGA data
normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}

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
hgncGenes = hgncGenes[-which(hgncGenes[,2]==""),]
hgncGenes = hgncGenes[ -which( duplicated(hgncGenes[,1])  ),]

skcmMatFinal = skcmMatFinal[rownames(skcmMatFinal) %in% hgncGenes[,1], ]
rownames(skcmMatFinal) = hgncGenes[match(rownames(skcmMatFinal), hgncGenes[,1]) , 2]

coldata = data.frame(Samp = colnames(skcmMatFinal))
rownames(coldata) = colnames(skcmMatFinal)

ddsCounts <- DESeqDataSetFromMatrix(countData = skcmMatFinal,
                                    colData = coldata,
                                    design = ~ Samp)

ntd = normTransform(ddsCounts)
normExpressionTCGA_ALL =  assay(ntd)

scaled = t(as.data.frame(apply(normExpressionTCGA_ALL, 1, normalize)))
clustGenes = clustGenes[clustGenes %in% rownames(normExpressionTCGA_ALL)]

out = pheatmap(scaled[clustGenes,])
clusts = cutree(out$tree_col, k=2)

## use Kmeans to define clusters
kmeanRes = kmeans(t(normExpressionTCGA_ALL[clustGenes,]), 2, iter.max = 20, nstart = 1 )
clusts2 = kmeanRes$cluster

allTestGenes = c(actualProGenes, actualInvGenes)
allTestGenes = allTestGenes[allTestGenes %in% rownames(normExpressionTCGA_ALL)]

# Figure 5b
qq = Rtsne(t(scaledCCLE[clustGenes, ]), labels = as.factor(clusts3), perplexity = 10, epoch=1000)

clusts3[clusts3==1] = "Invasive"
clusts3[clusts3==2] = "Proliferative"

clusterDatFrame = data.frame(tsne1 = qq$Y[,1], tsne2 = qq$Y[,2], Subtype = clusts3   )

tiff(file="samplesCLustering.tiff", width=6, height=4.2, units="in", res=300)

ggplot(clusterDatFrame, aes(x=tsne1, y=tsne2, fill=Subtype)) + xlim(-8,8) + ylim(-8,8) + geom_point(size=3, alpha=0.5, shape=21) + theme_bw() +
  theme(axis.text = element_text(size = 16), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size=14),
        legend.title = element_text(size=16))

dev.off()

## Fig 5c
## load CCLE expression data
## CCLE expression data from kallisto quantification
load("data/CCLE/DEseq2_norm_CCLE_Melanoma.Rdata")
normExpression_filt_CCLE = normExpressionCCLE[rowSums(normExpressionCCLE>0.5)>3,]

### ccle vs tcga expression
maxMITF = c("TSPAN10", "PMEL", "TUBB4A", "TRPM1", "CABLES1")
maxTead = c("AXL", "COL8A1", "THBS1", "COL5A2", "IGFBP3")

tcgaExprs = melt(normExpressionTCGA_ALL[c(maxMITF, maxTead), ])
ccleExprs = melt(normExpressionCCLE[c(maxMITF, maxTead),])

colnames(tcgaExprs) = c("Gene", "Sample", "Expression")
colnames(ccleExprs) = c("Gene", "Sample", "Expression")

ccleTcgaExprs = 
    data.frame(rbind(tcgaExprs, ccleExprs), 
               Dataset = c(rep("TCGA", nrow(tcgaExprs)), rep("CCLE", nrow(ccleExprs))) )

ccleTcgaExprs$Gene = factor(forcats::fct_rev(factor(ccleTcgaExprs$Gene)))

tiff(file="ccle_TCGA_comp.tiff", width=10, height=7, units="in", res=300)

e = ggplot(ccleTcgaExprs, aes(x = Gene, y = Expression, fill=Dataset))

e + geom_violin(position=position_dodge(0.8)) + theme_classic() + geom_boxplot(width=0.1, position=position_dodge(0.8))+
    theme(axis.text = element_text(size = 22), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size=18), legend.title = element_text(size=18),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylim(0,27) + scale_fill_brewer(palette="RdBu")

dev.off()

# Fig 5d
clusts = clusts2
y = factor(clusts)

set.seed(54328)
featList = c()

wantEpxrMat = normExpressionTCGA_ALL[allgenes[allgenes %in% rownames(normExpressionTCGA_ALL)],]

mads = apply(wantEpxrMat, 1, mad)
wantGenes = names(mads[mads>=1])

currEx = t(wantEpxrMat[wantGenes,])
bortOut = Boruta(currEx, y, maxRuns = 500)
importGenes = names(bortOut$finalDecision[bortOut$finalDecision=="Confirmed"])

importGenes = names(bortOut$finalDecision[bortOut$finalDecision=="Confirmed"])
importGenes = sub("\\.", "-", importGenes)

ccleMads = apply(normExpressionCCLE[importGenes[importGenes %in% rownames(normExpressionCCLE)],], 1, mad)
importGenes = names(ccleMads[ccleMads>1.2])

down = importGenes[importGenes %in% teadRegGenes]
up = importGenes[importGenes %in% MITFregGenes_filt]

miniMat = normExpressionTCGA_ALL[importGenes, ]
importCorMat = cor(t(miniMat))

LRexpr =  t(scaled[importGenes,names(clusts)])
LRexpr[LRexpr>0.66] = 1
LRexpr[LRexpr <= 0.66] = 0

scaled[importGenes,names(clusts)]

clusts[clusts=="Proliferative"] = 2
clusts[clusts=="Invasive"] = 1
clusts2 = names(clusts)
clusts2 = sub("-06A-.*", "", clusts2)
clusts2 = sub("-01A-.*", "", clusts2)

clusts3=clusts
names(clusts3) = clusts2

## rank by frequency observe in other class
topTEADx = sort(rowSums(t(LRexpr)[down, names(clusts[clusts==2])])/length(clusts[clusts==2]) )
topMITFx = sort(rowSums(t(LRexpr)[up, names(clusts[clusts==1])])/length(clusts[clusts==1]) )

topTEADy = sort(rowSums(t(LRexpr)[down, names(clusts[clusts==1])])/length(clusts[clusts==1]) )
topMITFy = sort(rowSums(t(LRexpr)[up, names(clusts[clusts==2])])/length(clusts[clusts==2]) )

maxExprTEAD = apply(normExpressionTCGA_ALL[names(topTEADx),], 1, function(x) quantile(x,0.9))
maxExprMITF = apply(normExpressionTCGA_ALL[names(topMITFx),], 1, function(x) quantile(x,0.9))

wantTEADNames = names(maxExprTEAD)[maxExprTEAD>10]
wantMITFNames = names(maxExprMITF)[maxExprMITF>10]

topTEADx = topTEADx[names(topTEADx) %in% wantTEADNames]
topTEADy = topTEADy[names(topTEADy) %in% wantTEADNames]

topMITFx = topMITFx[names(topMITFx) %in% wantMITFNames]
topMITFy = topMITFy[names(topMITFy) %in% wantMITFNames]

plot(topTEADx, topTEADy) +
  plot(topMITFx, topMITFy)

topMITFx = topMITFx[!(names(topMITFx) %in% c("ARNTL2", "ID2"))]

SelectedTEAD = rep("NO", length(topTEADx))
SelectedMITF = rep("NO", length(topMITFx))

names(SelectedTEAD) = names(topTEADx)
names(SelectedMITF) = names(topMITFx)

SelectedTEAD[names(SelectedTEAD) %in% maxTead] ="YES"
SelectedMITF[names(SelectedMITF) %in% maxMITF] = "YES"

df1 = data.frame(Gene = names(topTEADx), consist = topTEADy[names(topTEADx)]*100, inconsist = topTEADx*100, Selected=SelectedTEAD)
df2 = data.frame(Gene = names(topMITFx), consist = topMITFy[names(topMITFx)]*100, inconsist = topMITFx*100, Selected=SelectedMITF)

df3 = rbind(df1, df2)
df3= data.frame(df3, Class=c(rep("TEADs target genes", nrow(df1)), rep("MITF target genes", nrow(df2))) )

df3$Class = factor(df3$Class, level = c("TEADs target genes", "MITF target genes"))

tiff(file="consistentFig.tiff", width=5, height=4.7, units="in", res=300)

ggplot(df3, aes(consist, inconsist, label = Gene, fill=Selected ))+
    geom_point() + geom_text(size=2) + theme_bw() + scale_fill_manual(values=c("gray91", "yellow")) +
    theme(axis.text = element_text(size = 13), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size=11),
          legend.title = element_text(size=14), legend.position = "none", strip.text = element_text(size=13)) + 
    geom_label(size=1.5) + labs(x="% of Samples Consistent", y="% of Samples Inconsistent") + 
    facet_wrap(~Class, scale="free",nrow = 2)

dev.off()

# Fig 5e
ggplot2Dat = setNames(melt(t(LRexpr)[c(maxTead, maxMITF), c(names(clusts[clusts==1]), names(clusts[clusts==2])) ]), c('Genes', 'Samples', 'Expression'))

ggplot2Dat[,3] = as.factor(ggplot2Dat[,3])

mat = ggplot(ggplot2Dat, aes(Samples, Genes, fill= Expression)) + 
    geom_boxplot() + theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y = element_text(colour = "black"),
          axis.title.y=element_blank()) + 
    scale_fill_manual(values=c("gray88", "slateblue4"))

ggplot2Dat$Genes = forcats::fct_rev(factor(ggplot2Dat$Genes))

tiff(file="candidateGene.tiff", width=3.2, height=3, units="in", res=300)

ggplot(ggplot2Dat, aes(Samples, Genes, fill= Expression)) + 
    geom_tile() + theme(axis.title.x=element_blank(),
                           axis.text.x=element_blank(),
                           axis.ticks.x=element_blank(),
                           axis.text.y = element_text(colour = "black"),
                           axis.title.y=element_blank(), legend.position = "none") + 
    scale_fill_manual(values=c("gray88", "slateblue4"))

dev.off()

# Fig 6b
# naive bayes classifer
## CCLE expression data from counts data (more cell lines available) 
load("data/CCLE/normExpressionCCLE_final_Old_ALL.Rdata")

fibroCellLines = c("HS895T", "HS940T", "HS600T", "HS834T", "HS688AT", "HS839T", "HS934T")
wantCls = colnames(normExpressionCCLE_final_Old)[ !(colnames(normExpressionCCLE_final_Old) %in% fibroCellLines )]

scaledCCLE = t(as.data.frame(apply(normExpressionCCLE_final_Old[, wantCls], 1, normalize)))

nB = t(scaledCCLE[c(maxTead, maxMITF),])
nB2 = nB
nB[nB>0.66] = "1"
nB[nB<=0.66] = "0"

nB2[nB2>0.66] = 1
nB2[nB2<=0.66] = 0

kmeanRes = kmeans(t(scaledCCLE[clustGenesR,]), 2, iter.max = 100, nstart = 1 )
kmeanRes["A2058_SKIN"] = 1
clClusts = kmeanRes$cluster

nBexpr =  t(scaled[c(maxTead, maxMITF), names(clusts)])
nBexpr[nBexpr>0.66] = "1"
nBexpr[nBexpr <= 0.66] = "0"

trainMat = data.frame(class=clusts2, nBexpr[,c(maxTead, maxMITF)])
trainMat = data.frame(class=clusts2, nBexpr)

nb_laplace1 <- naiveBayes(class~., data=trainMat, laplace=1)
Predict <- predict(nb_laplace1, newdata = nB ) 

Predict = as.character(Predict)

Predict[Predict=="Invasive"] = "TEADs"
Predict[Predict=="Proliferative"] = "MITF"

Predict = factor(Predict)

clClusts[clClusts==1] = "MITF"
clClusts[clClusts==2] = "TEADs"

clClusts = as.factor(clClusts)

confusionMatrix(Predict, clClusts)

annotationSamples = as.character(Predict)
annotationSamples[annotationSamples=="Invasive"] = "TEADs"
annotationSamples[annotationSamples=="Proliferative"] = "MITF"

clClusts["invasive"] = "TEADs"
clClusts["proliferative"] = "MITF"

annotationSamples = data.frame(Prediction = annotationSamples, Actual = clClusts)
rownames(annotationSamples) = rownames(nB)

outTree = pheatmap(scaledCCLE[clustGenesR,c(names(clClusts[clClusts=="MITF"]), names(clClusts[clClusts=="TEADs"]))])
rorder = outTree$tree_row$labels

tiff(file="cellLinePredict.tiff", width=2.5, height=2.5, units="in", res=300)

  pheatmap(scaledCCLE[rev(rorder),c(names(clClusts[clClusts=="TEADs"]), names(clClusts[clClusts=="MITF"]))], 
           treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE, show_colnames = FALSE, 
           annotation_col = annotationSamples, cluster_rows = FALSE,
           color = magma(100),
           cluster_cols = FALSE,
           annotation_colors = list(Actual = c(MITF = "gold", TEADs = "darkorchid4"),
                                        Prediction = c(MITF = "gold", TEADs = "darkorchid4")), 
           fontsize = 8)

dev.off()

tiff(file="cellLineBinary.tiff", width=2, height=2.5, units="in", res=300)
  pheatmap(t(nB2)[,c(names(clClusts[clClusts=="TEADs"]), names(clClusts[clClusts=="MITF"]))], 
           treeheight_row = 0, treeheight_col = 0, show_colnames = FALSE, 
           annotation_col = annotationSamples, 
           color = inferno(2),
           cluster_cols = FALSE, 
           annotation_colors = list(Actual = c(MITF = "gold", TEADs = "darkorchid4"),
                                    Prediction = c(MITF = "gold", TEADs = "darkorchid4")), 
           legend = FALSE,
           fontsize = 8, annotation_legend = FALSE)

dev.off()

fibroCellLines = c("HS895T", "HS940T", "HS600T", "HS834T", "HS688AT", "HS839T", "HS934T")
wantCls = colnames(normExpressionCCLE_final_Old)[ !(colnames(normExpressionCCLE_final_Old) %in% fibroCellLines )]

scaledCCLE = t(as.data.frame(apply(normExpressionCCLE_final_Old[, wantCls], 1, normalize)))

kmeanRes = kmeans(t(scaledCCLE[clustGenesR,]), 2, iter.max = 100, nstart = 1 )
clClusts = kmeanRes$cluster

clClusts[clClusts==2] = "De-diff"
clClusts[clClusts==1] = "Diff"

# Fig 6c
brafSens = get(load("data/CCLE/cell_line_drug_info.Rdata"))
names(brafSens) = sub("_SKIN", "", names(brafSens))

drugSensDat = data.frame(Sensitivity=brafSens, Subtype=clClusts[names(brafSens)])
drugSensDat$Subtype = factor(drugSensDat$Subtype, levels = c("De-diff", "Diff"))

tiff(file="cellLines_sens.tiff", width=3.5, height=3, units="in", res=300)

ggplot(drugSensDat, aes(Subtype, Sensitivity)) + geom_boxplot() + theme_classic()+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
    scale_fill_manual(values=c("#999999", "#E69F00")) + 
    theme(axis.title.x = element_blank(), axis.text = element_text(size = 12),
          axis.title.y = element_text(size=14)) + ylab("BRAFi Sensitivity (AUC)") + 
    stat_compare_means(method = "t.test", tip.length = 0.01, size=5, label.y = 8.7,
                       label.x=1.1) + ylim(0, 9)

dev.off()

# Fig 6d
maxMITF2 = c("TSPAN10", "MITF", "PMEL", "TRIM63")
maxTeadGeneList = c("FMOD", "AXL", "COL8A1", "NRP1", "IGFBP3", "THBS1", "SLFN11", "COL5A2")
maxTeadGeneList = c("AXL", "COL8A1", "NRP1", "IGFBP3", "THBS1", "ADAMTS12", "COL5A2")
mat = read.csv("data/Tirosh/Tirosh_data.csv", header=TRUE, stringsAsFactors = FALSE)

wantGenes = c(maxMITF, maxTead)

wantMat = mat[mat[,2] %in% wantGenes, 3:ncol(mat)]

rownames(wantMat) = mat[mat[,2] %in% wantGenes, 2]

mat = as.matrix(wantMat)
mat = mat[wantGenes,]
mat = mat + 0.5
Fcs = cbind(mat[,2]-mat[,1], mat[,4]-mat[,3], mat[,6]-mat[,5], 
            mat[,8]-mat[,7], mat[,10]-mat[,9], mat[,12]-mat[,11])

qq = Fcs
qq = qq/mat[,c(1,3,5,7,9,11)]

qq[qq > 0.5] = 1
qq[qq < -0.5] = -1
qq[abs(qq)<=0.5] = 0

ptVem = melt(qq)
colnames(ptVem) = c("Gene", "Patient", "Direction")

ptVem[ptVem[,3]==1,3] = "Increase"
ptVem[ptVem[,3]==0,3] = "No change"
ptVem[ptVem[,3]==-1,3] = "Decrease"

Class = rep("MITF", nrow(ptVem))
Class[ptVem[,1] %in% maxTead] = "TEADs"

ptVem = data.frame(ptVem, Class=Class)

ptVem$Direction = factor(ptVem$Direction, levels=c("Increase", "No change", "Decrease"))
ptVem$Class = factor(ptVem$Class, levels=c("TEADs", "MITF"))
ptVem$Gene = forcats::fct_rev(factor(ptVem$Gene))

##patientVem
tiff(file="patient_vem.tiff", width=5.5, height=3, units="in", res=300)

ggplot(ptVem, aes(x=Patient, y=Gene)) +
    geom_point(aes(colour=Direction), size=4) + 
    scale_color_manual(values=c("gold", "lightgray", "black")) +
    # theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    scale_x_discrete(limits=1:6, 
                     labels = c("1","2","3","4","5","6"))
dev.off()

# Fig 6e
clClusts = clClusts[( names(clClusts) %in% wantCls)]

immuneScores = read.csv("data/Jerby-Arnon/ccle_immune_sig.csv", skip=2, stringsAsFactors = FALSE)
rownames(immuneScores) = immuneScores[,1]
wantCls = immuneScores[immuneScores[,1]%in% names(clClusts), 1]

drugImmune = data.frame(Scores=immuneScores[wantCls, 2], Subtype=clClusts[wantCls])
drugImmune$Subtype = factor(drugImmune$Subtype, levels = c("De-diff", "Diff"))

tiff(file="cellLines_immune_sig.tiff", width=3.1, height=3, units="in", res=300)

ggplot(drugImmune, aes(Subtype, Scores)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c("#999999", "#E69F00")) + theme_classic() +
  theme(axis.title.x = element_blank(), axis.text = element_text(size=12), 
        axis.title.y=element_text(size=14)) + ylab("Cancer cell ICI resistance score") + 
  stat_compare_means(method = "t.test", size=5, label.y = 1.15,
                     label.x = 1.15) + ylim(-1, 1.2)

dev.off()

# Fig 6f
rogerExp = read.csv("data/Lo/Roger_lo_pd1.csv", header=TRUE, stringsAsFactors = FALSE)
rownames(rogerExp) = rogerExp[,1]
rogerExp = rogerExp[,-1]
rogerExp = rogerExp[,grep("baseline", colnames(rogerExp))]
colnames(rogerExp) = sub("\\..*", "", colnames(rogerExp))

sampInfos = read.csv("data/Lo/SraRunTable.txt", header=TRUE, stringsAsFactors = FALSE)

responders = sampInfos[c(grep("Complete", sampInfos[,"anti.pd.1_response"]), grep("Partial", sampInfos[,"anti.pd.1_response"]) ), "patient_id"]
nonResponders = sampInfos[grep("Progressive", sampInfos[,"anti.pd.1_response"]), "patient_id"]

responders = c(responders, "Pt27A")

totReads = colSums(rogerExp)
totReadsMat = t(replicate(nrow(rogerExp), totReads))

wantPts = c(responders, nonResponders)
wantPts = unique(wantPts[wantPts %in% colnames(rogerExp)])

rogerExp = as.matrix(rogerExp)

rogerExp[c("AXL", "CDH1", "DCT", "PMEL", "TRIM63"), wantPts]

rogerExp[names(topMITFx), c(wantPts, "Pt27A")]
sort(rowMeans(rogerExp[names(topMITFx), wantPts]))

rogerExp[names(topTEADx), wantPts]
sort(rowMeans(rogerExp[names(topTEADx), wantPts]))

testExpr = rogerExp[c(maxMITF,maxTead), wantPts]

scaledRoger = t(as.data.frame(apply(rogerExp, 1, normalize)))+0.01

maxGenes = c(maxTead, maxMITF)
qq = scaledRoger[maxGenes[maxGenes %in% rownames(scaledRoger)], wantPts]

tiff(file="pd1heatmap.tiff", width=2.5, height=3, units="in", res=300)
    pheatmap(qq, 
         cluster_cols = FALSE, cluster_rows = FALSE, color=inferno(100), show_colnames = TRUE)

dev.off()

qq2 = qq
qq2[qq2>0.25] = 1
qq2[qq2<0.25] = 0

tiff(file="pd1heatmap.tiff", width=2.8, height=3, units="in", res=300)
    pheatmap(qq2[,rev(colnames(qq2))], cluster_rows = FALSE, cluster_cols = FALSE, color = inferno(2), legend = FALSE, fontsize = 8, show_colnames = TRUE)
dev.off()

mitfCounts = colSums(qq2[6:10,])
teadCounts = colSums(qq2[1:5,])

tead = data.frame(Counts = teadCounts[wantPts], Status = c(rep("responders",14), rep("non-responders",12))   )
mitf = data.frame(Counts = mitfCounts[wantPts], Status = c(rep("responders",14), rep("non-responders",12))   )

tead$Status = factor(tead$Status, levels = c("non-responders", "responders"))
mitf$Status = factor(mitf$Status, levels = c("non-responders", "responders"))

library(ggpubr)

p1 = ggplot(tead, aes(x=Status, y=Counts, fill=Status)) + 
    geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
    scale_fill_manual(values=c("#999999", "#E69F00")) + theme_classic() +
    theme(axis.text = element_text(size = 14), axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_text(size=16), 
          legend.text = element_text(size=14),
          legend.title = element_text(size=16)) + ylim(0,7) + 
    stat_compare_means(method = "t.test", tip.length = 0.01, size=5,
                       symnum.args=list(cutpoints=c(0,0.05,1), symbols=c("*","ns")), label.y = 6.1,
                       label.x = 1.2)
    
p2 = ggplot(mitf, aes(x=Status, y=Counts, fill=Status)) + geom_boxplot() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
    scale_fill_manual(values=c("#999999", "#E69F00")) + theme_classic() +
    theme(axis.text = element_text(size = 14), axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16)) + ylim(0,5) + 
    stat_compare_means(method = "t.test", tip.length = 0.01, size=5,
                       symnum.args=list(cutpoints=c(0,0.05,1), symbols=c("*","ns")), label.y = 4.35,
                       label.x = 1.2)

tiff(file="pd1comp.tiff", width=6, height=4, units="in", res=300)
  ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE)
dev.off()
