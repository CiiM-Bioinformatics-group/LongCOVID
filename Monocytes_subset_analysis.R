library(Seurat)
library(Signac)
library(miloR)
library(patchwork)
library(scater)
library(ggplot2)
library(scran)
library(dplyr)
library(destiny)
library(ggthemes)
library(ggbeeswarm)
library(tidyverse)
library(tradeSeq)
library(gridExtra)
library(slingshot)
library(mclust)
library(BiocParallel)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 8 


cd14Mono[["RNA"]] <- split(cd14Mono[["RNA"]], f = cd14Mono$poolID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
cd14Mono = subset(cd14Mono, subset = poolID != "pool4")
cd14Mono[["integrated.rpca"]] = NULL
cd14Mono[["umap.rpca"]] = NULL
# run standard anlaysis workflow
cd14Mono <- NormalizeData(cd14Mono)
cd14Mono <- FindVariableFeatures(cd14Mono)
cd14Mono <- ScaleData(cd14Mono)
cd14Mono <- RunPCA(cd14Mono)
cd14Mono <- IntegrateLayers(object = cd14Mono, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
cd14Mono[["RNA"]] <- JoinLayers(cd14Mono[["RNA"]])

cd14Mono <- FindNeighbors(cd14Mono, reduction = "integrated.rpca", dims = 1:10)
cd14Mono <- FindClusters(cd14Mono, resolution = 0.2)
cd14Mono <- RunUMAP(cd14Mono, dims = 1:10, reduction = "integrated.rpca")

#####cd14Mono = subset(cd14Mono,subset = RNA_snn_res.0.5 != 6) Since this was CD8+ Monocytes, to removed them###
cd14MonoMarkers = FindAllMarkers(cd14Mono,only.pos = T)
cd14MonoMarkers_signi = subset(cd14MonoMarkers, cd14MonoMarkers$p_val_adj < 0.05)

#### Milo analysis######
cd14Mono$cellCluster = ifelse(cd14Mono$RNA_snn_res.0.2 == 0, "MC1",
                              ifelse(cd14Mono$RNA_snn_res.0.2 == 1, "MC2", 
                                     ifelse(cd14Mono$RNA_snn_res.0.2 == 2, "MC3","MC4")))

sce_myeloid = as.SingleCellExperiment(cd14Mono)
reducedDim(sce_myeloid,"UMAP") = Embeddings(cd14Mono,reduction = "umap")
reducedDim(sce_myeloid,"PCA") = Embeddings(cd14Mono,reduction = "pca")
sce_myeloid_milo = Milo(sce_myeloid)

sce_myeloid_milo = buildGraph(sce_myeloid_milo,k=50,d=50)
sce_myeloid_milo = makeNhoods(sce_myeloid_milo,prop = 0.1, k = 50, d=50, refined = TRUE, refinement_scheme = "graph")
plotNhoodSizeHist(sce_myeloid_milo)
sce_myeloid_milo = countCells(sce_myeloid_milo, meta.data = data.frame(colData(sce_myeloid_milo)), sample = "donorID_TimePoint")
head(nhoodCounts(sce_myeloid_milo))

#### DA using negative binomial GLM

designMat = data.frame(colData(sce_myeloid_milo))[,c("donorID_TimePoint","AcuteCovidStatus_TimePoint2","Recov_LC"), drop=FALSE]
designMat = distinct(designMat)
rownames(designMat) = designMat$donorID_TimePoint
designMat = designMat[colnames(nhoodCounts(sce_myeloid_milo)), , drop=FALSE]

##### contrasts ######

contrasts.1 = c("AcuteCovidStatus_TimePoint2nonICU_T2 - AcuteCovidStatus_TimePoint2nonICU_T1")
contrasts.2 = c("AcuteCovidStatus_TimePoint2nonICU_T3 - AcuteCovidStatus_TimePoint2nonICU_T1")
contrasts.3 = c("AcuteCovidStatus_TimePoint2nonICU_T4 - AcuteCovidStatus_TimePoint2nonICU_T1")
contrasts.4 = c("AcuteCovidStatus_TimePoint2nonICU_T5 - AcuteCovidStatus_TimePoint2nonICU_T1")
contrasts.5 = c("Recov_LCnonICU_LC - Recov_LCR_Recov")


contrasts.6 = c("AcuteCovidStatus_TimePoint2ICU_T2 - AcuteCovidStatus_TimePoint2ICU_T1")
contrasts.7 = c("AcuteCovidStatus_TimePoint2ICU_T3 - AcuteCovidStatus_TimePoint2ICU_T1")
contrasts.8 = c("AcuteCovidStatus_TimePoint2ICU_T4 - AcuteCovidStatus_TimePoint2ICU_T1")
contrasts.9 = c("AcuteCovidStatus_TimePoint2ICU_T5 - AcuteCovidStatus_TimePoint2ICU_T1")
contrasts.10 = c("Recov_LCICU_LC - Recov_LCR_Recov")

perTimePointCompare = function(cont){
 da_results =  testNhoods(sce_myeloid_milo, design = ~0 + AcuteCovidStatus_TimePoint2, 
             design.df = designMat, 
             model.contrasts = cont,
             fdr.weighting = "graph-overlap")
 da_results = annotateNhoods(sce_myeloid_milo, da_results, coldata_col = "cellCluster")
 return(da_results)
}

compareRecov = function(cont){
  da_results =  testNhoods(sce_myeloid_milo, design = ~0 + Recov_LC, 
                           design.df = designMat, 
                           model.contrasts = cont,
                           fdr.weighting = "graph-overlap")
  da_results = annotateNhoods(sce_myeloid_milo, da_results, coldata_col = "cellCluster")
  return(da_results)
}

sce_myeloid_milo = buildNhoodGraph(sce_myeloid_milo)
nonICU_T2vsT1 = perTimePointCompare(contrasts.1)
nonICU_T3vsT1 = perTimePointCompare(contrasts.2)
nonICU_T4vsT1 = perTimePointCompare(contrasts.3)
nonICU_T5vsT1 = perTimePointCompare(contrasts.4)
nonICU_LCvsR = compareRecov(contrasts.5)

ICU_T2vsT1 = perTimePointCompare(contrasts.6)
ICU_T3vsT1 = perTimePointCompare(contrasts.7)
ICU_T4vsT1 = perTimePointCompare(contrasts.8)
ICU_T5vsT1 = perTimePointCompare(contrasts.9)
ICU_LCvsR = compareRecov(contrasts.10)


##### AUC score calculation
Monocounts = GetAssayData(cd14Mono,slot = "counts") ### can be any celltype R object.
Mono_cell_rankings <- AUCell_buildRankings(Monocounts)

#### optional for AUC of genesets #####
geneSets = getGmt("h.all.v2022.1.Hs.symbols.gmt")
geneSets <- subsetGeneSets(geneSets, rownames(Monocounts)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))
hallmark_cells_AUC = AUCell_calcAUC(geneSets, Mono_cell_rankings)

