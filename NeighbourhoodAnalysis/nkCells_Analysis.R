library(Seurat)
library(Signac)
library(miloR)
library(patchwork)
library(scater)
library(ggplot2)
library(edgeR)
library(scran)
library(dplyr)
library(DESeq2)

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$poolID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)


obj <- IntegrateLayers(object = obj, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:10)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, dims = 1:10, reduction = "integrated.rpca")
objMarkers = FindAllMarkers(obj,only.pos = T)
objMarkers_signi = subset(objMarkers, objMarkers$p_val_adj < 0.05)

DefaultAssay(obj) = "chromvar"

nk_nonICU_T2_Recov_motif = FindMarkers(obj,ident.1 = "nonICU_T2", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
nk_nonICU_T3_Recov_motif = FindMarkers(obj,ident.1 = "nonICU_T3", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
nk_nonICU_T4_Recov_motif = FindMarkers(obj,ident.1 = "nonICU_T4", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
nk_nonICU_T5_Recov_motif = FindMarkers(obj,ident.1 = "nonICU_T5", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")

nk_ICU_T2_Recov = FindMarkers(obj,ident.1 = "ICU_T2", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
nk_ICU_T3_Recov = FindMarkers(obj,ident.1 = "ICU_T3", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
nk_ICU_T4_Recov = FindMarkers(obj,ident.1 = "ICU_T4", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
nk_ICU_T5_Recov = FindMarkers(obj,ident.1 = "ICU_T5", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")


#####Milo analysis######
sce = as.SingleCellExperiment(obj)
reducedDim(sce,"UMAP") = Embeddings(obj,reduction = "umap")
reducedDim(sce,"PCA") = Embeddings(obj,reduction = "pca")
sce_milo = Milo(sce)

sce_milo = buildGraph(sce_milo,k=50,d=50)
sce_milo = makeNhoods(sce_milo,prop = 0.1, k = 50, d=50, refined = TRUE, refinement_scheme = "graph")
plotNhoodSizeHist(sce_milo)
sce_milo = countCells(sce_milo, meta.data = data.frame(colData(sce_milo)), sample = "donorID_TimePoint")
head(nhoodCounts(sce_milo))

designMat = data.frame(colData(sce_milo))[,c("donorID_TimePoint","TP_RecovLC","Recov_LC"), drop=FALSE]
#designMat$poolID = as.factor(designMat$poolID)
designMat = distinct(designMat)
rownames(designMat) = designMat$donorID_TimePoint
designMat = designMat[colnames(nhoodCounts(sce_milo)), , drop=FALSE]
contrast.0 = c("TP_RecovLCnonICU_T1 - TP_RecovLCHC_HC")
contrast.1 = c("TP_RecovLCnonICU_T2 - TP_RecovLCnonICU_T1")
contrast.2 = c("TP_RecovLCnonICU_T3 - TP_RecovLCnonICU_T1")
contrast.3 = c("TP_RecovLCnonICU_T4 - TP_RecovLCnonICU_T1")
contrast.4 = c("TP_RecovLCnonICU_T5 - TP_RecovLCnonICU_T1")
contrast.5 =  c("TP_RecovLCnonICU_T3 - TP_RecovLCR_Recov")
contrast.6 =  c("TP_RecovLCnonICU_T4 - TP_RecovLCR_Recov")
contrast.7 = c("TP_RecovLCICU_T3 - TP_RecovLCR_Recov")
contrast.8 = c("TP_RecovLCICU_T4 - TP_RecovLCR_Recov")
sce_milo = buildNhoodGraph(sce_milo)

perTimePointCompare = function(contrast){
  da_result <- testNhoods(sce_milo, design = ~ 0 + TP_RecovLC, design.df = designMat, model.contrasts = contrast, fdr.weighting="graph-overlap")
  da_result = groupNhoods(sce_milo,da_result,max.lfc.delta = 1,overlap = 10)
  da_result = annotateNhoods(sce_milo, da_result, coldata_col = "cellCluster")
  da_result = annotateNhoods(sce_milo, da_result, coldata_col = "RNA_snn_res.0.4")
  return(da_result)
}

nonICU_T1vsHC = perTimePointCompare(contrast.0)
nonICU_T2vsT1 = perTimePointCompare(contrast.1)
nonICU_T3vsT1 = perTimePointCompare(contrast.2)
nonICU_T4vsT1 = perTimePointCompare(contrast.3)
nonICU_T5vsT1 = perTimePointCompare(contrast.4)
nonICU_T3vsRecov = perTimePointCompare(contrast.5)
nonICU_T4vsRecov = perTimePointCompare(contrast.6)
ICU_T3vsRecov = perTimePointCompare(contrast.7)
ICU_T4vsRecov = perTimePointCompare(contrast.8)

