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
library(AUCell)

cd8T = subset(bigObj,subset = CellType == "CD8 T cells")
cd8T[["RNA"]] <- split(cd8T[["RNA"]], f = cd8T$poolID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
cd8T <- NormalizeData(cd8T)
cd8T <- FindVariableFeatures(cd8T)
cd8T <- ScaleData(cd8T)
cd8T <- RunPCA(cd8T)


cd8T <- IntegrateLayers(object = cd8T, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
cd8T[["RNA"]] <- JoinLayers(cd8T[["RNA"]])

cd8T <- FindNeighbors(cd8T, reduction = "integrated.rpca", dims = 1:10)###check if anything changes
cd8T <- FindClusters(cd8T, resolution = 0.4)
cd8T <- RunUMAP(cd8T, dims = 1:10, reduction = "integrated.rpca")
cd8TMarkers = FindAllMarkers(cd8T,only.pos = T)
cd8TMarkers_signi = subset(cd8TMarkers, cd8TMarkers$p_val_adj < 0.05)

DefaultAssay(cd8T) = "chromvar"

cd8T_nonICU_T2_Recov_motif = FindMarkers(cd8T,ident.1 = "nonICU_T2", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
cd8T_nonICU_T3_Recov_motif = FindMarkers(cd8T,ident.1 = "nonICU_T3", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
cd8T_nonICU_T4_Recov_motif = FindMarkers(cd8T,ident.1 = "nonICU_T4", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
cd8T_nonICU_T5_Recov_motif = FindMarkers(cd8T,ident.1 = "nonICU_T5", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")

cd8T_ICU_T2_Recov = FindMarkers(cd8T,ident.1 = "ICU_T2", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
cd8T_ICU_T3_Recov = FindMarkers(cd8T,ident.1 = "ICU_T3", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
cd8T_ICU_T4_Recov = FindMarkers(cd8T,ident.1 = "ICU_T4", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")
cd8T_ICU_T5_Recov = FindMarkers(cd8T,ident.1 = "ICU_T5", ident.2 = "R_Recov",mean.fxn = rowMeans,fc.name = "avg_diff")


DefaultAssay(cd8T) = "RNA"
sce_cd8T = as.SingleCellExperiment(cd8T)
reducedDim(sce_cd8T,"UMAP") = Embeddings(cd8T,reduction = "umap")
reducedDim(sce_cd8T,"PCA") = Embeddings(cd8T,reduction = "integrated.rpca")
sce_cd8T_milo = Milo(sce_cd8T)

sce_cd8T_milo = buildGraph(sce_cd8T_milo,k=50,d=50)
sce_cd8T_milo = makeNhoods(sce_cd8T_milo,prop = 0.1, k = 50, d=50, refined = TRUE, refinement_scheme = "graph")
plotNhoodSizeHist(sce_cd8T_milo)
sce_cd8T_milo = countCells(sce_cd8T_milo, meta.data = data.frame(colData(sce_cd8T_milo)), sample = "donorID_TimePoint")
head(nhoodCounts(sce_cd8T_milo))

designMat = data.frame(colData(sce_cd8T_milo))[,c("donorID_TimePoint","AcuteCovStatus2"), drop=FALSE]
#designMat$poolID = as.factor(designMat$poolID)
designMat = distinct(designMat)
rownames(designMat) = designMat$donorID_TimePoint
designMat = designMat[colnames(nhoodCounts(sce_cd8T_milo)), , drop=FALSE]
contrast.0.0 = c("AcuteCovStatus2T1 - AcuteCovStatus2HC_HC")
contrast.0.1 = c("AcuteCovStatus2nonICU_T2 - AcuteCovStatus2nonICU_T1")
contrast.0.2 = c("AcuteCovStatus2nonICU_T3 - AcuteCovStatus2nonICU_T1")
contrast.0.3 = c("AcuteCovStatus2nonICU_T4 - AcuteCovStatus2nonICU_T1")
contrast.0.4 = c("AcuteCovStatus2nonICU_T5 - AcuteCovStatus2nonICU_T1") 
contrast.0.5 = c("AcuteCovStatus2nonICU_T2 - AcuteCovStatus2R_Recov")
contrast.0.6 = c("AcuteCovStatus2nonICU_T3 - AcuteCovStatus2R_Recov")
contrast.0.7 = c("AcuteCovStatus2nonICU_T4 - AcuteCovStatus2R_Recov")
contrast.0.8 = c("AcuteCovStatus2nonICU_T5 - AcuteCovStatus2R_Recov")
contrast.1.5 = c("AcuteCovStatus2ICU_T2 - AcuteCovStatus2R_Recov")
contrast.1.6 = c("AcuteCovStatus2ICU_T3 - AcuteCovStatus2R_Recov")
contrast.1.7 = c("AcuteCovStatus2ICU_T4 - AcuteCovStatus2R_Recov")
contrast.1.8 = c("AcuteCovStatus2ICU_T5 - AcuteCovStatus2R_Recov")


sce_cd8T_milo = buildNhoodGraph(sce_cd8T_milo)

perTimePointCompare = function(contrast){
  da_result <- testNhoods(sce_cd8T_milo, design = ~ 0 + AcuteCovStatus2, design.df = designMat, model.contrasts = contrast, fdr.weighting="graph-overlap")
  da_result = groupNhoods(sce_cd8T_milo,da_result,max.lfc.delta = 1)
  da_result = annotateNhoods(sce_cd8T_milo, da_result, coldata_col = "RNA_snn_res.0.4")
  return(da_result)
}

cd8T_nonICU_T1vsHC = perTimePointCompare(contrast.0.0)
cd8T_nonICU_T2vsT1 = perTimePointCompare(contrast.0.1)
cd8T_nonICU_T3vsT1 = perTimePointCompare(contrast.0.2)
cd8T_nonICU_T4vsT1 = perTimePointCompare(contrast.0.3)
cd8T_nonICU_T5vsT1 = perTimePointCompare(contrast.0.4)
cd8T_nonICU_T2vsR = perTimePointCompare(contrast.0.5)
cd8T_nonICU_T3vsR = perTimePointCompare(contrast.0.6)
cd8T_nonICU_T4vsR = perTimePointCompare(contrast.0.7)
cd8T_nonICU_T5vsR = perTimePointCompare(contrast.0.8)

cd8T_ICU_T2vsR = perTimePointCompare(contrast.1.5)
cd8T_ICU_T3vsR = perTimePointCompare(contrast.1.6)
cd8T_ICU_T4vsR = perTimePointCompare(contrast.1.7)
cd8T_ICU_T5vsR = perTimePointCompare(contrast.1.8)


