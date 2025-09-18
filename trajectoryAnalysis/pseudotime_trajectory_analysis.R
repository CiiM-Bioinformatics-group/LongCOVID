library(Seurat)
library(Signac)
#library(miloR)
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
BPPARAM$workers <- 6 



######## all dge heatmaps and results are saved in "cd14Mono_pseudotimeAnalysis_15k_tradeseq.RData"

DefaultAssay(cd14Mono) = "RNA"
sce_LC = as.SingleCellExperiment(cd14Mono)
reducedDim(sce_LC,"UMAP") = Embeddings(cd14Mono,reduction = "umap")
reducedDim(sce_LC,"PCA") = Embeddings(cd14Mono,reduction = "integrated.rpca")
#reducedDim(sce_LC,"rpca") = Embeddings(cd14Mono,reduction = "integrated.rpca")
pca = reducedDim(sce_LC,"PCA")


# Make a diffusion map.
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
signs = find_sigmas(as.matrix(logcounts(sce_LC)))
dm <- DiffusionMap(sce_LC,n_pcs = 30,signs)

########adding DCs to singlecell object
DCs = as.matrix(eigenvectors(dm))
reducedDim(sce_LC,"DC") = DCs

#### clustering and trajectory analysis 
rd1 <- DCs[,1:15]
cl1 = Mclust(rd1)$classification
colData(sce_LC)$GMM <-cl1
cd14Mono$GMM = cl1

df <- data.frame(DC1 = DCs[,1], DC2 = DCs[,2],DC3 = DCs[,3],DC4 = DCs[,4],DC5 =  DCs[,5],
                 RecovStatus = colData(sce_LC)$Recov_LC,
                 CellType = colData(sce_LC)$cellCluster,
                 CovidStatus = colData(sce_LC)$AcuteCovidStatus,
                 AcuteCovTP = colData(sce_LC)$AcuteCovidStatus_TimePoint2,
                 mclustClusters = as.factor(colData(sce_LC)$GMM))

cd14Mono[["DC"]] = CreateDimReducObject(embeddings = as.matrix(df[,1:5]),key = "DC_",assay = DefaultAssay(cd14Mono))

sce_LC <-slingshot(sce_LC, clusterLabels = 'GMM',reducedDim = "DC",start.clus = 4,approx_points = 150)

crv1 <- getCurves(SlingshotDataSet(sce_LC))

curves = slingCurves(crv1,as.df = TRUE)

