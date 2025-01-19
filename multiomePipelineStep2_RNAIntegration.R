rm(list = ls())

library(ggsci)
library(future)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(SeuratData)

combined = readRDS( "merged_unprocessedObj.rds")

##### post evaluation of violinPlot, across all pools and as how it is described in the vignette, these thresholds were chosen #####
combined_filtered = subset(combined,subset = nCount_ATAC < 15000 & nCount_RNA < 6000 & percent.mt < 10 & percent.ribo < 10 & TSS.enrichment > 1 & nucleosome_signal < 2)

###### IMPORTANT!!!! filter based on the violinPlot before proceeding ######

seur_list = SplitObject(combined_filtered,split.by = "poolID")

print("RNA normalization starting")

seur_list <- lapply(X = seur_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

print("RNA normalization finished, starting integration")

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seur_list)

 #for (i in seur_list) {
 #  print(i)

 #  i = ScaleData(i,features = features,verbose = FALSE)
 #  i = RunPCA(i,features = features, verbose = FALSE)
#}
seur_list <- lapply(X = seur_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = seur_list, anchor.features = features, reduction = "rpca",dims = 1:50)


integrated <- IntegrateData(anchorset = immune.anchors, dims = 1:50)
integrated = ScaleData(integrated)
integrated = RunPCA(integrated,verbose = FALSE)
integrated = RunUMAP(integrated,reduction = "pca",dims = 1:50,return.model = TRUE)
integrated = FindNeighbors(integrated,reduction = "pca",dims = 1:50)
integrated = FindClusters(integrated,resolution=0.2)
print("Finished integration,saving object")
saveRDS(integrated,file="RNAIntegrated.rds")


DefaultAssay(integrated) = "RNA"
globalMarkers = FindAllMarkers(integrated,only.pos = T)
#write.csv(globalMarkers,file="/vol/projects/skumar/LongCovid/freeze2Data/multiomeObjects/freeze2_dgeRNA.csv",quote = F,row.names = F)

