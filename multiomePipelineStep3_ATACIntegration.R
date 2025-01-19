rm(list = ls())

library(ggsci)
library(future)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(SeuratData)


plan("multiprocess", workers = 10)
combined = readRDS("/merged_unprocessedObj.rds")

combined_filtered = subset(combined,subset = nCount_ATAC < 15000 & nCount_RNA < 6000 & percent.mt < 10 & percent.ribo < 10 & TSS.enrichment > 1 & nucleosome_signal < 2)

DefaultAssay(combined_filtered) = "ATAC"
combined_filtered = FindTopFeatures(combined_filtered, min.cutoff = 50)
combined_filtered =  RunTFIDF(combined_filtered)
combined_filtered <- RunSVD(combined_filtered)

seur_list = SplitObject(combined_filtered,split.by = "poolID")
print("Finding integration anchors")
integration.anchors <- FindIntegrationAnchors(object.list = seur_list,anchor.features = rownames(seur_list$pool52),reduction = "rlsi",dims = 2:30)
print("Integrate Embeddings")
integrated_atac = IntegrateEmbeddings(anchorset = integration.anchors,reductions = combined_filtered[["lsi"]], new.reduction.name = "integrated_lsi",dims.to.integrate=1:30)
print("RunUMAP, find neighbours and clusters")
integrated_atac = RunUMAP(integrated_atac, reduction = "integrated_lsi", dims = 2:30)
integrated_atac = FindNeighbors(integrated_atac,reduction = "integrated_lsi",dims = 2:30)
integrated_atac = FindClusters(integrated_atac,verbose = FALSE, algorithm = 3,resolution = 0.2)
#Idents(integrated_atac) = "ATAC_snn_res.0.2"

print("finished integration, now gene activity assay creation")
gene.activities <- GeneActivity(integrated_atac)
integrated_atac[["genes"]]  = CreateAssayObject(counts = gene.activities)
integrated_atac <- NormalizeData(
  object = integrated_atac,
  assay = 'genes',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated_atac$nCount_RNA)
)
saveRDS(integrated_atac,file="ATACIntegrated.rds")

