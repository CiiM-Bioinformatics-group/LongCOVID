library(Seurat)
library(ggplot2)
library(future)

sampleIDs =  c("pool-LR1","pool-LR2")
finalMetadata = read.csv("finalMetadata.csv",header = T)
for (i in sampleIDs) {
  poolSamples = subset(finalMetadata, finalMetadata$CorrectPool == i)
  temp = read.table(paste0("soupOrCell/",i,"/clusters.tsv"),header = T)
  temp = subset(temp,temp$status == "singlet")
  tempMerged = merge(poolSamples,temp,by.y="assignment",by.x="RNA_Cluster")
  tempMerged2 = tempMerged[,c(1:7)]
  write.csv(tempMerged2,file=paste0("soupOrCell/",i,"SingletsFile.csv"),quote = F,row.names = F)
}

seur_list = list()
for (i in sampleIDs) {
  
  print(paste0('Processing pool: ', i))
  
  # Reading in singlets pre-processed file#
  singlets <- read.csv(paste0("soupOrCell/",i,"SingletsFile.csv"),header = T)
  rownames(singlets) = singlets$barcode
  # Reading in 10X data which is a list of two data matrices: GEX and ATAC
  data <- Read10X_h5(paste0("countsFolder/", i, "/outs/filtered_feature_bc_matrix.h5"))
  seur <- CreateSeuratObject(data,
                             assay = 'RNA',
                             min.cells = 5,
                             meta.data = singlets)
  seur[["percent.mt"]] =  PercentageFeatureSet(seur, pattern = "^MT-")
  seur[["percent.ribo"]] = PercentageFeatureSet(seur,pattern = "^RP[SL][0-9]+$")
  seur_list[[i]] <- seur
  
}


mergedSeurObj = merge(x=seur_list$`pool-LR1`,y=seur_list$`pool-LR2`,add.cell.ids = c("LR1","LR2"))
mergedSeurObj = subset(mergedSeurObj,CorrectPool %in% c("pool-LR1","pool-LR2")) ### removing cells which were doublets or mess

VlnPlot(mergedSeurObj,features = c("nCount_RNA","nFeature_RNA","percent.ribo","percent.mt"),ncol = 3,pt.size=0,group.by = "CorrectPool")

mergedSeurObj_filtered = subset(mergedSeurObj, nCount_RNA < 8000 & nFeature_RNA < 3500 & percent.mt < 20) ## based on the violin plot above

mergedSeurObj_filtered <- NormalizeData(mergedSeurObj_filtered)
mergedSeurObj_filtered <- FindVariableFeatures(mergedSeurObj_filtered)
mergedSeurObj_filtered <- ScaleData(mergedSeurObj_filtered)
mergedSeurObj_filtered <- RunPCA(mergedSeurObj_filtered)
plan("multicore", workers = 8)
options(future.globals.maxSize = 30 * 1024 ^ 3)
mergedSeurObj_filtered <- IntegrateLayers(
  object = mergedSeurObj_filtered, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

mergedSeurObj_filtered <- FindNeighbors(mergedSeurObj_filtered, reduction = "integrated.rpca", dims = 1:30)
mergedSeurObj_filtered <- FindClusters(mergedSeurObj_filtered, resolution = 0.2, cluster.name = "rpca_clusters")
mergedSeurObj_filtered <- FindClusters(mergedSeurObj_filtered, resolution = 0.5)
mergedSeurObj_filtered <- RunUMAP(mergedSeurObj_filtered, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

mergedSeurObj_filtered <- JoinLayers(mergedSeurObj_filtered)
mergedSeurObj_filtered

mergedSeurObj_filtered$MajorCellType = ifelse(mergedSeurObj_filtered$rpca_clusters == 0, "CD14 Mono", 
                                              ifelse(mergedSeurObj_filtered$rpca_clusters %in% c(1,2,7), "CD4 T",
                                                     ifelse(mergedSeurObj_filtered$rpca_clusters == 3, "NK",
                                                            ifelse(mergedSeurObj_filtered$rpca_clusters == 4, "CD8 T",
                                                                   ifelse(mergedSeurObj_filtered$rpca_clusters == 5, "B",
                                                                          ifelse(mergedSeurObj_filtered$rpca_clusters == 6, "CD16 Mono",
                                                                                 "other"))))))
DimPlot(mergedSeurObj_filtered,group.by = "MajorCellType", label = T)

saveRDS(mergedSeurObj_filtered, file = "integrated_annotated.rds")

### subsetting and re-integrating CD14 Mono cells only####
cd14Mono = subset(mergedSeurObj_filtered, MajorCellType == "CD14 Mono")
cd14Mono[["RNA"]] <- split(cd14Mono[["RNA"]], f = cd14Mono$CorrectPool)
cd14Mono <- NormalizeData(cd14Mono)
cd14Mono <- FindVariableFeatures(cd14Mono)
cd14Mono <- ScaleData(cd14Mono)
cd14Mono <- RunPCA(cd14Mono)
cd14Mono <- IntegrateLayers(object = cd14Mono, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                            verbose = FALSE)
# re-join layers after integration
cd14Mono[["RNA"]] <- JoinLayers(cd14Mono[["RNA"]])
cd14Mono <- FindNeighbors(cd14Mono, reduction = "integrated.cca", dims = 1:15)
cd14Mono <- RunUMAP(cd14Mono,dims = 1:15)
DimPlot(cd14Mono,reduction = "umap")
cd14Mono = FindClusters(cd14Mono,resolution = 0.2)

cd14Mono$Clusters = ifelse(cd14Mono$RNA_snn_res.0.2 %in% c(0,4),"Clust0",ifelse(cd14Mono$RNA_snn_res.0.2 == 1, "Clust1","Clust2"))
Idents(cd14Mono) = "Clusters"


###### get AUC scores per cell #####
library(AUCell)
library(GSEABase)
myeloidcounts = GetAssayData(cd14Mono,layer = "counts")
myeloid_cell_rankings <- AUCell_buildRankings(myeloidcounts)
myeloidCells_MC4AUC = AUCell_calcAUC(allmarkers_MC4_subset$gene,myeloid_cell_rankings)
myeloidCells_MC4AUCscore = t(getAUC(myeloidCells_MC4AUC))
cd14Mono$MC4AUC =myeloidCells_cs4AUCscore[,1]


####plotting#####
DimPlot(cd14Mono,group.by = "Clusters")+NoAxes() + scale_color_manual(values = c("#FFbcb4","#C77CFF","#ffd64c"))

cd14MonoMeta_clust1 = subset(cd14MonoMeta,cd14MonoMeta$Clusters == "Clust1")

ggplot(cd14MonoMeta,aes(x = Clusters,y = MC4AUC,fill = Clusters))+
  geom_violin(scale = "width")+theme_classic()+geom_hline(yintercept = 0.08,linetype = "dashed")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12,color = "black",angle = 45, hjust = 1,vjust = 1))+
  xlab("")+ylab("MC4 AUC")+scale_fill_manual(values = c("#FFbcb4","#C77CFF","#ffd64c"))+
  stat_compare_means(comparisons = list(c("Clust1","Clust0"),c("Clust1","Clust2")),aes(label = after_stat(p.signif)))


ggplot(cd14MonoMeta_clust1,aes(x = PASCphenotype,y = MC4AUC,fill = PASCphenotype))+ggtitle("Clust1")+
  geom_violin(scale = "width")+theme_classic()+geom_hline(yintercept = 0.08,linetype = "dashed")+
  theme(axis.text.x = element_text(size = 12,color = "black",angle = 45, hjust = 1,vjust = 1))+
  xlab("")+ylab("MC4 AUC") + scale_fill_manual(values = c("nonResp-PASC" = "cyan4", "Resp-PASC" = "coral1"))+
  stat_compare_means(comparisons = list(c("nonResp-PASC","Resp-PASC")),aes(label = after_stat(p.signif)))
