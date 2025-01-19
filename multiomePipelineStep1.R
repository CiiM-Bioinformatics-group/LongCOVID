rm(list = ls())

library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(future)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(SeuratData)
#library(tidyr)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

######## this pipeline creates R object only for those donors which are demultiplexed in both RNA and ATAC####
atac_cluster_donor_batch3 = read.table("Batch3freeze2Atac.txt",header = T)


sampleIDs = unique(atac_cluster_donor_batch3$pool)
batch2 = c("pool3A","pool3B","pool4B","pool4C")
sampleIDsAll =  append(sampleIDs,c("pool3A","pool3B","pool4B","pool4C"))

##### creating granges and common peaks #####

gr_list = list()
for (i in sampleIDs) {

  peaks_bed = paste0("/batch3/Multiome/", i,"/outs/atac_peaks.bed")
  peaksFile = read.table(file = peaks_bed,col.names = c("chr", "start", "end"))
  grFile = makeGRangesFromDataFrame(peaksFile)
  gr_list[[i]] = grFile
}

for (i in batch2){
  peaks_bed = paste0("/batch2/", i,"/outs/atac_peaks.bed")
  peaksFile = read.table(file = peaks_bed,col.names = c("chr", "start", "end"))
  grFile = makeGRangesFromDataFrame(peaksFile)
  gr_list[[i]] = grFile
}


###### add all the granges objects of all the pools to reduce to 1 common peaks object ##### this includes both batches #######
combined.peaks = reduce(x = c(gr_list$pool51,gr_list$pool52,gr_list$pool53,gr_list$pool54,gr_list$pool55,gr_list$pool56,gr_list$pool5A,
                              gr_list$pool5B,gr_list$pool71,gr_list$pool73,gr_list$pool74,gr_list$pool75,gr_list$poolCV1,gr_list$poolCV2,
                              gr_list$poolCV3,gr_list$poolCV4,gr_list$poolCV5,gr_list$poolCV6,
                              gr_list$poolF1,gr_list$poolF2,gr_list$poolF3,gr_list$poolF4,gr_list$poolF5,gr_list$poolF6,
                              gr_list$poolFA,gr_list$poolFB,gr_list$poolS1,gr_list$poolS2,gr_list$poolS3,gr_list$poolS4,gr_list$poolS5,gr_list$poolS6,gr_list$poolT1,
                              gr_list$poolT2,gr_list$poolT3,gr_list$poolT4,gr_list$poolT5,gr_list$poolT6,gr_list$poolT7,gr_list$poolTA,
                              gr_list$poolTB,gr_list$poolTC,gr_list$pool3A,gr_list$pool3B,gr_list$pool4B,gr_list$pool4C))

peakwidths = width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
seur_list = list()
i <- 'pool51'
for (i in sampleIDs) {

  print(paste0('Processing pool: ', i))

  # Reading in singlets pre-processed file#
  singlets <- read.csv(paste0('cluster_donorInfor/', i, '_common_singlets.csv'),header = T)
  rownames(singlets) = singlets$cellID
  # Reading in 10X data which is a list of two data matrices: GEX and ATAC
  data <- Read10X_h5(paste0("/batch3/Multiome/", i, "/outs/filtered_feature_bc_matrix.h5"))
  gex <- data$`Gene Expression`
  atac <- data$Peaks
  gex = gex[,colnames(gex) %in% singlets$cellID]
  dim(gex)

  seur <- CreateSeuratObject(gex,
                             assay = 'RNA',
                             min.cells = 5,
                             meta.data = singlets)
  seur[["percent.mt"]] =  PercentageFeatureSet(seur, pattern = "^MT-")
  seur[["percent.ribo"]] = PercentageFeatureSet(seur,pattern = "^RP[SL][0-9]+$")


  # Same singletons as in GEX data
  atac= atac[,colnames(atac) %in% singlets$cellID]
  stopifnot(all(colnames(atac) == colnames(gex)))

  ########creating atac assay #####
  frags <- CreateFragmentObject(paste0("/batch3/Multiome/", i,"/outs/atac_fragments.tsv.gz"),cells = colnames(atac))
  counts <- FeatureMatrix(fragments = frags,features = combined.peaks,cells = colnames(atac))
  chromAssay <- CreateChromatinAssay(counts, fragments = frags)
  seur[["ATAC"]] <- chromAssay
  DefaultAssay(seur) <- "ATAC"
  Annotation(seur) = annotations
  seur = NucleosomeSignal(seur)
  seur = TSSEnrichment(seur)
  DefaultAssay(seur) <- "RNA"
  seur_list[[i]] <- seur

}
for (i in batch2) {
  print(paste0('Processing pool: ', i))

  # Reading in singlets pre-processed file#
  singlets <- read.csv(paste0('cluster_donorInfor/', i, '_common_singlets.csv'),header = T)
  rownames(singlets) = singlets$cellID
  # Reading in 10X data which is a list of two data matrices: GEX and ATAC
  data <- Read10X_h5(paste0("/batch2/", i, "/outs/filtered_feature_bc_matrix.h5"))
  gex <- data$`Gene Expression`
  atac <- data$Peaks
  gex = gex[,colnames(gex) %in% singlets$cellID]
  dim(gex)

  seur <- CreateSeuratObject(gex,
                             assay = 'RNA',
                             min.cells = 5,
                             meta.data = singlets)
  seur[["percent.mt"]] =  PercentageFeatureSet(seur, pattern = "^MT-")
  seur[["percent.ribo"]] = PercentageFeatureSet(seur,pattern = "^RP[SL][0-9]+$")


  # Same singletons as in GEX data
  atac= atac[,colnames(atac) %in% singlets$cellID]
  stopifnot(all(colnames(atac) == colnames(gex)))

  ######## creating atac assay #####
  frags <- CreateFragmentObject(paste0("/batch2/", i,"/outs/atac_fragments.tsv.gz"),cells = colnames(atac))
  counts <- FeatureMatrix(fragments = frags,features = combined.peaks,cells = colnames(atac))
  chromAssay <- CreateChromatinAssay(counts, fragments = frags)
  seur[["ATAC"]] <- chromAssay
  DefaultAssay(seur) <- "ATAC"
  Annotation(seur) = annotations
  seur = NucleosomeSignal(seur)
  seur = TSSEnrichment(seur)
  #VlnPlot(object = seur,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)
  #ggsave(paste0("/vol/projects/skumar/LongCovid/freeze2Data/",i,"_VlnPlot2.png"),width = 50,height = 30, units = "cm")
  DefaultAssay(seur) <- "RNA"
  seur_list[[i]] <- seur
}

########merging seurat object list to for QC #######

mergedSeurObj = merge(x=seur_list$pool51,y=c(seur_list$pool52,seur_list$pool53,seur_list$pool54,seur_list$pool55,seur_list$pool56,seur_list$pool5A,
                                             seur_list$pool5B,seur_list$pool71,seur_list$pool73,seur_list$pool74,seur_list$pool75,seur_list$poolCV1,seur_list$poolCV2,
                                             seur_list$poolCV3,seur_list$poolCV4,seur_list$poolCV5,seur_list$poolCV6,
                                             seur_list$poolF1,seur_list$poolF2,seur_list$poolF3,seur_list$poolF4,seur_list$poolF5,seur_list$poolF6,
                                             seur_list$poolFA,seur_list$poolFB,seur_list$poolS1,seur_list$poolS2,seur_list$poolS3,seur_list$poolS4,seur_list$poolS5,seur_list$poolS6,seur_list$poolT1,
                                             seur_list$poolT2,seur_list$poolT3,seur_list$poolT4,seur_list$poolT5,seur_list$poolT6,seur_list$poolT7,seur_list$poolTA,
                                             seur_list$poolTB,seur_list$poolTC,seur_list$pool3A,seur_list$pool3B,seur_list$pool4B,seur_list$pool4C),add.cell.ids = sampleIDsAll)


saveRDS(mergedSeurObj,file = "/merged_unprocessedObj.rds")
