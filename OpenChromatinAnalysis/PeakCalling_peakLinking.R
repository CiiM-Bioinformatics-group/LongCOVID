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
library(tidyr)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)



DefaultAssay(bigObj) = "ATAC"
plan("multicore", workers = 6)
options(future.globals.maxSize = 100 * 1024 ^ 3) # for 100 Gb RAM
####### running with macs2 conda environment activated.
macs2peaks <- CallPeaks(object = bigObj,group.by = "CellType",macs2.path="miniconda3/envs/macs3/bin/macs3")
macs2peaks <- keepStandardChromosomes(macs2peaks, pruning.mode = "coarse")
macs2peaks <- subsetByOverlaps(x = macs2peaks, ranges = blacklist_hg38_unified, invert = TRUE)

######## creating counts matrix from macs2 peaks called #######
macs2_counts = FeatureMatrix(fragments = Fragments(bigObj),features = macs2peaks,cells = colnames(bigObj))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(bigObj) = annotations
bigObj[["peaks"]] = CreateChromatinAssay(counts = macs2_counts, fragments = Fragments(bigObj),annotation = annotations)
DefaultAssay(bigObj) <- "peaks"
bigObj <- RunTFIDF(bigObj)
bigObj <- FindTopFeatures(bigObj, min.cutoff = 20)
### linking peaks
bigObj =  LinkPeaks(bigObj,peak.assay = "peaks",expression.assay = "RNA")