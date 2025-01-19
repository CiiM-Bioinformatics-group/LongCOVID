library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(future)

plan("multicore", workers = 6)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

########calculating all dpa, dtfa within cd14 Mono clusters####
Idents(cd14Mono) = "cellCluster"
DefaultAssay(cd14Mono) = "chromvar"
tfActivity <- FindAllMarkers(cd14Mono,mean.fxn = rowMeans,fc.name = "avg_diff")
tfActivity_signi = subset(tfActivity, tfActivity$p_val_adj < 0.05)


Idents(cd14Mono) = "cellCluster"
DefaultAssay(cd14Mono) = "peaks"
clusterPeaks = FindAllMarkers(cd14Mono,only.pos = TRUE,test.use = 'LR',latent.vars = 'nCount_peaks')
clusterPeaks_signi = subset(clusterPeaks, clusterPeaks$p_val_adj < 0.05)
motifNames = ConvertMotifID(cd14Mono,id = tfActivity_signi$gene)
tfActivity_signi$motifname = motifNames

######getting linked peaks genes matrix#####
DefaultAssay(bigObj) = "peaks"
links = Links(bigObj)
library(GenomicRanges)
linksDF = as.data.frame(links)

Idents(cd14Mono) = "cellCluster"
DefaultAssay(cd14Mono) = "chromvar"
MC4_tfActivity <- FindMarkers(cd14Mono,ident.1 = "MC4",ident.2 = c("MC1","MC2","MC3","MC4"),mean.fxn = rowMeans,fc.name = "avg_diff")
MC4_tfActivity = subset(MC4_tfActivity, MC4_tfActivity$p_val_adj < 0.05)

DefaultAssay(cd14Mono) = "peaks"
motifNames = ConvertMotifID(cd14Mono,id = MC4_tfActivity$gene)
MC4_tfActivity$MotifName = motifNames

MC4Peaks = FindMarkers(cd14Mono,ident.1 = "MC4",ident.2 = c("MC1","MC2","MC3"),
                       only.pos = TRUE,test.use = 'LR',min.pct = 0.05,latent.vars = "nCount_peaks")

MC4Peak_linkedGenes = subset(linksDF,linksDF$peak %in% MC4Peaks$gene)

MC4Peak_linkedGenes_MC4_DEgenes = subset(MC4_dge,MC4_dge$gene %in% MC4Peak_linkedGenes$gene)

####checking if motif is present in the peaks that are linked to genes for MC4 only . 
#########the motifs used are those that are significant for MC4 based on chromvar avg_diff > 0.2
library(Signac)
library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(cd14Mono),
  pwm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

topMC4Peaks = rownames(MC4Peaks[MC4Peaks$p_val_adj < 0.005 & MC4Peaks$pct.1 > 0.1, ])
# find peaks open in MC4 cluster
open.peaks <- AccessiblePeaks(cd14Mono, idents = "MC4")
topMC4PeaksToKeep = intersect(topMC4Peaks,open.peaks)
# match the overall GC content in the peak set
meta.feature <- GetAssayData(cd14Mono, assay = "peaks", layer = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[topMC4PeaksToKeep, ],
  n = 50000)

# test enrichment
MC4enriched.motifs <- FindMotifs(
  object = cd14Mono,
  features =topMC4PeaksToKeep, background = peaks.matched
)

#####correlating gene expression of TFs and genes regulating ###
topMC4Peaks = rownames(MC4Peaks_signi[MC4Peaks_signi$pct.1 > 0.1, ])
open.peaks <- AccessiblePeaks(cd14Mono, idents = "MC4")
peakLinkedGenes_down = c("TNFRSF1B","SLC11A1","NACC2","ARAP1","FCHSD2","GRN","LDLRAD4","MAFB","SULF2","TYMP","ODF3B")

topMC4PeaksToKeep = Reduce(intersect,list(topMC4Peaks,linksDF$peak,open.peaks))
topMC4PeakstoKeep_linkedGenes = subset(linksDF,linksDF$peak %in% topMC4PeaksToKeep)
topMC4PeakstoKeep_linkedGenes = subset(topMC4PeakstoKeep_linkedGenes,!topMC4PeakstoKeep_linkedGenes$gene %in% peakLinkedGenes_down)

topMC4PeaksToKeep = unique(topMC4PeakstoKeep_linkedGenes$peak)
# match the overall GC content in the peak set
meta.feature <- GetAssayData(cd14Mono, assay = "peaks", layer = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[topMC4PeaksToKeep, ],
  n = 50000)

# test enrichment
MC4enriched.motifs <- FindMotifs(
  object = cd14Mono,
  features =topMC4PeaksToKeep, background = peaks.matched
)

MC4enriched.motifs_signi = subset(MC4enriched.motifs,MC4enriched.motifs$p.adjust < 0.05)

MC4enriched.motifs_merged = merge(MC4enriched.motifs,MC4_tfActivity,by.x= "motif",by.y = "gene")
MC4enriched.motifs_merged$motif.name = gsub("JUN\\(.*","JUN",MC4enriched.motifs_merged$motif.name)
MC4enriched.motifs_merged = subset(MC4enriched.motifs_merged,!MC4enriched.motifs_merged$motif == "MA0488.1")
MC4enriched.motifs_merged$neglog10PvalEnriched = -(log10(MC4enriched.motifs_merged$p.adjust)) 


MC4_peaksEnrichedTFs = MC4enriched.motifs_signi$motif.name
MC4_peaksEnrichedTFs = MC4_peaksEnrichedTFs[!grepl("::",MC4_peaksEnrichedTFs)]
MC4_peaksEnrichedTFs = gsub("\\(var.*","",MC4_peaksEnrichedTFs)
MC4_peaksEnrichedTFs = append(MC4_peaksEnrichedTFs,c("SMAD2","SMAD3","SMAD4","MAF","MYC","MAX","NFIX","NR2F1",
                                                     "TFAP2B","TLX1","JUN","TFAP2C","TFAP2A","NEUROG2","THRB",
                                                     "ATOH1","MZF1","NFIC","RXRB","RXRG","RARA","BHLHA15","RARG"
                                                     ,"RARA"))
MC4_peaksEnrichedTFs = unique(MC4_peaksEnrichedTFs)

MC4_tfActivity_subset = subset(MC4_tfActivity, MC4_tfActivity$avg_diff > 0.2)
MC4enriched.motifs_signi_subset = subset(MC4enriched.motifs_signi,MC4enriched.motifs_signi$fold.enrichment > 2)
AllMotifsIDsToConsider = union(MC4_tfActivity_subset$motifID, MC4enriched.motifs_signi_subset$motif)

######## getting TF motifs in the peaks accessible to linked genes for MC4
motifs_peaks = motif.matrix[topMC4PeaksToKeep,AllMotifsIDsToConsider]

motifs = list()
for(i in colnames(motifs_peaks)){
  print(i)
  allPeaks = motifs_peaks[,i]
  truePeaks = subset(allPeaks,allPeaks == TRUE)
  motifs[[i]] = truePeaks
}

peaks = list()
for(i in rownames(motifs_peaks)){
  print(i)
  allPeaks = motifs_peaks[i,]
  #print(allPeaks)
  truePeaks = subset(allPeaks,allPeaks == TRUE)
  peaks[[i]] = truePeaks
}

peaks_motif_df = list()
for(i in names(peaks)){
  print(i)
  mots = names(peaks[[i]])
  if(length(mots) > 0){
    mots = as.data.frame(mots)
    mots$peak = i
    peaks_motif_df[[i]] = mots
  }
}
allFound = bind_rows(peaks_motif_df)
allFound_WithGenesLinked = merge(allFound,topMC4PeakstoKeep_linkedGenes,by = "peak")
allFound_WithGenesLinked$MotifName = ConvertMotifID(cd14Mono,id = allFound_WithGenesLinked$mots)
motifsTouseForGeneExp = unique(allFound_WithGenesLinked$MotifName)
motifsTouseForGeneExp =gsub("\\(var.*","",motifsTouseForGeneExp)
motifswithdelim = motifsTouseForGeneExp[grepl("::",motifsTouseForGeneExp)]
motifswithdelim = unlist(strsplit(motifswithdelim,"::"))
motifsTouseForGeneExp = unique(append(motifsTouseForGeneExp,motifswithdelim))
motifsTouseForGeneExp = motifsTouseForGeneExp[!grepl("::",motifsTouseForGeneExp)]

allGenesForGeneMat = unique(c(motifsTouseForGeneExp,topMC4PeakstoKeep_linkedGenes$gene))

motifsGeneMat = subset(geneMat, rownames(geneMat) %in% motifsTouseForGeneExp)

allGenesGeneMat = subset(geneMat, rownames(geneMat) %in% allGenesForGeneMat)

motifsGeneMat_MC4CellsOnly = motifsGeneMat[,rownames(MC4CellIDs)]
allGenesGeneMat_MC4CellsOnly = allGenesGeneMat[,rownames(MC4CellIDs)]

motifsGeneMat_MC4CellsOnly = t(motifsGeneMat_MC4CellsOnly)
identical(rownames(motifsGeneMat_MC4CellsOnly),rownames(genesWithLinkedPeaks_MC4Mat_MC4CellsOnly))

motif_gene_cor = corSparse(motifsGeneMat_MC4CellsOnly,genesWithLinkedPeaks_MC4Mat_MC4CellsOnly)

colnames(motif_gene_cor) = colnames(genesWithLinkedPeaks_MC4Mat_MC4CellsOnly)
rownames(motif_gene_cor) = colnames(motifsGeneMat_MC4CellsOnly)




#####chromvar activity compared to recovered.


DefaultAssay(cd14Mono) = "chromvar"
Idents(cd14Mono) = "AC_RecovStatus"

cd14Mono_nonICU_T2_Recov_motif = FindMarkers(cd14Mono,ident.1 = "nonICU_T2", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")
cd14Mono_nonICU_T3_Recov_motif = FindMarkers(cd14Mono,ident.1 = "nonICU_T3", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")
cd14Mono_nonICU_T4_Recov_motif = FindMarkers(cd14Mono,ident.1 = "nonICU_T4", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")
cd14Mono_nonICU_T5_Recov_motif = FindMarkers(cd14Mono,ident.1 = "nonICU_T5", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")

cd14Mono_ICU_T2_Recov = FindMarkers(cd14Mono,ident.1 = "ICU_T2", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")
cd14Mono_ICU_T3_Recov = FindMarkers(cd14Mono,ident.1 = "ICU_T3", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")
cd14Mono_ICU_T4_Recov = FindMarkers(cd14Mono,ident.1 = "ICU_T4", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")
cd14Mono_ICU_T5_Recov = FindMarkers(cd14Mono,ident.1 = "ICU_T5", ident.2 = "Recovered",mean.fxn = rowMeans,fc.name = "avg_diff")

DefaultAssay(cd14Mono) = "peaks"

cd14Mono_nonICU_T2_Recov_motif$motifName = ConvertMotifID(cd14Mono,id = rownames(cd14Mono_nonICU_T2_Recov_motif))
cd14Mono_nonICU_T3_Recov_motif$motifName = ConvertMotifID(cd14Mono,id = rownames(cd14Mono_nonICU_T3_Recov_motif))
cd14Mono_nonICU_T4_Recov_motif$motifName = ConvertMotifID(cd14Mono,id = rownames(cd14Mono_nonICU_T4_Recov_motif))
cd14Mono_nonICU_T5_Recov_motif$motifName = ConvertMotifID(cd14Mono,id = rownames(cd14Mono_nonICU_T5_Recov_motif))

cd14Mono_nonICU_T2_Recov_motif_signi = subset(cd14Mono_nonICU_T2_Recov_motif, cd14Mono_nonICU_T2_Recov_motif$p_val_adj < 0.05)
cd14Mono_nonICU_T3_Recov_motif_signi = subset(cd14Mono_nonICU_T3_Recov_motif, cd14Mono_nonICU_T3_Recov_motif$p_val_adj < 0.05)
cd14Mono_nonICU_T4_Recov_motif_signi = subset(cd14Mono_nonICU_T4_Recov_motif, cd14Mono_nonICU_T4_Recov_motif$p_val_adj < 0.05)

cd14Mono_nonICU_T2_Recov_motif_signi$Comparison = "T2vsR"
cd14Mono_nonICU_T3_Recov_motif_signi$Comparison = "T3vsR"
cd14Mono_nonICU_T4_Recov_motif_signi$Comparison = "T4vsR"

cd14Mono_nonIC_vsR_motifs = rbind(cd14Mono_nonICU_T2_Recov_motif_signi,cd14Mono_nonICU_T3_Recov_motif_signi,cd14Mono_nonICU_T4_Recov_motif_signi)
cd14Mono_nonIC_vsR_motifs$negLog10Padj = -(log10(cd14Mono_nonIC_vsR_motifs$p_val_adj))
cd14Mono_nonIC_vsR_motifs$CellType = "CD14 Mono"
