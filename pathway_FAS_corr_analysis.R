library(Seurat)
library(Signac)
library(AUCell)
library(GSEABase)

cd8T = subset(bigObj,subset = CellType == "CD8 T cells")
nkCells = subset(bigObj,subset = CellType == "NK cells")
cd14Mono = subset(bigObj,subset = CellType == "CD14 Monocytes")

rm(bigObj)

cd8T_pb = AggregateExpression(cd8T,assays = "RNA",group.by = "donorID_TimePoint")
nk_pb = AggregateExpression(nkCells,assays = "RNA",group.by = "donorID_TimePoint")
Mono_pb = AggregateExpression(cd14Mono,assays = "RNA",group.by = "donorID_TimePoint")

cd8T_pb = cd8T_pb$RNA
nk_pb = nk_pb$RNA
Mono_pb = Mono_pb$RNA


REACTOME_genesets = getGmt('/home/skumar/bulkRNAseq-influenza/analysis/filesForGSEA/c2.cp.reactome.v7.5.1.symbols.gmt')
reactPathways = names(REACTOME_genesets)
reactPathwaysTOKeep = reactPathways[grepl("TOLL",reactPathways)]
reactPathwaysTOKeep = append(reactPathwaysTOKeep,reactPathways[grepl("ZBP1",reactPathways)])
reactPathwaysTOKeep = append(reactPathwaysTOKeep,reactPathways[grepl("MYD88_INDEPENDENT_TLR4",reactPathways)])
reactPathwaysTOKeep = append(reactPathwaysTOKeep,reactPathways[grepl("INTERFERON_ALPHA_BETA",reactPathways)])
reactPathwaysTOKeep = append(reactPathwaysTOKeep,reactPathways[grepl("SIGNALING_BY_WNT",reactPathways)])
reactPathwaysTOKeep = reactPathwaysTOKeep[c(1:4,6:7,9)]
REACTOME_genesets = subset(REACTOME_genesets,names(REACTOME_genesets) %in% reactPathwaysTOKeep)

HALLMARK_genesets = getGmt("/vol/projects/skumar/LongCovid/freeze2Data_workingAnalysis2/h.all.v2022.1.Hs.symbols.gmt")
HM_pathways = names(HALLMARK_genesets)
HM_pathwaysTOKeep = HM_pathways[grepl("TGF_BETA",HM_pathways)]
HM_pathwaysTOKeep = append(HM_pathwaysTOKeep,HM_pathways[grepl("TNFA",HM_pathways)])
HM_pathwaysTOKeep = append(HM_pathwaysTOKeep,HM_pathways[grepl("WNT",HM_pathways)])
HM_pathwaysTOKeep = append(HM_pathwaysTOKeep,HM_pathways[grepl("HALLMARK_NOTCH_SIGNALING",HM_pathways)])
HM_pathwaysTOKeep = append(HM_pathwaysTOKeep,HM_pathways[grepl("HALLMARK_IL2_STAT5_SIGNALING",HM_pathways)])
HALLMARK_genesets = subset(HALLMARK_genesets,names(HALLMARK_genesets) %in% HM_pathwaysTOKeep)


pathwaysAUC = function(pseudobulkFile,cellType){
  pseudobulkFile_ranked = AUCell_buildRankings(pseudobulkFile)
  REACTOME_genesets = subsetGeneSets(REACTOME_genesets,rownames(pseudobulkFile))
  pb_REACT_geneSetsAUC = AUCell_calcAUC(REACTOME_genesets,pseudobulkFile_ranked)
  pb_REACT_geneSetsAUCScore = t(getAUC(pb_REACT_geneSetsAUC))
  pb_REACT_geneSetsAUCScore = as.data.frame(pb_REACT_geneSetsAUCScore)
  
  HALLMARK_genesets = subsetGeneSets(HALLMARK_genesets,rownames(pseudobulkFile))
  pb_HM_geneSetsAUC = AUCell_calcAUC(HALLMARK_genesets,pseudobulkFile_ranked)
  pb_HM_geneSetsAUCScore = t(getAUC(pb_HM_geneSetsAUC))
  pb_HM_geneSetsAUCScore = as.data.frame(pb_HM_geneSetsAUCScore)
  pb_AUCscoreFile = cbind(pb_HM_geneSetsAUCScore,pb_REACT_geneSetsAUCScore)
  rownames(pb_AUCscoreFile) = gsub("g","",rownames(pb_AUCscoreFile))
  rownames(pb_AUCscoreFile) = gsub("-","_",rownames(pb_AUCscoreFile))
  metadataFile = cellType@meta.data
  rownames(metadataFile) = NULL
  metadataFile = metadataFile[,c(23,20,32,34)]
  metadataFile = metadataFile[!duplicated(metadataFile),]
  mergedFile = merge(metadataFile,pb_AUCscoreFile,by.x="donorID_TimePoint",by.y=0)
  mergedFile2 = merge(phenotypeFile,mergedFile,by.x="PatientID_TP",by.y="donorID_TimePoint")
  mergeFileMelted= melt(mergedFile2,id.vars = colnames(mergedFile2)[1:23])
  return(mergeFileMelted)
}

cd8TFileFinal = pathwaysAUC(cd8T_pb,cd8T)
nkFileFinal = pathwaysAUC(nk_pb,nkCells)
monoFileFinal = pathwaysAUC(Mono_pb,cd14Mono)
