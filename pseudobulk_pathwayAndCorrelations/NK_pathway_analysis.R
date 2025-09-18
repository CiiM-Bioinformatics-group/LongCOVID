library(Seurat)
library(Signac)
library(DESeq2)

nk_pseudobulk = AggregateExpression(nkcells,assays = "RNA",group.by = "donorID_TimePoint")
nk_RNA_pseudobulk = nk_pseudobulk$RNA
colnames(nk_RNA_pseudobulk) = gsub("-","_",colnames(nk_RNA_pseudobulk))
colnames(nk_RNA_pseudobulk) = gsub("g","",colnames(nk_RNA_pseudobulk))
nk_RNA_pseudobulk = as.data.frame(nk_RNA_pseudobulk)
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("^AL[0-6].+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("^AC[0-6].+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("^AP[0-6].+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("^MT-",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("IGH.+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("IGK.+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("IGL.+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("RPL.+",rownames(nk_RNA_pseudobulk)))
nk_RNA_pseudobulk = subset(nk_RNA_pseudobulk,!grepl("RPS.+",rownames(nk_RNA_pseudobulk)))


metadata = nkcells@meta.data
metadataForPseudobulk = metadata[,c(4,5,17,23,24,30)]

rownames(metadataForPseudobulk) = NULL
metadataForPseudobulk = metadataForPseudobulk[!duplicated(metadataForPseudobulk),]
metadataForPseudobulk$cellCluster = NULL

metadataForPseudobulk = metadataForPseudobulk[match(colnames(nk_RNA_pseudobulk),metadataForPseudobulk$donorID_TimePoint),]

####### deseq2 analysis####
dds <- DESeqDataSetFromMatrix(countData = nk_RNA_pseudobulk,
                              colData = metadataForPseudobulk,
                              design= ~ Runs + AcuteCovidStatus_TimePoint2)

smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$AcuteCovidStatus_TimePoint2 <- factor(dds$AcuteCovidStatus_TimePoint2, levels = c("HC_HC","nonICU_T1","nonICU_T2","nonICU_T3","nonICU_T4",
                                                                                      "nonICU_T5","ICU_T1","ICU_T2","ICU_T3","ICU_T4","ICU_T5","R_Recov"))
dds <- DESeq(dds)

res_nk_nonICU_T2vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T2","nonICU_T1"))
res_nk_nonICU_T3vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T3","nonICU_T1"))
res_nk_nonICU_T4vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T4","nonICU_T1"))
res_nk_nonICU_T5vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T5","nonICU_T1"))


res_nk_nonICU_T2vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T2","R_Recov"))
res_nk_nonICU_T3vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T3","R_Recov"))
res_nk_nonICU_T4vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T4","R_Recov"))
res_nk_nonICU_T5vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","nonICU_T5","R_Recov"))


res_nk_ICU_T2vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T2","ICU_T1"))
res_nk_ICU_T3vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T3","ICU_T1"))
res_nk_ICU_T4vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T4","ICU_T1"))
res_nk_ICU_T5vsT1 = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T5","ICU_T1"))


res_nk_ICU_T2vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T2","R_Recov"))
res_nk_ICU_T3vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T3","R_Recov"))
res_nk_ICU_T4vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T4","R_Recov"))
res_nk_ICU_T5vsRecov = results(dds, contrast=c("AcuteCovidStatus_TimePoint2","ICU_T5","R_Recov"))


library(fgsea)
fgsea::gmtPathways('/home/skumar/bulkRNAseq-influenza/analysis/filesForGSEA/c2.cp.reactome.v7.5.1.symbols.gmt') -> genesets
fgsea::gmtPathways("/vol/projects/skumar/LongCovid/freeze2Data_workingAnalysis2/h.all.v2022.1.Hs.symbols.gmt") -> genesets_hallmark

enrichmentFunc = function(resFile){
  resFile = as.data.frame(resFile)
  resFile$GeneName =rownames(resFile)
  resFile$StatCalc = sign(resFile$log2FoldChange)*(-(log10(resFile$pvalue)))
  ranked_list = resFile$StatCalc
  names(ranked_list) = resFile$GeneName
  ranked_list = na.omit(ranked_list)
  ranked_list <- ranked_list[!duplicated(names(ranked_list))]
  sort(ranked_list, decreasing = T) -> ranked_list
  fgseaRes_reactome <- fgsea(pathways = genesets, 
                             stats    = ranked_list, nPermSimple = 5000)
  
  fgseaRes_hallmark <- fgsea(pathways = genesets_hallmark, 
                             stats    = ranked_list, nPermSimple = 5000)
  fgseaBoth = rbind(fgseaRes_reactome,fgseaRes_hallmark)
  fgseaBoth = subset(fgseaBoth, fgseaBoth$padj < 0.1)
  return(fgseaBoth)
}
# fgsea_nonICU_T1vsHC = enrichmentFunc(res_nonICU_T1_HC)
 fgsea_nk_nonICU_T2vsT1 = enrichmentFunc(res_nk_nonICU_T2vsT1)
 fgsea_nk_nonICU_T3vsT1 = enrichmentFunc(res_nk_nonICU_T3vsT1)
 fgsea_nk_nonICU_T4vsT1 = enrichmentFunc(res_nk_nonICU_T4vsT1)
 fgsea_nk_nonICU_T5vsT1 = enrichmentFunc(res_nk_nonICU_T5vsT1)


fgsea_nk_nonICU_T2vsRecov = enrichmentFunc(res_nk_nonICU_T2vsRecov)
fgsea_nk_nonICU_T3vsRecov = enrichmentFunc(res_nk_nonICU_T3vsRecov)
fgsea_nk_nonICU_T4vsRecov = enrichmentFunc(res_nk_nonICU_T4vsRecov)
fgsea_nk_nonICU_T5vsRecov = enrichmentFunc(res_nk_nonICU_T5vsRecov)

fgsea_nk_ICU_T2vsT1 = enrichmentFunc(res_nk_ICU_T2vsT1)
fgsea_nk_ICU_T3vsT1 = enrichmentFunc(res_nk_ICU_T3vsT1)
fgsea_nk_ICU_T4vsT1 = enrichmentFunc(res_nk_ICU_T4vsT1)
fgsea_nk_ICU_T5vsT1 = enrichmentFunc(res_nk_ICU_T5vsT1)


fgsea_nk_ICU_T2vsRecov = enrichmentFunc(res_nk_ICU_T2vsRecov)
fgsea_nk_ICU_T3vsRecov = enrichmentFunc(res_nk_ICU_T3vsRecov)
fgsea_nk_ICU_T4vsRecov = enrichmentFunc(res_nk_ICU_T4vsRecov)
fgsea_nk_ICU_T5vsRecov = enrichmentFunc(res_nk_ICU_T5vsRecov)

