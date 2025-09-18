library(Seurat)
library(Signac)
library(DESeq2)

load("cd14MonoSubset_pseudobulk_deseq2Analysis.RData")

cd14Mono_pseudobulk = AggregateExpression(cd14Mono,assays = c("RNA","peaks"),group.by = "donorID_TimePoint")
cd14Mono_RNA_pseudobulk = cd14Mono_pseudobulk$RNA
cd14Mono_peaks_pseudobulk = cd14Mono_pseudobulk$peaks
colnames(cd14Mono_RNA_pseudobulk) = gsub("-","_",colnames(cd14Mono_RNA_pseudobulk))
colnames(cd14Mono_RNA_pseudobulk) = gsub("g","",colnames(cd14Mono_RNA_pseudobulk))
cd14Mono_RNA_pseudobulk = as.data.frame(cd14Mono_RNA_pseudobulk)
cd14Mono_RNA_pseudobulk = subset(cd14Mono_RNA_pseudobulk,!grepl("^AL[0-6].+",rownames(cd14Mono_RNA_pseudobulk)))
cd14Mono_RNA_pseudobulk = subset(cd14Mono_RNA_pseudobulk,!grepl("^AC[0-6].+",rownames(cd14Mono_RNA_pseudobulk)))
cd14Mono_RNA_pseudobulk = subset(cd14Mono_RNA_pseudobulk,!grepl("^AP[0-6].+",rownames(cd14Mono_RNA_pseudobulk)))
cd14Mono_RNA_pseudobulk = subset(cd14Mono_RNA_pseudobulk,!grepl("^MT-",rownames(cd14Mono_RNA_pseudobulk)))
cd14Mono_RNA_pseudobulk = subset(cd14Mono_RNA_pseudobulk,!grepl("IGH.+",rownames(cd14Mono_RNA_pseudobulk)))


metadata = cd14Mono@meta.data
metadataForPseudobulk = metadata[,c(4,5,17,22,23,24,31,33)]

rownames(metadataForPseudobulk) = NULL
metadataForPseudobulk = metadataForPseudobulk[!duplicated(metadataForPseudobulk),]
metadataForPseudobulk$cellCluster = NULL

metadataForPseudobulk = metadataForPseudobulk[match(colnames(cd14Mono_RNA_pseudobulk),metadataForPseudobulk$donorID_TimePoint),]

####### deseq2 analysis####
dds <- DESeqDataSetFromMatrix(countData = cd14Mono_RNA_pseudobulk,
                              colData = metadataForPseudobulk,
                              design= ~ Runs + TP_RecovLC)

smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds$Recov_LC <- factor(dds$Recov_LC, levels = c("HC_HC","nonICU_T1","nonICU_LC","ICU_T1","ICU_LC","R_Recov"))

dds$TP_RecovLC = factor(dds$TP_RecovLC,levels = c("R_Recov_R","HC_HC_HC","nonICU_T1_T1","nonICU_LC_T2","nonICU_LC_T3","nonICU_LC_T4","nonICU_LC_T5","ICU_T1_T1",
                                                    "ICU_LC_T2","ICU_LC_T3","ICU_LC_T4","ICU_LC_T5"))

dds <- DESeq(dds)

res_nonICU_T2vsT1 = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T2","nonICU_T1_T1"))
res_nonICU_T3vsT1 = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T3","nonICU_T1_T1"))
res_nonICU_T4vsT1 = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T4","nonICU_T1_T1"))
res_nonICU_T5vsT1 = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T5","nonICU_T1_T1"))


res_nonICU_T2vsRecov = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T2","R_Recov_R"))
res_nonICU_T3vsRecov = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T3","R_Recov_R"))
res_nonICU_T4vsRecov = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T4","R_Recov_R"))
res_nonICU_T5vsRecov = results(dds, contrast=c("TP_RecovLC","nonICU_LC_T5","R_Recov_R"))

res_ICU_T2vsT1 = results(dds, contrast=c("TP_RecovLC","ICU_LC_T2","ICU_T1_T1"))
res_ICU_T3vsT1 = results(dds, contrast=c("TP_RecovLC","ICU_LC_T3","ICU_T1_T1"))
res_ICU_T4vsT1 = results(dds, contrast=c("TP_RecovLC","ICU_LC_T4","ICU_T1_T1"))
res_ICU_T5vsT1 = results(dds, contrast=c("TP_RecovLC","ICU_LC_T5","ICU_T1_T1"))

res_ICU_T2vsRecov = results(dds, contrast=c("TP_RecovLC","ICU_LC_T2","R_Recov_R"))
res_ICU_T3vsRecov = results(dds, contrast=c("TP_RecovLC","ICU_LC_T3","R_Recov_R"))
res_ICU_T4vsRecov = results(dds, contrast=c("TP_RecovLC","ICU_LC_T4","R_Recov_R"))
res_ICU_T5vsRecov = results(dds, contrast=c("TP_RecovLC","ICU_LC_T5","R_Recov_R"))


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
fgsea_nonICU_T2vsT1 = enrichmentFunc(res_nonICU_T2vsT1)
fgsea_nonICU_T3vsT1 = enrichmentFunc(res_nonICU_T3vsT1)
fgsea_nonICU_T4vsT1 = enrichmentFunc(res_nonICU_T4vsT1)
fgsea_nonICU_T5vsT1 = enrichmentFunc(res_nonICU_T5vsT1)

fgsea_nonICU_T2vsRecov = enrichmentFunc(res_nonICU_T2vsRecov)
fgsea_nonICU_T3vsRecov = enrichmentFunc(res_nonICU_T3vsRecov)
fgsea_nonICU_T4vsRecov = enrichmentFunc(res_nonICU_T4vsRecov)
fgsea_nonICU_T5vsRecov = enrichmentFunc(res_nonICU_T5vsRecov)


fgsea_ICU_T2vsT1 = enrichmentFunc(res_ICU_T2vsT1)
fgsea_ICU_T3vsT1 = enrichmentFunc(res_ICU_T3vsT1)
fgsea_ICU_T4vsT1 = enrichmentFunc(res_ICU_T4vsT1)
fgsea_ICU_T5vsT1 = enrichmentFunc(res_ICU_T5vsT1)

fgsea_ICU_T2vsRecov = enrichmentFunc(res_ICU_T2vsRecov)
fgsea_ICU_T3vsRecov = enrichmentFunc(res_ICU_T3vsRecov)
fgsea_ICU_T4vsRecov = enrichmentFunc(res_ICU_T4vsRecov)
fgsea_ICU_T5vsRecov = enrichmentFunc(res_ICU_T5vsRecov)
