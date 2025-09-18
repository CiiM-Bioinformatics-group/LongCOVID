library(Seurat)
library(Signac)

DefaultAssay(bigObj_withStim) = "RNA"
monocytes = subset(bigObj_withStim,subset = CellTypeAnnot == "Monocytes" & TimePoint %in% c("T2","T3","T4","T5") & Stimulation %in% c("C0","C2"))

monocytes[["RNA"]] <- split(monocytes[["RNA"]], f = monocytes$poolID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
monocytes = subset(monocytes, subset = poolID %in% c("pool3B","pool4B","pool51","pool54","pool55","pol56","pool5A",
                                                     "pool5B","pool71","pool73","pool74","pool75","poolF1","poolF3",
                                                     "poolF4","poolF6","poolS1","poolS2","poolS3","poolS4","poolS5","poolS6",
                                                     "poolS6","poolT1","poolT2","poolT5","poolT6","poolT7","poolTA",
                                                     "poolTB","poolTC"))
monocytes <- NormalizeData(monocytes)
monocytes <- FindVariableFeatures(monocytes)
monocytes <- ScaleData(monocytes)
monocytes <- RunPCA(monocytes)
monocytes <- IntegrateLayers(object = monocytes, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 50)

# re-join layers after integration
monocytes[["RNA"]] <- JoinLayers(monocytes[["RNA"]])

monocytes <- FindNeighbors(monocytes, reduction = "integrated.rpca", dims = 1:10)
monocytes <- FindClusters(monocytes, resolution = 0.2)
monocytes <- RunUMAP(monocytes, dims = 1:10, reduction = "integrated.rpca")
monocytesMarkers = FindAllMarkers(monocytes,only.pos = T)
monocytesMarkers_signi = subset(monocytesMarkers, monocytesMarkers$p_val_adj < 0.05 & monocytesMarkers$pct.1 > 0.1)

monocytes$CellTypeAnnot = ifelse(monocytes$RNA_snn_res.0.2 %in% c(0,1,4), "CD14", ifelse(monocytes$RNA_snn_res.0.2 == 3, "CD16",ifelse(monocytes$RNA_snn_res.0.2 == 5, "moDC", "CD8Amix")))

cd14MonoOnly = subset(monocytes, subset = CellTypeAnnot == "CD14")
metdata = cd14MonoOnly@meta.data
export = CreateSeuratObject(counts = LayerData(cd14MonoOnly,assay = "RNA",layer = "counts"),project = "Export",meta.data = metdata)
export[["RNA"]] <- split(export[["RNA"]], f = export$poolID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
export <- NormalizeData(export)
export <- FindVariableFeatures(export)
export <- ScaleData(export)
export <- RunPCA(export)
export <- IntegrateLayers(object = export, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 40)

# re-join layers after integration
export[["RNA"]] <- JoinLayers(export[["RNA"]])

export <- FindNeighbors(export, reduction = "integrated.rpca", dims = 1:10)
export <- FindClusters(export, resolution = 0.2)
export <- RunUMAP(export, dims = 1:10, reduction = "integrated.rpca")
exportMarkers = FindAllMarkers(export,only.pos = T)
exportMarkers_signi = subset(exportMarkers, exportMarkers$p_val_adj < 0.05 & exportMarkers$pct.1 > 0.10 )

MC4vsother_stim_genes2 <- FindMarkers(object = export,ident.1 = 'MC4Donor_C2',ident.2 = 'no_C2',group.by = "MC4Donors_Stim")

####motif analysis ####
DefaultAssay(export) = "chromvar"
Idents(export) = "MC4Donors_Stim"
MC4vsother_unstim_tf = FindMarkers(export,ident.1 = 'MC4Donor_C0',ident.2 = 'no_C0',only.pos = TRUE,mean.fxn = rowMeans,fc.name = "avg_diff")

MC4vsother_stim_tf = FindMarkers(export,ident.1 = 'MC4Donor_C2',ident.2 = 'no_C2',only.pos = TRUE,mean.fxn = rowMeans,fc.name = "avg_diff")

stimvsunstim_tf = FindMarkers(export,ident.1 = "C2",ident.2 = "C0",group.by = "Stimulation",only.pos = TRUE,mean.fxn  = rowMeans,fc.name = "avg_diff")

####pathway analysis#####
export_pb = AggregateExpression(export,assays = "RNA",group.by = "DonorID_TimePoint_Stim")
exportMeta  = export@meta.data
exportMeta = exportMeta[,c(26,30)] #### column numbers may differ a bit
exportMeta = exportMeta[!duplicated(exportMeta),]
export_pb = as.matrix(export_pb$RNA)
colnames(export_pb) = gsub("g","",colnames(export_pb))
colnames(export_pb) = gsub("-","_",colnames(export_pb))
export_pb = subset(export_pb,!grepl("^AL[0-6].+",rownames(export_pb)))
export_pb = subset(export_pb,!grepl("^AC[0-6].+",rownames(export_pb)))
export_pb = subset(export_pb,!grepl("^AP[0-6].+",rownames(export_pb)))
export_pb = subset(export_pb,!grepl("^MT-",rownames(export_pb)))
export_pb = subset(export_pb,!grepl("IGH.+",rownames(export_pb)))

exportMeta = exportMeta[match(colnames(export_pb),exportMeta$DonorID_TimePoint_Stim),]
dds <- DESeqDataSetFromMatrix(countData = export_pb,
                              colData = exportMeta,
                              design= ~ 0+MC4Donors_Stim)

smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]


dds <- DESeq(dds)


pb_MC4vsOther_C2 = results(dds, contrast=c("MC4Donors_Stim","MC4Donor_C2","no_C2"))
pb_MC4vsOther_C0 = results(dds, contrast=c("MC4Donors_Stim","MC4Donor_C0","no_C0"))

library(fgsea)
fgsea::gmtPathways('c2.cp.reactome.v7.5.1.symbols.gmt') -> genesets
fgsea::gmtPathways("h.all.v2022.1.Hs.symbols.gmt") -> genesets_hallmark


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

fgsea_pb_MC4vsOther_C2 = enrichmentFunc(pb_MC4vsOther_C2)
fgsea_pb_MC4vsOther_C0 = enrichmentFunc(pb_MC4vsOther_C0)


